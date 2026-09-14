# /// script
# requires-python = ">=3.10"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0", "numpy>=1.26"]
# ///
"""
nvidia-bio-mcp — servidor MCP para os NIM de biologia da NVIDIA (BioNeMo).

Por que este servidor existe
----------------------------
O BioNeMo Agent Toolkit da NVIDIA entrega *Skills* (markdown + scripts) para o
Claude Code. Skills sao otimas para procedimento, mas cada chamada vira um
`curl` reinventado pelo modelo: sem retry, sem polling assincrono, sem controle
de orcamento de tokens quando a resposta e um mmCIF de 2 MB.

Este servidor cobre essa lacuna. Ele expoe os NIM como ferramentas MCP tipadas,
trata o protocolo assincrono do NVCF (202 + nvcf-reqid + polling), e **grava
estruturas em disco em vez de despejar no contexto**, devolvendo so o caminho e
um resumo. Isso e a diferenca entre 300 tokens e 400 mil tokens por predicao.

Dois modos, alternados por variavel de ambiente:

  hospedado (padrao)   NVIDIA_NIM_MODE=hosted
                       https://health.api.nvidia.com/v1/biology/<org>/<model>/<action>
                       Auth: Bearer $NVIDIA_API_KEY
                       Sem VRAM local. Esta e a rota certa para GPU de 8 GB.

  local                NVIDIA_NIM_MODE=local
                       http://localhost:<porta>/biology/<org>/<model>/<action>
                       Containers NIM do nvcr.io. Exige >= 24 GB de VRAM para
                       os modelos de estrutura.

Uso academico: todas as ferramentas aqui sao de *analise e predicao* sobre
sequencias e estruturas — dobramento, docking, embeddings, busca de MSA. Cada
chamada e registrada em runs/provenance.jsonl para caderno de laboratorio.
"""

from __future__ import annotations

import asyncio
import base64
import hashlib
import io
import json
import math
import os
import re
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal

import httpx
import numpy as np

# O SDK Python do MCP renomeou FastMCP -> MCPServer na versao 2.0. Suportar as
# duas evita que a stack quebre dependendo de qual versao o uv resolver.
try:
    from mcp.server.mcpserver import MCPServer as _Server  # mcp >= 2.0
except ImportError:  # pragma: no cover - caminho legado
    from mcp.server.fastmcp import FastMCP as _Server  # mcp 1.x

# ---------------------------------------------------------------------------
# Configuracao
# ---------------------------------------------------------------------------

MODE: str = os.environ.get("NVIDIA_NIM_MODE", "hosted").lower()
API_KEY: str = os.environ.get("NVIDIA_API_KEY", "").strip()
HOSTED_BASE: str = os.environ.get("NVIDIA_NIM_HOSTED_BASE", "https://health.api.nvidia.com/v1")
LOCAL_BASE: str = os.environ.get("NVIDIA_NIM_LOCAL_BASE", "http://localhost:8000")

# Onde as estruturas/resultados pesados sao gravados.
OUT_DIR = Path(os.environ.get("BIO_STACK_OUT", "./runs/nim")).resolve()
PROV_LOG = Path(os.environ.get("BIO_STACK_PROVENANCE", "./runs/provenance.jsonl")).resolve()

# Predicoes de estrutura podem levar minutos. Nao seja impaciente.
REQUEST_TIMEOUT = float(os.environ.get("NVIDIA_NIM_TIMEOUT", "900"))
POLL_INTERVAL = 5.0
POLL_MAX_SECONDS = REQUEST_TIMEOUT

mcp = _Server("nvidia-bio")


# ---------------------------------------------------------------------------
# Registro de modelos
#
# Os slugs de org/modelo mudam conforme a NVIDIA promove modelos no catalogo.
# Em vez de fingir certeza, cada entrada carrega `confidence` e existe a
# ferramenta `nim_probe_endpoints` que verifica ao vivo quais respondem.
# Para qualquer NIM fora desta lista, use `nim_invoke` com o caminho cru.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NimModel:
    key: str
    org: str
    model: str
    action: str
    task: str
    summary: str
    confidence: Literal["confirmed", "likely"] = "likely"
    local_port: int = 8000
    notes: str = ""
    # Preenchido quando a NVIDIA desliga o endpoint (HTTP 410). O modelo fica
    # no registro para o erro explicar o que houve, mas nao e mais chamavel.
    retired: str = ""

    @property
    def path(self) -> str:
        # Varios NIM hospedados nao tem segmento de acao; barra final da 404.
        tail = f"/{self.action}" if self.action else ""
        return f"biology/{self.org}/{self.model}{tail}"


# Camada de embedding padrao do Evo 2. Nao ha "nenhuma camada": o /forward
# recusa output_layers vazio, entao omitir o parametro precisa de um default
# concreto. blocks.28.mlp.l3 e a camada intermediaria que a propria NVIDIA usa
# nos exemplos; devolve (1, n_bases, 8192) float64.
EVO2_DEFAULT_LAYER = "blocks.28.mlp.l3"


REGISTRY: dict[str, NimModel] = {
    m.key: m
    for m in [
        NimModel(
            key="boltz2",
            org="mit",
            model="boltz2",
            action="predict",
            task="cofolding",
            summary="Estrutura de complexo proteina-ligante + predicao de afinidade (IC50).",
            confidence="confirmed",
            notes="Aceita polymers (protein/dna/rna) e ligands (SMILES/CCD) na mesma chamada.",
        ),
        NimModel(
            key="openfold2",
            org="openfold",
            model="openfold2",
            action="predict-structure-from-msa-and-template",
            task="folding",
            summary="Dobramento de proteina a partir de MSA pre-computado. Precisa de msa_search antes.",
        ),
        NimModel(
            key="openfold3",
            org="openfold",
            model="openfold3",
            action="predict",
            task="folding",
            summary="Sucessor do OpenFold2, suporta complexos biomoleculares.",
        ),
        NimModel(
            key="esmfold",
            org="meta",
            model="esmfold",
            action="",
            task="folding",
            summary="Dobramento sem MSA, rapido. Menos preciso que AF2 em alvos sem homologos.",
            retired=(
                "endpoint hospedado desligado pela NVIDIA — HTTP 410 Gone em "
                "biology/meta/esmfold, verificado em 2026-09-11. Nao ha slug "
                "substituto; use boltz2, openfold2 ou alphafold2."
            ),
        ),
        NimModel(
            key="alphafold2",
            org="deepmind",
            model="alphafold2",
            action="",
            task="folding",
            summary="AlphaFold2 monomerico, pipeline completo incluindo MSA.",
        ),
        NimModel(
            key="alphafold2_multimer",
            org="deepmind",
            model="alphafold2-multimer",
            action="",
            task="folding",
            summary="AlphaFold2-Multimer: estrutura de complexo proteico (varias cadeias).",
            confidence="confirmed",
            notes=(
                "campo obrigatorio e 'sequences' (lista), nao 'sequence' como no "
                "alphafold2 monomerico. sem segmento de acao. rota verificada com "
                "chamada real em 2026-09-12 (422 pedindo 'sequences')."
            ),
        ),
        NimModel(
            key="msa_search",
            org="colabfold",
            model="msa-search",
            action="predict",
            task="msa",
            summary="Busca de alinhamento multiplo (UniRef30, ColabFoldDB). Entrada do OpenFold2.",
        ),
        NimModel(
            key="diffdock",
            org="mit",
            model="diffdock",
            action="",
            task="docking",
            summary="Docking cego proteina-ligante por difusao. Devolve poses ranqueadas.",
            notes=(
                "sem segmento de acao: biology/mit/diffdock/generate devolve 404, "
                "o caminho sem acao devolve 422/500. verificado com chamada real em 2026-09-12."
            ),
        ),
        NimModel(
            key="rfdiffusion",
            org="ipd",
            model="rfdiffusion",
            action="generate",
            task="design",
            summary="Gera backbone de proteina (de novo ou binder condicionado a um alvo).",
            confidence="confirmed",
            notes=(
                "Saida e so backbone (chave output_pdb); encadeie com proteinmpnn "
                "para obter sequencia. verificado com chamada real em 2026-09-12."
            ),
        ),
        NimModel(
            key="proteinmpnn",
            org="ipd",
            model="proteinmpnn",
            action="predict",
            task="design",
            summary="Desenha sequencia de aminoacidos compativel com um backbone dado.",
            confidence="confirmed",
            notes=(
                "Devolve mfasta (multi-FASTA), scores e probs "
                "[n_seq][n_residuo][21]. verificado com chamada real em 2026-09-12."
            ),
        ),
        NimModel(
            key="molmim",
            org="nvidia",
            model="molmim",
            action="generate",
            task="chemistry",
            summary="Geracao de moleculas pequenas otimizadas por propriedade (QED, logP).",
            confidence="confirmed",
            notes=(
                "Devolve 'molecules' como JSON dentro de uma string, e chama o SMILES "
                "de 'sample' em vez de 'smiles'. verificado com chamada real em 2026-09-12."
            ),
        ),
        NimModel(
            key="genmol",
            org="nvidia",
            model="genmol",
            action="generate",
            task="chemistry",
            summary="Geracao fragmento-a-molecula por difusao discreta mascarada.",
            confidence="confirmed",
            notes=(
                "Entrada precisa da mascara [*{N-N}] indicando onde crescer. "
                "verificado com chamada real em 2026-09-12."
            ),
        ),
        NimModel(
            key="evo2",
            org="arc",
            model="evo2-40b",
            action="forward",
            task="genomics",
            summary="Modelo de linguagem de DNA: embeddings de camadas internas.",
            confidence="confirmed",
            notes=(
                "/forward devolve SO as camadas pedidas em output_layers, nunca "
                "logits. output_layers e obrigatorio e precisa de pelo menos um "
                "item — lista vazia da 422. Para verossimilhanca use evo2_likelihood."
            ),
        ),
        NimModel(
            key="evo2_generate",
            org="arc",
            model="evo2-40b",
            action="generate",
            task="genomics",
            summary="Evo 2 em modo autoregressivo; com enable_sampled_probs devolve prob por token.",
            confidence="confirmed",
            notes=(
                "sampled_probs e o caminho de log-likelihood para comparar "
                "referencia contra variante. verificado com chamada real em 2026-09-12."
            ),
        ),
    ]
}


# ---------------------------------------------------------------------------
# Infraestrutura HTTP
# ---------------------------------------------------------------------------


class NimError(RuntimeError):
    """Falha previsivel de chamada ao NIM, com mensagem util para o modelo."""


def _base_url() -> str:
    return LOCAL_BASE if MODE == "local" else HOSTED_BASE


def _headers() -> dict[str, str]:
    h = {"Content-Type": "application/json", "Accept": "application/json"}
    if MODE == "hosted":
        if not API_KEY:
            raise NimError(
                "NVIDIA_API_KEY nao esta definida. Rode bootstrap/01_setup_keys.ps1 "
                "e reinicie o Claude Code."
            )
        h["Authorization"] = f"Bearer {API_KEY}"
    return h


def _url_for(path: str) -> str:
    return f"{_base_url().rstrip('/')}/{path.lstrip('/')}"


async def _poll_nvcf(client: httpx.AsyncClient, req_id: str) -> httpx.Response:
    """
    NVCF devolve 202 para trabalhos longos. O resultado sai por polling no
    endpoint de status ate virar 200.
    """
    status_url = f"{HOSTED_BASE.rstrip('/')}/status/{req_id}"
    deadline = time.monotonic() + POLL_MAX_SECONDS
    attempt = 0

    while time.monotonic() < deadline:
        await asyncio.sleep(POLL_INTERVAL)
        attempt += 1
        r = await client.get(status_url, headers=_headers(), timeout=60)
        if r.status_code == 200:
            return r
        if r.status_code != 202:
            raise NimError(
                f"polling falhou com HTTP {r.status_code} apos {attempt} tentativas: {r.text[:400]}"
            )

    raise NimError(
        f"tempo esgotado ({POLL_MAX_SECONDS:.0f}s) esperando o job {req_id}. "
        "Aumente NVIDIA_NIM_TIMEOUT ou reduza o tamanho da entrada."
    )


CACHE_DIR = OUT_DIR / "cache"
CACHE_ON = os.environ.get("BIO_STACK_NIM_CACHE", "1") not in ("0", "false", "off")


def _cache_key(path: str, payload: dict[str, Any]) -> str:
    """
    Endereca pelo conteudo: mesma rota + mesmos parametros = mesma chave.

    `sort_keys` importa mais do que parece — sem ele, dois dicts identicos com
    ordem de insercao diferente geram chaves diferentes e o cache nunca acerta.
    """
    canonico = json.dumps({"path": path, "payload": payload}, sort_keys=True, default=str)
    return hashlib.sha256(canonico.encode()).hexdigest()


def _cache_read(chave: str) -> dict[str, Any] | None:
    if not CACHE_ON:
        return None
    arq = CACHE_DIR / f"{chave}.json"
    try:
        envelope = json.loads(arq.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    corpo = envelope.get("response")
    return corpo if isinstance(corpo, dict) else None


def _cache_write(chave: str, path: str, resposta: dict[str, Any]) -> None:
    if not CACHE_ON:
        return
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        # O envelope guarda quando e de onde veio. Uma predicao reaproveitada
        # seis meses depois continua valida como numero, mas a banca pode
        # querer saber que ela nao foi recalculada — e uma busca de MSA contra
        # um banco que foi atualizado nesse meio tempo talvez devesse ser.
        (CACHE_DIR / f"{chave}.json").write_text(
            json.dumps(
                {
                    "gravado_em": datetime.now(UTC).isoformat(timespec="seconds"),
                    "rota": path,
                    "response": resposta,
                },
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
    except OSError:
        pass  # cache e otimizacao; falhar em gravar nunca pode derrubar a chamada


async def _post(path: str, payload: dict[str, Any]) -> dict[str, Any]:
    # Predicao de estrutura e a operacao cara desta stack. Dobrar a mesma
    # sequencia duas vezes custava duas vezes; numa varredura de mutantes com
    # prefixo compartilhado, a busca de MSA era recomprada a cada variante.
    chave = _cache_key(path, payload)
    guardado = _cache_read(chave)
    if guardado is not None:
        _log_provenance("_cache_hit", {"rota": path, "chave": chave[:12]}, {"reaproveitado": True})
        return guardado

    url = _url_for(path)
    async with httpx.AsyncClient(timeout=REQUEST_TIMEOUT, follow_redirects=True) as client:
        try:
            r = await client.post(url, headers=_headers(), json=payload)
        except httpx.ConnectError as exc:
            hint = (
                f"nao consegui conectar em {url}. Modo local: o container NIM esta rodando? "
                "Cheque com `docker ps`."
                if MODE == "local"
                else f"nao consegui conectar em {url}. Cheque a rede."
            )
            raise NimError(hint) from exc

        if r.status_code == 202:
            req_id = r.headers.get("nvcf-reqid")
            if not req_id:
                raise NimError("recebi 202 sem o header nvcf-reqid; nao da para fazer polling.")
            r = await _poll_nvcf(client, req_id)

        if r.status_code == 401:
            raise NimError("HTTP 401: NVIDIA_API_KEY invalida ou expirada.")
        if r.status_code == 403:
            raise NimError(
                "HTTP 403: chave valida mas sem permissao para este modelo. "
                "Verifique se sua licenca NVIDIA AI Enterprise cobre este NIM."
            )
        if r.status_code == 404:
            raise NimError(
                f"HTTP 404 em {url}. O slug do modelo provavelmente mudou. "
                "Rode nim_probe_endpoints para descobrir quais respondem."
            )
        if r.status_code == 429:
            raise NimError("HTTP 429: limite de taxa atingido. Espere e tente de novo.")
        if r.status_code >= 400:
            raise NimError(f"HTTP {r.status_code} em {url}: {r.text[:600]}")

        try:
            data = r.json()
        except json.JSONDecodeError as exc:
            raise NimError(f"resposta nao era JSON: {r.text[:400]}") from exc
        if not isinstance(data, dict):
            raise NimError(
                f"resposta JSON nao era objeto (tipo {type(data).__name__}): {r.text[:300]}"
            )
        # So sucesso entra no cache: todo caminho de erro acima levanta NimError
        # e nunca chega aqui. Guardar um 429 seria transformar uma falha
        # temporaria em falha permanente.
        _cache_write(chave, path, data)
        return data


# ---------------------------------------------------------------------------
# Persistencia e proveniencia
# ---------------------------------------------------------------------------


def _slug(text: str, n: int = 12) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:n]


def _save(content: str, stem: str, ext: str) -> Path:
    """
    Grava em OUT_DIR. `stem` vem de parametro de ferramenta e `ext` pode vir da
    resposta da API, entao os dois sao sanitizados: sem isso, um `name` como
    "../../.ssh/authorized_keys" escreveria fora do diretorio.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^A-Za-z0-9_.-]", "_", str(stem)).strip("._")[:64] or "output"
    ext = re.sub(r"[^A-Za-z0-9]", "", str(ext))[:8] or "txt"

    ts = datetime.now(UTC).strftime("%Y%m%dT%H%M%S_%f")
    path = (OUT_DIR / f"{ts}_{stem}.{ext}").resolve()
    if OUT_DIR.resolve() not in path.parents:
        raise NimError(f"caminho de saida escaparia de {OUT_DIR}: {path}")

    path.write_text(content, encoding="utf-8")
    return path


def _log_provenance(tool: str, params: dict[str, Any], result: dict[str, Any]) -> None:
    """Registro append-only para reprodutibilidade academica."""
    try:
        PROV_LOG.parent.mkdir(parents=True, exist_ok=True)
        safe = {k: v for k, v in params.items() if "key" not in k.lower()}
        entry = {
            "ts": datetime.now(UTC).isoformat(),
            "server": "nvidia-bio-mcp",
            "tool": tool,
            "mode": MODE,
            "endpoint": _base_url(),
            "params": safe,
            "outputs": result.get("files", []),
        }
        with PROV_LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except OSError:
        pass  # log nunca deve derrubar a ferramenta


AA = set("ACDEFGHIKLMNPQRSTVWY")
NT = set("ACGTUN")


def _validate_protein(seq: str) -> str:
    s = re.sub(r"\s+", "", seq).upper()
    if not s:
        raise NimError("sequencia vazia.")
    bad = set(s) - AA - {"X"}
    if bad:
        raise NimError(
            f"caracteres invalidos para proteina: {sorted(bad)}. "
            "Esperado apenas o alfabeto de 20 aminoacidos (X permitido)."
        )
    if len(s) > 4000:
        raise NimError(f"sequencia com {len(s)} residuos excede o limite pratico de 4000.")
    return s


def _validate_nucleotide(seq: str) -> str:
    s = re.sub(r"\s+", "", seq).upper()
    if not s:
        raise NimError("sequencia vazia.")
    bad = set(s) - NT
    if bad:
        raise NimError(f"caracteres invalidos para acido nucleico: {sorted(bad)}")
    return s


# ---------------------------------------------------------------------------
# Ferramentas: diagnostico
# ---------------------------------------------------------------------------


@mcp.tool()
def nim_config() -> dict[str, Any]:
    """
    Mostra como este servidor esta configurado: modo, endpoint base, se a chave
    esta presente, e onde os resultados sao gravados.

    Chame isto primeiro quando qualquer outra ferramenta falhar.
    """
    return {
        "mode": MODE,
        "base_url": _base_url(),
        "api_key_present": bool(API_KEY),
        "api_key_prefix": (API_KEY[:8] + "...") if API_KEY else None,
        "output_dir": str(OUT_DIR),
        "provenance_log": str(PROV_LOG),
        "timeout_seconds": REQUEST_TIMEOUT,
        "models_registered": sorted(REGISTRY),
        "cache": _cache_stats(),
    }


def _cache_stats() -> dict[str, Any]:
    """
    Quantas respostas ja estao guardadas e quanto ocupam.

    Serve para responder "esta varredura vai custar chamada nova?" antes de
    dispara-la, e para decidir quando limpar. Predicao nao expira sozinha:
    o modelo e deterministico, entao a resposta de ontem vale hoje.
    """
    if not CACHE_ON:
        return {"ativo": False}
    try:
        arquivos = list(CACHE_DIR.glob("*.json"))
    except OSError:
        return {"ativo": True, "entradas": 0, "bytes": 0, "dir": str(CACHE_DIR)}
    return {
        "ativo": True,
        "entradas": len(arquivos),
        "bytes": sum(a.stat().st_size for a in arquivos),
        "dir": str(CACHE_DIR),
        "desligar_com": "BIO_STACK_NIM_CACHE=0",
    }


@mcp.tool()
def nim_list_models(task: str | None = None) -> list[dict[str, Any]]:
    """
    Lista os NIM de biologia conhecidos por este servidor.

    Args:
        task: filtra por categoria — folding, cofolding, docking, design,
              chemistry, genomics, msa. Omita para listar tudo.
    """
    out = []
    for m in REGISTRY.values():
        if task and m.task != task:
            continue
        out.append(
            {
                "key": m.key,
                "task": m.task,
                "path": m.path,
                "summary": m.summary,
                "slug_confidence": m.confidence,
                "notes": m.notes or None,
            }
        )
    return sorted(out, key=lambda d: (d["task"], d["key"]))


@mcp.tool()
async def nim_probe_endpoints(keys: list[str] | None = None) -> dict[str, Any]:
    """
    Verifica ao vivo quais endpoints respondem com a chave atual.

    Faz um POST com corpo vazio em cada rota: 4xx de validacao significa que a
    rota EXISTE (so rejeitou o payload), enquanto 404 significa slug errado e
    401/403 significa problema de credencial. Use isto quando um modelo der 404.

    Args:
        keys: subconjunto de modelos a testar. Omita para testar todos.
    """
    targets = [REGISTRY[k] for k in (keys or REGISTRY.keys()) if k in REGISTRY]
    results: dict[str, str] = {}

    async with httpx.AsyncClient(timeout=45) as client:

        async def probe(m: NimModel) -> tuple[str, str]:
            try:
                r = await client.post(_url_for(m.path), headers=_headers(), json={})
            except httpx.HTTPError as exc:
                return m.key, f"rede: {type(exc).__name__}"
            code = r.status_code
            if code == 404:
                return m.key, "404 rota inexistente — slug provavelmente mudou"
            if code == 410:
                return m.key, "410 endpoint aposentado — removido do catalogo, sem slug novo"
            if code in (401, 403):
                return m.key, f"{code} credencial sem acesso a este modelo"
            if code in (200, 202):
                return m.key, f"{code} OK (aceitou ate corpo vazio)"
            if 400 <= code < 500:
                return m.key, f"{code} rota existe (erro de validacao, esperado)"
            if code >= 500:
                # Varios NIM devolvem 500 (em vez de 422) a corpo vazio. Isso nao
                # prova rota quebrada: molmim e diffdock respondem 200 com payload
                # real. Inconclusivo, portanto — nunca "inalcancavel".
                return m.key, f"{code} inconclusivo — teste com payload real"
            return m.key, f"{code} inesperado"

        for key, verdict in await asyncio.gather(*(probe(m) for m in targets)):
            results[key] = verdict

    reachable = {k for k, v in results.items() if "rota existe" in v or "OK" in v}
    inconclusive = {k for k, v in results.items() if "inconclusivo" in v}
    return {
        "mode": MODE,
        "base_url": _base_url(),
        "results": results,
        "reachable": sorted(reachable),
        "inconclusive": sorted(inconclusive),
        "unreachable": sorted(set(results) - reachable - inconclusive),
    }


# ---------------------------------------------------------------------------
# Ferramentas: estrutura
# ---------------------------------------------------------------------------


@mcp.tool()
async def fold_protein(
    sequence: str,
    model: str = "boltz2",
    name: str = "target",
    recycling_steps: int = 3,
    sampling_steps: int = 50,
    diffusion_samples: int = 1,
) -> dict[str, Any]:
    """
    Prediz a estrutura 3D de uma proteina a partir da sequencia de aminoacidos.

    A estrutura e gravada em disco (mmCIF) e o retorno traz apenas o caminho e
    metricas — coordenadas atomicas nunca entram no contexto.

    Args:
        sequence: sequencia de aminoacidos, alfabeto de 20 letras.
        model: boltz2 (padrao, melhor geral) | openfold2 | openfold3 |
               alphafold2. (esmfold foi aposentado pela NVIDIA.)
        name: rotulo curto usado no nome do arquivo.
        recycling_steps: passos de reciclagem; mais = melhor e mais lento.
        sampling_steps: passos de difusao.
        diffusion_samples: quantas estruturas candidatas gerar.
    """
    seq = _validate_protein(sequence)
    m = REGISTRY.get(model)
    if m is None or m.task not in ("folding", "cofolding"):
        raise NimError(
            f"modelo '{model}' nao serve para dobramento. "
            f"Opcoes: {[k for k, v in REGISTRY.items() if v.task in ('folding', 'cofolding')]}"
        )
    if m.retired:
        tarefas = ("folding", "cofolding")
        vivos = [k for k, v in REGISTRY.items() if v.task in tarefas and not v.retired]
        raise NimError(f"modelo '{model}' aposentado: {m.retired} Alternativas: {vivos}")

    # Cada NIM tem um schema proprio - nao existe payload universal.
    if model == "boltz2":
        payload: dict[str, Any] = {
            "polymers": [{"id": "A", "molecule_type": "protein", "sequence": seq}],
            "recycling_steps": recycling_steps,
            "sampling_steps": sampling_steps,
            "diffusion_samples": diffusion_samples,
            "step_scale": 1.638,
            "output_format": "mmcif",
        }
    elif model == "openfold3":
        payload = {
            "inputs": [
                {
                    "molecules": [{"type": "protein", "sequence": seq, "id": "A"}],
                    "diffusion_samples": diffusion_samples,
                    "output_format": "cif",
                }
            ]
        }
    else:
        # esmfold, alphafold2, openfold2
        payload = {"sequence": seq}

    data = await _post(m.path, payload)
    return _handle_structure_response(
        data,
        stem=f"{name}_{model}_{_slug(seq)}",
        tool="fold_protein",
        params={"model": model, "length": len(seq), "name": name},
    )


@mcp.tool()
async def cofold_complex(
    protein_sequence: str,
    ligand_smiles: str | None = None,
    partner_sequence: str | None = None,
    name: str = "complex",
    predict_affinity: bool = True,
) -> dict[str, Any]:
    """
    Prediz a estrutura de um complexo com Boltz-2: proteina + ligante pequeno,
    ou proteina + proteina. Opcionalmente estima afinidade de ligacao.

    Uso tipico academico: avaliar se um composto conhecido se encaixa no sitio
    ativo de uma enzima que voce esta estudando.

    Args:
        protein_sequence: cadeia principal (aminoacidos).
        ligand_smiles: SMILES do ligante, se for complexo com molecula pequena.
        partner_sequence: segunda cadeia proteica, se for complexo proteina-proteina.
        name: rotulo para o arquivo de saida.
        predict_affinity: pede tambem a estimativa de afinidade do Boltz-2.
    """
    seq = _validate_protein(protein_sequence)
    polymers: list[dict[str, Any]] = [{"id": "A", "molecule_type": "protein", "sequence": seq}]
    if partner_sequence:
        polymers.append(
            {"id": "B", "molecule_type": "protein", "sequence": _validate_protein(partner_sequence)}
        )

    payload: dict[str, Any] = {
        "polymers": polymers,
        "recycling_steps": 3,
        "sampling_steps": 50,
        "diffusion_samples": 1,
        "step_scale": 1.638,
        "output_format": "mmcif",
    }

    if ligand_smiles:
        # A flag mora dentro do ligante. No nivel superior ela e ignorada em
        # silencio - voce recebe a estrutura e nenhuma afinidade.
        payload["ligands"] = [
            {
                "id": "L",
                "smiles": ligand_smiles,
                "predict_affinity": predict_affinity,
            }
        ]

    if not ligand_smiles and not partner_sequence:
        raise NimError("informe ligand_smiles ou partner_sequence — senao use fold_protein.")

    data = await _post(REGISTRY["boltz2"].path, payload)
    return _handle_structure_response(
        data,
        stem=f"{name}_boltz2_{_slug(seq + (ligand_smiles or '') + (partner_sequence or ''))}",
        tool="cofold_complex",
        params={
            "protein_len": len(seq),
            "ligand_smiles": ligand_smiles,
            "has_partner": bool(partner_sequence),
            "name": name,
        },
    )


@mcp.tool()
async def msa_search(
    sequence: str,
    databases: list[str] | None = None,
    max_sequences: int = 1000,
    name: str = "query",
) -> dict[str, Any]:
    """
    Busca alinhamento multiplo de sequencias (MSA) — passo obrigatorio antes de
    OpenFold2 e que melhora muito a qualidade de qualquer predicao de estrutura.

    Args:
        sequence: sequencia de aminoacidos da consulta.
        databases: bases a consultar. Padrao ["Uniref30_2302"].
        max_sequences: teto de sequencias homologas a retornar.
        name: rotulo para o arquivo .a3m de saida.
    """
    seq = _validate_protein(sequence)
    payload = {
        "sequence": seq,
        "databases": databases or ["Uniref30_2302"],
        "e_value": 0.0001,
        "iterations": 1,
        "max_msa_sequences": max_sequences,
        "output_alignment_formats": ["a3m"],
    }
    data = await _post(REGISTRY["msa_search"].path, payload)

    files: list[str] = []
    n_hits = 0
    alignments = data.get("alignments")

    if isinstance(alignments, dict):
        for db, fmts in alignments.items():
            if not isinstance(fmts, dict):
                continue
            for fmt, obj in fmts.items():
                if fmt not in ("a3m", "sto", "fasta"):
                    continue  # nome de arquivo nao deve vir livre da resposta remota
                body = obj.get("alignment") if isinstance(obj, dict) else obj
                if isinstance(body, str) and body.strip():
                    safe_db = re.sub(r"[^A-Za-z0-9_.-]", "_", str(db))[:32]
                    files.append(str(_save(body, f"{name}_{safe_db}_{_slug(seq)}", fmt)))
                    # cada cabecalho '>' no a3m/fasta e uma sequencia homologa
                    n_hits = max(n_hits, body.count("\n>") + body.startswith(">"))

    if not files:
        files.append(str(_save(json.dumps(data, indent=2), f"{name}_msa_raw", "json")))

    result = {
        "files": files,
        "homologs_found": n_hits,
        "databases": payload["databases"],
        "next_step": "passe o caminho do .a3m para o OpenFold2 via nim_invoke, "
        "ou use fold_protein com boltz2 que ja lida com isso internamente.",
    }
    _log_provenance("msa_search", {"length": len(seq), "databases": payload["databases"]}, result)
    return result


@mcp.tool()
async def dock_ligand(
    protein_pdb_path: str,
    ligand_smiles: str,
    num_poses: int = 10,
    steps: int = 18,
    name: str = "docking",
) -> dict[str, Any]:
    """
    Docking cego proteina-ligante com DiffDock. Devolve poses ranqueadas por
    confianca, gravadas em disco.

    Args:
        protein_pdb_path: caminho de um arquivo .pdb local (receptor).
        ligand_smiles: SMILES do ligante.
        num_poses: quantas poses gerar.
        steps: passos de difusao inversa.
        name: rotulo para os arquivos de saida.
    """
    pdb_path = Path(protein_pdb_path).expanduser().resolve()
    if not pdb_path.is_file():
        raise NimError(f"arquivo PDB nao encontrado: {pdb_path}")
    pdb_text = pdb_path.read_text(encoding="utf-8")
    if "ATOM" not in pdb_text:
        raise NimError(f"{pdb_path} nao parece um PDB valido (sem registros ATOM).")

    payload = {
        "protein": pdb_text,
        "ligand": ligand_smiles,
        # A API aceita so mol2/sdf/txt. SMILES vai como texto puro: "txt".
        "ligand_file_type": "txt",
        "num_poses": num_poses,
        "time_divisions": 20,
        "steps": steps,
        "save_trajectory": False,
        "is_staged": False,
    }
    data = await _post(REGISTRY["diffdock"].path, payload)

    poses = data.get("ligand_positions") or []
    raw_scores = list(data.get("position_confidence") or [])
    # Parear ANTES de filtrar: se alguma pose nao for string, fatiar os scores
    # por comprimento os associaria as poses erradas.
    raw_scores += [None] * (len(poses) - len(raw_scores))
    pairs = [(p, sc) for p, sc in zip(poses, raw_scores, strict=False) if isinstance(p, str)]
    files = [str(_save(p, f"{name}_pose{i:02d}", "sdf")) for i, (p, _) in enumerate(pairs)]
    scores = [sc for _, sc in pairs if sc is not None]

    result = {
        "files": files,
        "num_poses": len(files),
        "confidence_scores": scores,
        "best_score": max(scores) if scores else None,
        "interpretation": "confianca DiffDock: >0 boa, -1.5 a 0 incerta, <-1.5 provavelmente ruim.",
    }
    _log_provenance("dock_ligand", {"receptor": str(pdb_path), "smiles": ligand_smiles}, result)
    return result


# ---------------------------------------------------------------------------
# Ferramentas: genomica
# ---------------------------------------------------------------------------


# O Evo 2 e byte-level: o vocabulario de 512 posicoes e indexado por ord() do
# caractere, entao o logit da base 'A' fica no indice 65. Medido em 2026-09-12:
# a massa de probabilidade em {A,C,G,T} soma 1,0, o que confirma a indexacao.
EVO2_LOGIT_LAYER = "unembed"


async def _evo2_loglik(seq: str) -> tuple[float, list[float]]:
    """
    Log-likelihood media por base de `seq` sob o Evo 2, em UMA chamada.

    A camada `unembed` do /forward devolve os logits por posicao da entrada —
    shape (1, len(seq), 512). Com log-softmax e um gather na base realmente
    observada na posicao seguinte, sai log P(base_i+1 | prefixo).

    Por que nao /generate: la o `sampled_probs` cobre os tokens GERADOS, nao a
    sequencia de entrada. Usar aquilo pontuaria a continuacao inventada pelo
    modelo, nao a sequencia que se quer avaliar.
    """
    data = await _post(
        REGISTRY["evo2"].path, {"sequence": seq, "output_layers": [EVO2_LOGIT_LAYER]}
    )

    blob = data.get("data")
    if not isinstance(blob, str):
        raise NimError(
            "resposta do Evo 2 sem o campo 'data' base64 — contrato do endpoint mudou. "
            "Rode nim_probe_endpoints."
        )

    npz = np.load(io.BytesIO(base64.b64decode(blob)))
    key = f"{EVO2_LOGIT_LAYER}.output"
    if key not in npz.files:
        raise NimError(f"NPZ sem a chave {key}; veio {npz.files}.")

    logits = np.asarray(npz[key][0], dtype=np.float64)  # (len(seq), vocab)
    if logits.shape[0] != len(seq):
        raise NimError(f"logits com {logits.shape[0]} posicoes para {len(seq)} bases.")

    shifted = logits - logits.max(-1, keepdims=True)
    logp = shifted - np.log(np.exp(shifted).sum(-1, keepdims=True))

    ids = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).astype(np.int64)
    # Posicao i preve a base i+1, entao a ultima base nao tem alvo.
    per_base = logp[np.arange(len(seq) - 1), ids[1:]]
    return float(per_base.mean()), [round(float(x), 5) for x in per_base]


@mcp.tool()
async def evo2_score(
    sequence: str,
    output_layers: list[str] | None = None,
    name: str = "dna",
) -> dict[str, Any]:
    """
    Passa uma sequencia de DNA pelo Evo 2 e devolve embeddings/logits.

    Uso academico central: comparar o log-likelihood de uma sequencia de
    referencia contra uma variante. Queda grande de verossimilhanca sugere que
    a variante e disruptiva — util para priorizar variantes de significado
    incerto em trabalho de genetica.

    Args:
        sequence: sequencia de DNA (A/C/G/T/N).
        output_layers: camadas internas a extrair. Omita para a camada padrao.
        name: rotulo do arquivo de saida.
    """
    seq = _validate_nucleotide(sequence)
    if len(seq) > 8192:
        raise NimError(
            f"sequencia de {len(seq)} bp excede a janela pratica de 8192 bp por chamada. "
            "Quebre em janelas deslizantes."
        )

    # output_layers e obrigatorio E precisa de pelo menos um item: lista vazia
    # devolve 422 "List should have at least 1 item after validation, not 0".
    # O comentario anterior aqui dizia que lista vazia significava "so logits";
    # esta errado — /forward nunca devolve logits, so as camadas pedidas.
    # Medido com chamada real em 2026-09-12.
    payload: dict[str, Any] = {
        "sequence": seq,
        "output_layers": output_layers or [EVO2_DEFAULT_LAYER],
    }

    data = await _post(REGISTRY["evo2"].path, payload)

    # Embeddings sao enormes: gravar, nao devolver.
    path = _save(json.dumps(data), f"{name}_evo2_{_slug(seq)}", "json")
    result = {
        "files": [str(path)],
        "sequence_length": len(seq),
        "response_keys": sorted(data.keys()),
        "note": "tensores gravados em disco. Carregue com json/numpy para analise.",
    }
    _log_provenance("evo2_score", {"length": len(seq), "name": name}, result)
    return result


@mcp.tool()
async def evo2_likelihood(
    sequence: str,
    variant_sequence: str | None = None,
    name: str = "dna",
) -> dict[str, Any]:
    """
    Log-likelihood do Evo 2 para uma sequencia de DNA, e a comparacao
    referencia vs. variante — o uso academico central do modelo.

    Queda de verossimilhanca na variante sugere que ela e disruptiva. Serve
    para priorizar variantes de significado incerto antes de gastar tempo com
    curadoria manual.

    Nao confunda com evo2_score, que extrai embeddings de camada interna e nao
    devolve probabilidade nenhuma.

    Args:
        sequence: sequencia de referencia (A/C/G/T/N).
        variant_sequence: a variante. Omita para pontuar so a referencia.
        name: rotulo do arquivo de saida.
    """
    ref = _validate_nucleotide(sequence)
    if len(ref) > 8192:
        raise NimError(
            f"sequencia de {len(ref)} bp excede a janela pratica de 8192 bp por chamada. "
            "Quebre em janelas deslizantes."
        )

    ref_ll, ref_per_base = await _evo2_loglik(ref)
    result: dict[str, Any] = {
        "reference": {
            "length": len(ref),
            "mean_log_likelihood": round(ref_ll, 6),
            "perplexity": round(math.exp(-ref_ll), 4),
        }
    }
    var_per_base: list[float] = []

    if variant_sequence:
        var = _validate_nucleotide(variant_sequence)
        if len(var) > 8192:
            raise NimError(f"variante de {len(var)} bp excede 8192 bp.")
        var_ll, var_per_base = await _evo2_loglik(var)
        result["variant"] = {
            "length": len(var),
            "mean_log_likelihood": round(var_ll, 6),
            "perplexity": round(math.exp(-var_ll), 4),
        }
        result["delta_log_likelihood"] = round(var_ll - ref_ll, 6)
        result["reading"] = (
            "delta negativo = variante menos provavel que a referencia sob o "
            "modelo, ou seja, candidata a disruptiva. E sinal de priorizacao, "
            "nao criterio ACMG — confirme com variant-curator, frequencia "
            "populacional e literatura."
        )
        if len(ref) != len(var):
            result["aviso"] = (
                "comprimentos diferentes: a media por base ja normaliza, mas "
                "indel muda o registro de leitura e o delta fica menos interpretavel "
                "que numa troca de base unica."
            )

    path = _save(
        json.dumps({"ref_per_base": ref_per_base, "var_per_base": var_per_base}),
        f"{name}_evo2ll_{_slug(ref)}",
        "json",
    )
    result["files"] = [str(path)]
    result["note"] = (
        "log-likelihood por base gravado em disco; aqui volta so o resumo. "
        "O indice i do vetor e a predicao da base i+1 (1-based: posicao i+2), "
        "porque a posicao i so ve o prefixo ate ela. A maior queda entre ref e "
        "variante costuma cair no sitio trocado."
    )
    _log_provenance(
        "evo2_likelihood",
        {"ref_length": len(ref), "has_variant": bool(variant_sequence), "name": name},
        result,
    )
    return result


# ---------------------------------------------------------------------------
# Ferramentas: desenho de proteina
#
# RFdiffusion gera o esqueleto; ProteinMPNN veste esse esqueleto com uma
# sequencia. Os dois em sequencia sao o fluxo canonico de desenho in silico, e
# e por isso que aparecem em par aqui.
# ---------------------------------------------------------------------------


def _read_pdb(path_str: str) -> tuple[Path, str]:
    """Le um PDB local. O texto vai para o payload, nunca para o contexto."""
    p = Path(path_str).expanduser().resolve()
    if not p.is_file():
        raise NimError(f"arquivo PDB nao encontrado: {p}")
    text = p.read_text(encoding="utf-8")
    if "ATOM" not in text:
        raise NimError(f"{p} nao parece um PDB valido (sem registros ATOM).")
    return p, text


def _parse_mfasta(text: str) -> list[dict[str, Any]]:
    """
    Quebra o multi-FASTA do ProteinMPNN.

    O primeiro registro e sempre a sequencia nativa do PDB de entrada
    (cabecalho ">input"); os seguintes sao os desenhos, cujo cabecalho traz
    score e seq_recovery.
    """
    records: list[dict[str, Any]] = []
    header: str | None = None
    seq: list[str] = []

    def flush() -> None:
        if header is None:
            return
        fields: dict[str, float] = {}
        for part in header.split(","):
            k, _, v = part.partition("=")
            k, v = k.strip(), v.strip()
            if re.fullmatch(r"-?\d+(\.\d+)?", v):
                fields[k] = float(v)
        records.append(
            {
                "is_native": header.startswith("input"),
                "sequence": "".join(seq),
                "score": fields.get("score"),
                "seq_recovery": fields.get("seq_recovery"),
            }
        )

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header, seq = line[1:].strip(), []
        else:
            seq.append(line)
    flush()
    return records


def _handle_mpnn_response(
    data: dict[str, Any], stem: str, params: dict[str, Any]
) -> dict[str, Any]:
    """Normaliza a resposta do ProteinMPNN: mfasta + scores + probs."""
    mfasta = data.get("mfasta")
    if not isinstance(mfasta, str) or not mfasta.strip():
        raise NimError(
            f"resposta do ProteinMPNN sem 'mfasta'. Chaves recebidas: {sorted(data.keys())}"
        )

    files = [str(_save(mfasta, stem, "fasta"))]

    # probs e [n_seq][n_residuo][21] — dezenas de milhares de floats mesmo para
    # uma proteina pequena. Vai para disco, sempre.
    probs = data.get("probs")
    if probs:
        files.append(str(_save(json.dumps(probs), f"{stem}_probs", "json")))

    records = _parse_mfasta(mfasta)
    native = next((r for r in records if r["is_native"]), None)
    designs = [r for r in records if not r["is_native"]]

    result: dict[str, Any] = {
        "files": files,
        "num_designs": len(designs),
        "native_sequence": (native or {}).get("sequence"),
        "designs": [
            {"sequence": d["sequence"], "score": d["score"], "seq_recovery": d["seq_recovery"]}
            for d in designs
        ],
        "interpretation": (
            "score = perplexidade media do desenho, menor e melhor. "
            "seq_recovery = fracao de residuos identicos a nativa; 0.4-0.5 e "
            "tipico e saudavel, perto de 1.0 indica que o modelo apenas copiou "
            "a entrada em vez de desenhar."
        ),
    }
    if probs:
        result["probs_note"] = (
            "matriz de probabilidade por residuo (n_seq x n_residuo x 21) gravada em "
            "disco; a ultima dimensao e o alfabeto ACDEFGHIKLMNPQRSTVWYX."
        )

    _log_provenance("design_sequence", params, result)
    return result


@mcp.tool()
async def design_backbone(
    input_pdb_path: str,
    contigs: str,
    hotspot_res: list[str] | None = None,
    diffusion_steps: int = 15,
    name: str = "backbone",
) -> dict[str, Any]:
    """
    RFdiffusion: gera um backbone de proteina condicionado a uma estrutura de
    referencia. A saida e so o esqueleto, sem sequencia — encadeie com
    `design_sequence` para obter os aminoacidos.

    Uso academico tipico: desenhar um scaffold ou um binder contra um alvo
    modelo (lisozima, GFP, crambina) como exercicio de desenho estrutural.

    Args:
        input_pdb_path: caminho do .pdb de referencia (o alvo).
        contigs: especificacao RFdiffusion do que preservar e do que gerar.
                 "A1-46/0 20-20" preserva os residuos 1-46 da cadeia A e
                 desenha 20 residuos novos numa cadeia separada;
                 "A1-46" sozinho faz difusao parcial sobre a cadeia inteira.
        hotspot_res: residuos do alvo que o binder deve tocar, ex. ["A23","A27"].
        diffusion_steps: passos de difusao; mais = melhor e mais lento.
        name: rotulo curto usado no nome do arquivo.
    """
    pdb_path, pdb_text = _read_pdb(input_pdb_path)
    contigs = str(contigs).strip()
    if not contigs:
        raise NimError(
            "contigs e obrigatorio. Ex.: 'A1-46/0 20-20' preserva A1-46 e desenha "
            "20 residuos novos."
        )

    payload: dict[str, Any] = {
        "input_pdb": pdb_text,
        "contigs": contigs,
        "hotspot_res": list(hotspot_res or []),
        "diffusion_steps": int(diffusion_steps),
    }
    data = await _post(REGISTRY["rfdiffusion"].path, payload)

    result = _handle_structure_response(
        data,
        stem=f"{name}_rfdiff_{_slug(contigs)}",
        tool="design_backbone",
        params={
            "receptor": str(pdb_path),
            "contigs": contigs,
            "hotspot_res": list(hotspot_res or []),
            "diffusion_steps": diffusion_steps,
        },
    )
    result["next_step"] = (
        "backbone sem sequencia. Rode design_sequence neste arquivo para obter os aminoacidos."
    )
    return result


@mcp.tool()
async def design_sequence(
    input_pdb_path: str,
    chains: list[str] | None = None,
    num_sequences: int = 4,
    sampling_temp: float = 0.1,
    use_soluble_model: bool = False,
    ca_only: bool = False,
    name: str = "design",
) -> dict[str, Any]:
    """
    ProteinMPNN: desenha sequencias de aminoacidos compativeis com um backbone.

    Serve tanto para um backbone recem-gerado pelo `design_backbone` quanto para
    redesenhar uma proteina experimental existente. O retorno traz as sequencias
    e o seq_recovery; a matriz de probabilidade por residuo vai para disco.

    Args:
        input_pdb_path: caminho do .pdb com o backbone.
        chains: cadeias a desenhar, ex. ["A"]. Lista vazia = todas.
        num_sequences: quantas sequencias amostrar.
        sampling_temp: temperatura em (0, 1]. 0.1 conservador (recupera muito da
                       nativa), 0.3 diversifica, acima disso degrada.
        use_soluble_model: pesos treinados so em proteinas soluveis. Ligue quando
                           o alvo nao for de membrana.
        ca_only: usar apenas os carbonos alfa do backbone.
        name: rotulo curto usado no nome do arquivo.
    """
    pdb_path, pdb_text = _read_pdb(input_pdb_path)
    if not 0.0 < sampling_temp <= 1.0:
        raise NimError(
            f"sampling_temp={sampling_temp} fora de (0, 1]. Use 0.1 para conservador "
            "ou 0.3 para diversificar."
        )
    if num_sequences < 1:
        raise NimError("num_sequences deve ser >= 1.")

    payload: dict[str, Any] = {
        "input_pdb": pdb_text,
        "input_pdb_chains": list(chains or []),
        "ca_only": bool(ca_only),
        "use_soluble_model": bool(use_soluble_model),
        "num_seq_per_target": int(num_sequences),
        "sampling_temp": [float(sampling_temp)],
    }
    data = await _post(REGISTRY["proteinmpnn"].path, payload)
    return _handle_mpnn_response(
        data,
        stem=f"{name}_mpnn_{_slug(pdb_text)}",
        params={
            "backbone": str(pdb_path),
            "chains": list(chains or []),
            "num_sequences": num_sequences,
            "sampling_temp": sampling_temp,
            "use_soluble_model": use_soluble_model,
        },
    )


# ---------------------------------------------------------------------------
# Ferramentas: quimica generativa
# ---------------------------------------------------------------------------


def _normalize_molecules(data: dict[str, Any]) -> list[dict[str, Any]]:
    """
    MolMIM e GenMol devolvem a mesma ideia em formatos diferentes: o MolMIM manda
    `molecules` como JSON *dentro de uma string* e chama o SMILES de "sample"; o
    GenMol manda lista de verdade e chama de "smiles".
    """
    raw = data.get("molecules")
    if isinstance(raw, str):
        try:
            raw = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise NimError(f"campo 'molecules' nao era JSON valido: {raw[:200]}") from exc
    if not isinstance(raw, list):
        raise NimError(f"resposta sem lista 'molecules'. Chaves recebidas: {sorted(data.keys())}")

    out: list[dict[str, Any]] = []
    for item in raw:
        if isinstance(item, str) and item:
            out.append({"smiles": item, "score": None})
        elif isinstance(item, dict):
            smi = item.get("smiles") or item.get("sample")
            if isinstance(smi, str) and smi:
                out.append({"smiles": smi, "score": item.get("score")})
    return out


def _handle_molecule_response(
    data: dict[str, Any], stem: str, tool: str, params: dict[str, Any]
) -> dict[str, Any]:
    mols = _normalize_molecules(data)
    path = _save(json.dumps(mols, indent=2), stem, "json")

    unique = {m["smiles"] for m in mols}
    result: dict[str, Any] = {
        "files": [str(path)],
        "num_generated": len(mols),
        "num_unique": len(unique),
        "score_type": data.get("score_type") or params.get("scoring"),
        # SMILES sao curtos: cabem inline, e e isso que se quer olhar.
        "molecules": mols[:50],
    }
    if len(mols) > 50:
        result["note"] = f"mostrando 50 de {len(mols)}; lista completa em {path}."
    if len(unique) < len(mols):
        result["duplicates_note"] = (
            "ha SMILES repetidos — normal quando o otimizador converge. "
            "Conte por num_unique, nao por num_generated."
        )

    _log_provenance(tool, params, result)
    return result


@mcp.tool()
async def generate_molecules(
    smiles: str,
    property_name: str = "QED",
    num_molecules: int = 10,
    minimize: bool = False,
    min_similarity: float = 0.7,
    particles: int = 30,
    iterations: int = 10,
    name: str = "molmim",
) -> dict[str, Any]:
    """
    MolMIM: parte de uma molecula semente e busca analogos que melhorem uma
    propriedade, mantendo semelhanca quimica com a semente (CMA-ES no espaco
    latente do modelo).

    Uso academico tipico: exercicio de otimizacao de lead — pegar um composto
    conhecido e ver que vizinhanca quimica o modelo propoe.

    Args:
        smiles: SMILES da molecula semente.
        property_name: propriedade a otimizar. "QED" (drug-likeness) ou "plogP".
        num_molecules: quantas moleculas devolver.
        minimize: False maximiza a propriedade, True minimiza.
        min_similarity: piso de similaridade de Tanimoto com a semente (0-1).
                        Alto = analogos conservadores; baixo = saltos maiores.
        particles: tamanho da populacao do CMA-ES; precisa ser >= num_molecules.
        iterations: iteracoes de otimizacao.
        name: rotulo curto usado no nome do arquivo.
    """
    smi = str(smiles).strip()
    if not smi:
        raise NimError("smiles vazio.")
    if not 0.0 <= min_similarity <= 1.0:
        raise NimError(f"min_similarity={min_similarity} fora de [0, 1].")
    if particles < num_molecules:
        raise NimError(
            f"particles ({particles}) precisa ser >= num_molecules ({num_molecules}): "
            "o CMA-ES nao devolve mais moleculas do que tem particulas."
        )

    payload: dict[str, Any] = {
        "smi": smi,
        "algorithm": "CMA-ES",
        "property_name": property_name,
        "num_molecules": int(num_molecules),
        "minimize": bool(minimize),
        "min_similarity": float(min_similarity),
        "particles": int(particles),
        "iterations": int(iterations),
    }
    data = await _post(REGISTRY["molmim"].path, payload)
    return _handle_molecule_response(
        data,
        stem=f"{name}_{_slug(smi)}",
        tool="generate_molecules",
        params={
            "seed_smiles": smi,
            "property_name": property_name,
            "minimize": minimize,
            "min_similarity": min_similarity,
            "iterations": iterations,
        },
    )


@mcp.tool()
async def generate_from_fragment(
    fragment_smiles: str,
    num_molecules: int = 10,
    temperature: float = 1.0,
    noise: float = 0.0,
    step_size: int = 1,
    scoring: str = "QED",
    unique: bool = True,
    name: str = "genmol",
) -> dict[str, Any]:
    """
    GenMol: completa um fragmento molecular por difusao discreta mascarada.

    A entrada usa a sintaxe de mascara do GenMol: `[*{N-N}]` marca onde crescer e
    quantos atomos adicionar. Ex.: "c1ccccc1[*{25-25}]" pede um benzeno com 25
    atomos novos pendurados. Sem mascara o modelo so devolve variacoes do
    proprio fragmento.

    Uso academico tipico: crescimento de fragmento (fragment-based design)
    partindo de um core conhecido.

    Args:
        fragment_smiles: fragmento em SMILES, com a mascara [*{N-N}].
        num_molecules: quantas moleculas gerar.
        temperature: aleatoriedade da amostragem; mais alto diversifica.
        noise: ruido de difusao.
        step_size: passos por iteracao de decodificacao.
        scoring: propriedade usada para ranquear a saida ("QED").
        unique: descartar SMILES repetidos. Com True o retorno pode vir com menos
                moleculas do que num_molecules — e esperado, nao e erro.
        name: rotulo curto usado no nome do arquivo.
    """
    frag = str(fragment_smiles).strip()
    if not frag:
        raise NimError("fragment_smiles vazio.")
    if "[*" not in frag:
        raise NimError(
            "fragment_smiles sem mascara de crescimento. Adicione [*{N-N}] onde a "
            "molecula deve crescer, ex.: 'c1ccccc1[*{25-25}]'."
        )

    payload: dict[str, Any] = {
        "smiles": frag,
        "num_molecules": int(num_molecules),
        "temperature": float(temperature),
        "noise": float(noise),
        "step_size": int(step_size),
        "scoring": scoring,
        "unique": bool(unique),
    }
    data = await _post(REGISTRY["genmol"].path, payload)
    return _handle_molecule_response(
        data,
        stem=f"{name}_{_slug(frag)}",
        tool="generate_from_fragment",
        params={
            "fragment": frag,
            "num_molecules": num_molecules,
            "temperature": temperature,
            "scoring": scoring,
            "unique": unique,
        },
    )


# ---------------------------------------------------------------------------
# Escape hatch generico
# ---------------------------------------------------------------------------


@mcp.tool()
async def nim_invoke(
    path: str, payload: dict[str, Any], save_as: str | None = None
) -> dict[str, Any]:
    """
    Chama qualquer endpoint NIM de biologia pelo caminho cru. Use quando o
    modelo nao estiver no registro, ou quando um slug tiver mudado.

    Args:
        path: caminho depois da base, ex. "biology/ipd/rfdiffusion/generate".
        payload: corpo JSON da requisicao.
        save_as: extensao para gravar a resposta em arquivo (ex. "json", "pdb").
                 Omita para devolver a resposta inline (cuidado com o tamanho).
    """
    # Validar a string crua nao basta: o httpx normaliza '..' depois, entao
    # "biology/../v1/keys" escaparia do prefixo. Normalize primeiro.
    parts = [seg for seg in str(path).split("/") if seg not in ("", ".")]
    if ".." in parts or not parts or parts[0] != "biology":
        raise NimError("por seguranca, path deve comecar com 'biology/' e nao conter '..'.")
    path = "/".join(parts)

    data = await _post(path, payload)

    if save_as:
        body = data if isinstance(data, str) else json.dumps(data, indent=2)
        p = _save(body, f"invoke_{_slug(path + json.dumps(payload, sort_keys=True))}", save_as)
        result = {
            "files": [str(p)],
            "response_keys": sorted(data.keys()) if isinstance(data, dict) else None,
        }
    else:
        text = json.dumps(data)
        if len(text) > 20_000:
            p = _save(json.dumps(data, indent=2), f"invoke_{_slug(path)}", "json")
            result = {
                "files": [str(p)],
                "truncated": True,
                "note": f"resposta de {len(text)} chars gravada em disco em vez de inline.",
                "response_keys": sorted(data.keys()) if isinstance(data, dict) else None,
            }
        else:
            result = {"response": data, "files": []}

    _log_provenance("nim_invoke", {"path": path}, result)
    return result


# ---------------------------------------------------------------------------
# Helper de resposta de estrutura
# ---------------------------------------------------------------------------


def _handle_structure_response(
    data: dict[str, Any], stem: str, tool: str, params: dict[str, Any]
) -> dict[str, Any]:
    """
    Normaliza as varias formas que os NIM usam para devolver estrutura e grava
    em disco. Retorna so caminhos + metricas de confianca.
    """
    files: list[str] = []
    ext = "cif"

    for key in (
        "structures",
        "structures_in_ranked_order",
        "pdbs",
        "output_pdb",
        "structure",
        "output",
    ):
        val = data.get(key)
        if isinstance(val, list) and val:
            for i, item in enumerate(val):
                body = item.get("structure") if isinstance(item, dict) else item
                if isinstance(body, str) and body.strip():
                    ext = "pdb" if body.lstrip().startswith(("ATOM", "HEADER")) else "cif"
                    files.append(str(_save(body, f"{stem}_m{i}", ext)))
            break
        if isinstance(val, str) and val.strip():
            ext = "pdb" if val.lstrip().startswith(("ATOM", "HEADER")) else "cif"
            files.append(str(_save(val, stem, ext)))
            break

    if not files:
        files.append(str(_save(json.dumps(data, indent=2), f"{stem}_raw", "json")))

    confidence = {
        k: data[k]
        for k in (
            "confidence_scores",
            "plddt",
            "ptm",
            "iptm",
            "affinity",
            "affinity_probability",
            "affinities",
        )
        if k in data
    }

    result: dict[str, Any] = {
        "files": files,
        "format": ext,
        "num_models": len(files),
        "confidence": confidence or None,
    }
    if "plddt" in confidence or "confidence_scores" in confidence:
        result["interpretation"] = (
            "pLDDT: >90 muito confiavel, 70-90 confiavel, 50-70 baixo, <50 provavelmente desordenado."
        )

    _log_provenance(tool, params, result)
    return result


if __name__ == "__main__":
    mcp.run()
