# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=1.2.0"]
# ///
"""
vcc-mcp — servidor MCP para o Virtual Cell Challenge (Arc Institute).

Por que este servidor existe
----------------------------
O VCC distribui uma CLI (`vcc`, pacote `vcc-cli`) e um *skill* para agentes
(`vcc skill install`). O skill e um guia em prosa: explica cada passo e lista
guardrails. O mais grave deles: na fase final, submeter a saida de
`vcc sample` substitui a entrada real da equipe, gasta uma das duas submissoes
do dia, e como os placares ficam ocultos ate o anuncio dos vencedores, nada
mostra que aconteceu.

Prosa nao segura isso. Este servidor transforma os guardrails em recusa de
codigo:

  - `vcc_submit` recusa, na fase final, qualquer arquivo gerado por
    `vcc_sample`. O reconhecimento e pelo sha256, entao renomear ou mover o
    arquivo nao engana.
  - submissao na fase final exige `confirmar_final=True`: so a ultima entrada
    final conta para premio, e cada uma pode substituir uma melhor.
  - `--skip-limit-check`, `--force` e `--endpoint` nao sao expostos. A propria
    CLI checa o limite diario antes do upload, sem gastar submissao. O skill
    oficial trata redirecionar o host como phishing, e por isso VCC_ENDPOINT e
    removida do ambiente de todo subprocesso.
  - so `.vcc` e aceito para submissao. O `.h5ad` cru passa antes por
    `vcc_prep`, que valida genes, contextos, alvos e celulas por perturbacao.

E segue a regra dos outros servidores do projeto: o que e grande vai para o
disco (controles em data/raw/vcc/, predicoes em runs/vcc/) e a ferramenta
devolve caminho e resumo. Toda chamada e registrada em runs/provenance.jsonl.

A credencial nunca passa por aqui. A CLI le a chave do chaveiro do sistema,
gravada uma unica vez por `vcc login --token-stdin`.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

# O SDK Python do MCP renomeou FastMCP -> MCPServer na versao 2.0. Suportar as
# duas evita que a stack quebre dependendo de qual versao o uv resolver.
try:
    from mcp.server.mcpserver import MCPServer as _Server  # mcp >= 2.0
    from mcp.server.mcpserver.exceptions import ToolError as _ToolError
except ImportError:  # pragma: no cover - caminho legado
    from mcp.server.fastmcp import FastMCP as _Server  # mcp 1.x
    from mcp.server.fastmcp.exceptions import ToolError as _ToolError

# ---------------------------------------------------------------------------
# Configuracao
# ---------------------------------------------------------------------------

ROOT = Path(os.environ.get("BIO_STACK_ROOT", Path(__file__).resolve().parents[2]))
PROV_LOG = Path(os.environ.get("BIO_STACK_PROVENANCE", ROOT / "runs" / "provenance.jsonl"))

# Controles sao dado bruto: entram em data/raw e dali em diante o guard_data.py
# impede remocao. Predicoes e amostras sao derivados e ficam em runs/vcc/,
# ignorado pelo git porque um .vcc completo tem ~240 MB.
DATA_DIR = ROOT / "data" / "raw" / "vcc"
RUNS_DIR = ROOT / "runs" / "vcc"
SAMPLE_DIR = RUNS_DIR / "sample"
REGISTRO_AMOSTRAS = SAMPLE_DIR / "registro.jsonl"
LIVRO_SUBMISSOES = RUNS_DIR / "submissoes.jsonl"

TIMEOUT_CURTO = 120.0
# ~630 MiB de controles ou ~240 MB de upload. As duas transferencias sao
# retomaveis pela CLI, entao estourar o prazo nao perde o progresso.
TIMEOUT_TRANSFERENCIA = 3600.0
# prep e sample processam 360 000 celulas x 18 533 genes.
TIMEOUT_PROCESSAMENTO = 1800.0

_PROV_AVISADO = False

mcp = _Server("vcc")


def _achar_vcc() -> str:
    """Caminho do executavel. Cliente MCP nem sempre herda o PATH do shell."""
    explicito = os.environ.get("VCC_BIN")
    if explicito:
        return explicito
    achado = shutil.which("vcc")
    if achado:
        return achado
    nome = "vcc.exe" if os.name == "nt" else "vcc"
    return str(Path.home() / ".local" / "bin" / nome)


VCC_BIN = _achar_vcc()

# ---------------------------------------------------------------------------
# Execucao da CLI
# ---------------------------------------------------------------------------


def _ultimo_json(texto: str) -> Any:
    """JSON do stdout. Se houver linhas de progresso, vale a ultima parseavel."""
    texto = texto.strip()
    if not texto:
        return None
    try:
        return json.loads(texto)
    except json.JSONDecodeError:
        pass
    for linha in reversed(texto.splitlines()):
        linha = linha.strip()
        if linha.startswith(("{", "[")):
            try:
                return json.loads(linha)
            except json.JSONDecodeError:
                continue
    return None


def _vcc(args: list[str], timeout: float, tolerar_falha: bool = False) -> dict[str, Any]:
    """
    Roda `vcc --json <args>` e devolve o JSON.

    Erro da CLI vira ToolError com o motivo que ela deu. Com `tolerar_falha`,
    o codigo de saida volta no resultado em vez de virar erro — e o caso do
    `status --wait`, em que falha de pontuacao e informacao, nao defeito.
    """
    ambiente = {k: v for k, v in os.environ.items() if k != "VCC_ENDPOINT"}
    try:
        proc = subprocess.run(
            [VCC_BIN, "--json", *args],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            env=ambiente,
            check=False,
        )
    except FileNotFoundError as exc:
        raise _ToolError(
            f"vcc nao encontrado em {VCC_BIN}. Instale com `uv tool install vcc-cli` "
            "ou aponte VCC_BIN para o executavel."
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise _ToolError(
            f"`vcc {args[0]}` passou de {timeout:.0f} s. Download e upload sao "
            "retomaveis: chame a ferramenta de novo e a CLI continua de onde parou."
        ) from exc

    dados = _ultimo_json(proc.stdout)
    erro_cli = isinstance(dados, dict) and "error" in dados
    if (proc.returncode != 0 or erro_cli) and not tolerar_falha:
        motivo = dados.get("error") if isinstance(dados, dict) else None
        codigo = dados.get("code") if isinstance(dados, dict) else None
        detalhe = motivo or proc.stderr.strip()[-800:] or proc.stdout.strip()[-800:]
        raise _ToolError(f"`vcc {args[0]}` falhou (saida {proc.returncode}, {codigo}): {detalhe}")
    if dados is None:
        if tolerar_falha:
            return {"saida": proc.returncode, "stderr": proc.stderr.strip()[-800:]}
        raise _ToolError(f"`vcc {args[0]}` nao devolveu JSON: {proc.stdout.strip()[-400:]}")
    if not isinstance(dados, dict):
        dados = {"resultado": dados}
    if tolerar_falha:
        dados["saida"] = proc.returncode
    return dados


def _identidade() -> dict[str, Any]:
    """Estado da conta no servidor do VCC. A fase vem dali, nunca de suposicao."""
    return _vcc(["whoami"], timeout=TIMEOUT_CURTO).get("identity", {})


# ---------------------------------------------------------------------------
# Pacote de controles
# ---------------------------------------------------------------------------


def _nome_pacote(final: bool) -> str:
    return "controls-final" if final else "controls"


def _arquivos_do_pacote(final: bool) -> tuple[Path, Path]:
    """gene_names.csv e pert_counts.csv do pacote da fase atual."""
    nome = _nome_pacote(final)
    base = DATA_DIR / nome
    perts = next(base.rglob("pert_counts.csv"), None) if base.exists() else None
    genes = next(base.rglob("gene_names.csv"), None) if base.exists() else None
    if perts is None or genes is None:
        raise _ToolError(
            f"Pacote `{nome}` ausente ou nao extraido em {base}. "
            f"Rode vcc_datasets_download(dataset='{nome}') primeiro."
        )
    # O prep le os contextos (D/E/F na final) do manifest.json AO LADO do
    # pert_counts.csv. Sem ele, valida contra os contextos errados.
    if not (perts.parent / "manifest.json").exists():
        raise _ToolError(
            f"manifest.json nao esta ao lado de {perts}. O prep precisa dele para "
            "saber os contextos da fase. Reextraia o pacote."
        )
    return genes, perts


def _extrair(zip_path: Path, destino: Path) -> dict[str, list[str]]:
    """
    Extrai sem sobrescrever. Dado bruto, uma vez em data/raw, e somente-leitura:
    arquivo existente com o mesmo tamanho e pulado; com tamanho diferente vira
    conflito reportado, nunca substituicao silenciosa.
    """
    extraidos: list[str] = []
    pulados: list[str] = []
    conflitos: list[str] = []
    raiz = destino.resolve()
    with zipfile.ZipFile(zip_path) as zf:
        for membro in zf.infolist():
            if membro.is_dir():
                continue
            alvo = (destino / membro.filename).resolve()
            # Zip slip: um membro com "../" escaparia de data/raw/vcc.
            if raiz not in alvo.parents:
                conflitos.append(f"{membro.filename} (caminho fora do destino, recusado)")
                continue
            rel = str(alvo.relative_to(raiz))
            if alvo.exists():
                if alvo.stat().st_size == membro.file_size:
                    pulados.append(rel)
                else:
                    conflitos.append(f"{rel} (ja existe com outro tamanho, mantido)")
                continue
            alvo.parent.mkdir(parents=True, exist_ok=True)
            with zf.open(membro) as origem, alvo.open("wb") as saida:
                shutil.copyfileobj(origem, saida, length=1 << 20)
            extraidos.append(rel)
    return {"extraidos": extraidos, "pulados": pulados, "conflitos": conflitos}


# ---------------------------------------------------------------------------
# Registros locais: amostras e submissoes
# ---------------------------------------------------------------------------


def _sha256(caminho: Path) -> str:
    h = hashlib.sha256()
    with caminho.open("rb") as fh:
        for bloco in iter(lambda: fh.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def _ler_jsonl(caminho: Path) -> list[dict[str, Any]]:
    if not caminho.exists():
        return []
    linhas: list[dict[str, Any]] = []
    for linha in caminho.read_text(encoding="utf-8").splitlines():
        try:
            linhas.append(json.loads(linha))
        except json.JSONDecodeError:
            continue
    return linhas


def _anexar_jsonl(caminho: Path, entrada: dict[str, Any]) -> None:
    caminho.parent.mkdir(parents=True, exist_ok=True)
    with caminho.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entrada, ensure_ascii=False) + "\n")


def _hashes_de_amostra() -> set[str]:
    return {e["sha256"] for e in _ler_jsonl(REGISTRO_AMOSTRAS) if "sha256" in e}


def _e_amostra(caminho: Path, sha: str) -> bool:
    try:
        caminho.resolve().relative_to(SAMPLE_DIR.resolve())
        return True
    except ValueError:
        return sha in _hashes_de_amostra()


def _submissoes_de_hoje() -> list[dict[str, Any]]:
    hoje = datetime.now(UTC).date().isoformat()
    return [e for e in _ler_jsonl(LIVRO_SUBMISSOES) if e.get("dia_utc") == hoje]


def _carimbo() -> str:
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def _log_provenance(tool: str, params: dict[str, Any], result: dict[str, Any]) -> None:
    """Registro append-only para reprodutibilidade academica."""
    try:
        _anexar_jsonl(
            PROV_LOG,
            {
                "ts": datetime.now(UTC).isoformat(),
                "server": "vcc-mcp",
                "tool": tool,
                "params": params,
                "outputs": result.get("arquivos", []),
            },
        )
    except OSError as exc:
        # O log nunca derruba a ferramenta, mas tambem nao some calado. stderr e
        # seguro: o protocolo MCP stdio usa so o stdout.
        global _PROV_AVISADO
        if not _PROV_AVISADO:
            print(
                f"[vcc] AVISO: procedencia NAO esta sendo gravada em {PROV_LOG}: {exc}",
                file=sys.stderr,
                flush=True,
            )
            _PROV_AVISADO = True


# ---------------------------------------------------------------------------
# Ferramentas
# ---------------------------------------------------------------------------


@mcp.tool()
def vcc_whoami() -> dict[str, Any]:
    """
    Conta, equipe, fase, contextos, prazo e se a equipe pode submeter agora.

    Chame antes de qualquer trabalho de dados ou submissao: a fase decide o
    pacote de controles (`controls` na validacao, `controls-final` na final) e
    as regras de submissao.
    """
    dados = _vcc(["whoami"], timeout=TIMEOUT_CURTO)
    ident = dados.get("identity", {})
    final = bool(ident.get("final_phase"))
    resultado: dict[str, Any] = {
        "conta": ident.get("display_name"),
        "organizacao": ident.get("organization"),
        "equipe": ident.get("team_name"),
        "fase": "final" if final else "validacao",
        "contextos": ident.get("contexts"),
        "painel": ident.get("panel_id"),
        "aceitando_submissoes": ident.get("accepting_submissions"),
        "pode_submeter": ident.get("can_submit"),
        "bloqueios": ident.get("blockers"),
        "credencial": dados.get("token_source"),
        "pacote_de_controles": _nome_pacote(final),
    }
    prazo = ident.get("deadline")
    if prazo:
        fim = datetime.fromisoformat(prazo.replace("Z", "+00:00"))
        resultado["prazo_utc"] = prazo
        resultado["horas_restantes"] = round((fim - datetime.now(UTC)).total_seconds() / 3600, 1)
    finais_hoje = sum(1 for s in _submissoes_de_hoje() if s.get("fase") == "final")
    resultado["finais_hoje_por_este_servidor"] = finais_hoje
    resultado["regra_de_cota"] = (
        "Fase final: 2 submissoes por equipe por dia UTC, e so a ULTIMA conta para "
        "premio. A contagem acima so ve o que passou por este servidor; a CLI checa "
        "o limite real antes do upload, sem gastar submissao."
    )
    _log_provenance("vcc_whoami", {}, resultado)
    return resultado


@mcp.tool()
def vcc_datasets_list() -> dict[str, Any]:
    """Pacotes de referencia que a conta pode baixar nesta fase."""
    resultado = _vcc(["datasets", "list"], timeout=TIMEOUT_CURTO)
    _log_provenance("vcc_datasets_list", {}, resultado)
    return resultado


@mcp.tool()
def vcc_datasets_download(dataset: str = "controls-final", extrair: bool = True) -> dict[str, Any]:
    """
    Baixa um pacote para data/raw/vcc/<dataset>/ e extrai.

    `controls-final` e o pacote da fase final (contextos D/E/F); `controls` e o
    da validacao (A/B/C). Cada um tem ~630 MiB: celulas controle, gene_names.csv,
    pert_counts.csv e manifest.json. Transferencia retomavel e com checksum
    verificado — se estourar o tempo, chame de novo.

    Nunca sobrescreve arquivo ja extraido.

    Args:
        dataset: identificador, como em vcc_datasets_list.
        extrair: extrair o zip depois do download.
    """
    ident = _identidade()
    final = bool(ident.get("final_phase"))
    destino = DATA_DIR / dataset
    destino.mkdir(parents=True, exist_ok=True)

    download = _vcc(["datasets", "download", dataset, "-d", str(destino)], timeout=TIMEOUT_TRANSFERENCIA)
    zips = sorted(destino.glob("*.zip"))
    resultado: dict[str, Any] = {
        "dataset": dataset,
        "diretorio": str(destino),
        "zips": [{"arquivo": z.name, "mb": round(z.stat().st_size / 2**20, 1)} for z in zips],
        "cli": download,
        "arquivos": [str(z) for z in zips],
    }
    if final and dataset == "controls":
        resultado["aviso"] = (
            "A fase atual e a FINAL. `controls` e o pacote da validacao (A/B/C). Um .vcc "
            "preparado com ele e recusado pelo `vcc submit` antes do upload. Use "
            "`controls-final`."
        )
    if extrair:
        resumo = {"extraidos": [], "pulados": [], "conflitos": []}
        for z in zips:
            parcial = _extrair(z, destino)
            for chave in resumo:
                resumo[chave].extend(parcial[chave])
        resultado["extracao"] = {
            "extraidos": len(resumo["extraidos"]),
            "pulados": len(resumo["pulados"]),
            "conflitos": resumo["conflitos"],
        }
        for nome in ("gene_names.csv", "pert_counts.csv", "manifest.json"):
            achado = next(destino.rglob(nome), None)
            resultado[nome] = str(achado) if achado else None
    _log_provenance("vcc_datasets_download", {"dataset": dataset, "extrair": extrair}, resultado)
    return resultado


@mcp.tool()
def vcc_prep(
    arquivo_h5ad: str,
    nome_saida: str | None = None,
    dry_run: bool = True,
    pert_col: str | None = None,
    context_col: str | None = None,
    celltype_col: str | None = None,
) -> dict[str, Any]:
    """
    Valida uma predicao .h5ad e, sem dry_run, empacota num .vcc submetivel.

    Valida contra o pacote da fase atual: 18 533 genes, contextos (D/E/F na
    final), alvos da lista oficial, exatamente 400 celulas por perturbacao em
    cada contexto, contagens cruas inteiras nao negativas, nenhuma celula
    controle. Rode primeiro com dry_run=True (o padrao): valida sem gravar.

    Args:
        arquivo_h5ad: caminho da predicao.
        nome_saida: nome do .vcc em runs/vcc/. Padrao: nome do h5ad + carimbo.
        dry_run: so validar, sem empacotar.
        pert_col: coluna de perturbacao em obs, se nao for a padrao.
        context_col: coluna de contexto em obs, se nao for a padrao.
        celltype_col: coluna de tipo celular em obs, se nao for a padrao.
    """
    entrada = Path(arquivo_h5ad).expanduser().resolve()
    if not entrada.exists():
        raise _ToolError(f"{entrada} nao existe.")
    if entrada.suffix.lower() != ".h5ad":
        raise _ToolError(f"vcc_prep espera um .h5ad, recebeu {entrada.name}.")

    final = bool(_identidade().get("final_phase"))
    genes, perts = _arquivos_do_pacote(final)
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    base = (nome_saida or f"{entrada.stem}_{_carimbo()}").removesuffix(".vcc")
    saida = RUNS_DIR / f"{base}.vcc"
    # Sem expor --force: nome repetido e recusado, nao sobrescrito.
    if saida.exists() and not dry_run:
        raise _ToolError(f"{saida} ja existe. Escolha outro nome_saida.")

    args = ["prep", str(entrada), "-g", str(genes), "--perts", str(perts), "-o", str(saida)]
    for flag, valor in (("--pert-col", pert_col), ("--context-col", context_col),
                        ("--celltype-col", celltype_col)):
        if valor:
            args += [flag, valor]
    if dry_run:
        args.append("--dry-run")

    cli = _vcc(args, timeout=TIMEOUT_PROCESSAMENTO)
    resultado: dict[str, Any] = {
        "fase": "final" if final else "validacao",
        "pacote": _nome_pacote(final),
        "dry_run": dry_run,
        "cli": cli,
        "arquivos": [],
    }
    if not dry_run and saida.exists():
        resultado["vcc"] = str(saida)
        resultado["mb"] = round(saida.stat().st_size / 2**20, 1)
        resultado["sha256"] = _sha256(saida)
        resultado["arquivos"] = [str(saida)]
        resultado["proximo_passo"] = "vcc_submit(arquivo_vcc=..., nome_modelo=...)"
    params = {"arquivo_h5ad": str(entrada), "nome_saida": nome_saida, "dry_run": dry_run,
              "pert_col": pert_col, "context_col": context_col, "celltype_col": celltype_col}
    _log_provenance("vcc_prep", params, resultado)
    return resultado


@mcp.tool()
def vcc_sample(nome: str = "teste") -> dict[str, Any]:
    """
    Gera um .vcc aleatorio e valido, so para exercitar o pipeline local.

    Ruido por construcao: pontua mal e nao diz nada sobre modelo. Fica em
    runs/vcc/sample/ e o sha256 entra num registro — vcc_submit recusa esse
    arquivo na fase final, mesmo renomeado ou movido, porque la ele
    substituiria a entrada real da equipe sem deixar rastro visivel.

    Args:
        nome: prefixo do arquivo.
    """
    final = bool(_identidade().get("final_phase"))
    genes, perts = _arquivos_do_pacote(final)
    SAMPLE_DIR.mkdir(parents=True, exist_ok=True)
    saida = SAMPLE_DIR / f"{nome}_{_carimbo()}.vcc"

    cli = _vcc(["sample", "-g", str(genes), "-p", str(perts), "-o", str(saida)],
               timeout=TIMEOUT_PROCESSAMENTO)
    sha = _sha256(saida)
    _anexar_jsonl(REGISTRO_AMOSTRAS, {"ts": datetime.now(UTC).isoformat(),
                                       "arquivo": str(saida), "sha256": sha})
    resultado: dict[str, Any] = {
        "vcc": str(saida),
        "mb": round(saida.stat().st_size / 2**20, 1),
        "sha256": sha,
        "cli": cli,
        "arquivos": [str(saida)],
        "aviso": (
            "Arquivo de TESTE, registrado como amostra. Na fase final vcc_submit o "
            "recusa. Serve para testar prep e leitura local, nao para submeter."
        ),
    }
    _log_provenance("vcc_sample", {"nome": nome}, resultado)
    return resultado


@mcp.tool()
def vcc_submit(
    arquivo_vcc: str,
    nome_modelo: str,
    descricao: str | None = None,
    confirmar_final: bool = False,
) -> dict[str, Any]:
    """
    Envia um .vcc e inicia a pontuacao.

    Recusas, nesta ordem, antes de qualquer upload:
      1. arquivo que nao e .vcc — o .h5ad passa antes por vcc_prep;
      2. na fase final, arquivo gerado por vcc_sample (reconhecido pelo sha256);
      3. conta que o servidor diz nao poder submeter agora;
      4. na fase final, sem confirmar_final=True.

    Fase final: 2 submissoes por equipe por dia UTC, uma em andamento por vez
    (a segunda da HTTP 409), e so a ULTIMA entrada final conta para premio —
    enviar substitui a anterior, mesmo que a anterior fosse melhor. Os placares
    ficam ocultos ate o anuncio: `published` com todos os scores nulos e o
    resultado correto.

    Args:
        arquivo_vcc: caminho do .vcc (saida de vcc_prep).
        nome_modelo: rotulo da submissao (obrigatorio na CLI).
        descricao: ate 2000 caracteres.
        confirmar_final: obrigatorio na fase final.
    """
    caminho = Path(arquivo_vcc).expanduser().resolve()
    if not caminho.exists():
        raise _ToolError(f"{caminho} nao existe.")
    if caminho.suffix.lower() != ".vcc":
        raise _ToolError(
            f"vcc_submit so aceita .vcc, recebeu {caminho.name}. Rode vcc_prep no .h5ad "
            "primeiro: ele valida contextos, alvos e celulas por perturbacao."
        )
    if not nome_modelo.strip():
        raise _ToolError("nome_modelo e obrigatorio.")
    if descricao and len(descricao) > 2000:
        raise _ToolError(f"descricao tem {len(descricao)} caracteres; o limite e 2000.")

    sha = _sha256(caminho)
    ident = _identidade()
    final = bool(ident.get("final_phase"))
    amostra = _e_amostra(caminho, sha)

    if final and amostra:
        raise _ToolError(
            "RECUSADO: este arquivo e uma amostra de vcc_sample (esta em runs/vcc/sample/ "
            "ou tem o sha256 no registro de amostras). Na "
            "fase final ele viraria a entrada final da equipe, SUBSTITUINDO a real, "
            "gastaria uma das 2 submissoes do dia, e com os placares ocultos nada "
            "mostraria o estrago. Submeta uma predicao de verdade preparada por vcc_prep."
        )
    if not ident.get("accepting_submissions") or not ident.get("can_submit"):
        raise _ToolError(
            f"O servidor do VCC diz que a equipe nao pode submeter agora. "
            f"aceitando={ident.get('accepting_submissions')}, "
            f"pode_submeter={ident.get('can_submit')}, bloqueios={ident.get('blockers')}."
        )
    if final and not confirmar_final:
        anteriores = [s for s in _submissoes_de_hoje() if s.get("fase") == "final"]
        raise _ToolError(
            "Fase final: confirme com confirmar_final=True. Esta submissao passa a ser "
            "a entrada final da equipe e SUBSTITUI a anterior, mesmo que ela fosse "
            f"melhor — e os placares estao ocultos. Este servidor ja registrou "
            f"{len(anteriores)} submissao(oes) final(is) hoje (UTC); o limite e 2."
        )

    args = ["submit", str(caminho), "-m", nome_modelo]
    if descricao:
        args += ["-d", descricao]
    cli = _vcc(args, timeout=TIMEOUT_TRANSFERENCIA)

    entrada_id = next(
        (cli.get(k) for k in ("entry_id", "id", "entry", "submission_id") if cli.get(k)), None
    )
    registro = {
        "ts": datetime.now(UTC).isoformat(),
        "dia_utc": datetime.now(UTC).date().isoformat(),
        "fase": "final" if final else "validacao",
        "entry_id": entrada_id,
        "nome_modelo": nome_modelo,
        "arquivo": str(caminho),
        "sha256": sha,
        "amostra": amostra,
    }
    _anexar_jsonl(LIVRO_SUBMISSOES, registro)
    resultado: dict[str, Any] = {
        "entry_id": entrada_id,
        "fase": registro["fase"],
        "sha256": sha,
        "cli": cli,
        "arquivos": [str(caminho)],
        "proximo_passo": "vcc_status(entry_id=..., esperar=True)",
    }
    if amostra:
        resultado["aviso"] = "Arquivo de amostra (fase de validacao): pontua mal por construcao."
    _log_provenance(
        "vcc_submit",
        {"arquivo_vcc": str(caminho), "nome_modelo": nome_modelo, "confirmar_final": confirmar_final},
        resultado,
    )
    return resultado


@mcp.tool()
def vcc_status(entry_id: str, esperar: bool = False, timeout_espera_s: int = 900) -> dict[str, Any]:
    """
    Estado (e placares, quando publicados) de uma submissao.

    Tres resultados que parecem falha e nao sao, na fase final: `published` com
    todos os scores nulos (fase cega); `superseded` (outra entrada terminou de
    pontuar depois); e "nao encontrado" para uma entrada antiga, que e removida
    quando uma mais nova e pontuada.

    Args:
        entry_id: identificador devolvido por vcc_submit.
        esperar: bloquear ate estado terminal.
        timeout_espera_s: teto da espera.
    """
    args = ["status", entry_id]
    timeout = TIMEOUT_CURTO
    if esperar:
        args += ["--wait", "--wait-timeout", str(timeout_espera_s)]
        timeout = timeout_espera_s + TIMEOUT_CURTO
    resultado = _vcc(args, timeout=timeout, tolerar_falha=True)
    resultado["leitura"] = (
        "Fase final: scores nulos, `superseded` e entrada antiga nao encontrada sao "
        "esperados. Nao reenvie para ver placar — gastaria submissao e substituiria a entrada."
    )
    _log_provenance("vcc_status", {"entry_id": entry_id, "esperar": esperar}, resultado)
    return resultado


@mcp.tool()
def vcc_cancel(entry_id: str) -> dict[str, Any]:
    """
    Abandona uma submissao em andamento e libera a vaga da equipe.

    So uma submissao por equipe roda por vez; cancelar e o caminho para uma que
    travou no upload. Para continuar um upload interrompido em vez de abandonar,
    use a CLI direto: `vcc submit --resume <entry_id>`.

    Args:
        entry_id: identificador da submissao.
    """
    resultado = _vcc(["cancel", entry_id], timeout=TIMEOUT_CURTO)
    _log_provenance("vcc_cancel", {"entry_id": entry_id}, resultado)
    return resultado


if __name__ == "__main__":
    mcp.run()
