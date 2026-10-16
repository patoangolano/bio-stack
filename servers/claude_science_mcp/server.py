# /// script
# requires-python = ">=3.11"
# dependencies = ["mcp>=1.2.0", "httpx>=0.27.0"]
# ///
"""
claude-science-mcp — opera o daemon local do Claude Science pela API dele.

Por que este servidor existe
----------------------------
O Claude Science roda como daemon na maquina, com UI web em 127.0.0.1. O
`claude-science --help` nao tem subcomando para rodar uma conversa: serve,
open, url, status, logs, stop, update, import, sandbox, install. Operar o app,
entao, exigiria um humano digitando no navegador.

Mas o daemon expoe uma API REST completa — `/api/frames`, `/api/frames/:id/
message`, `/api/projects`, `/api/agents` — e essa API aceita a mesma sessao que
a UI usa. Este servidor fecha a ponte: o Claude Code manda trabalho para o
Claude Science e le a resposta, sem navegador no meio.

A divisao de trabalho que isso habilita
---------------------------------------
O Claude Science tem o que esta sessao nao tem: sandbox de rede com lista de
dominios aprovada, 24 conectores MCP nativos de biologia (biomart, variants,
clinical-genomics, structures-interactions, literature, genomes, chemistry...)
e um orquestrador com agente revisor. Serve para varredura longa e para
trabalho que nao deve consumir o contexto desta sessao.

O caminho inverso tambem existe e ja esta montado: o conector `nvidia-bio`
registrado em ~/.claude-science/mcp/local-mcp.json da ao Claude Science as 14
ferramentas dos NIM hospedados da NVIDIA.

Autenticacao, e por que ela e assim
-----------------------------------
A API nao aceita bearer arbitrario: a tentativa devolve
`{"detail":"invalid bearer token"}`. O fluxo real, medido em 2026-10-11 contra
a versao 0.1.60, tem quatro passos e cada um foi necessario:

  1. `claude-science url` emite um nonce de uso unico, valido ~3 min.
  2. `GET /?nonce=<n>` devolve a pagina de login e semeia o cookie
     `operon_csrf`.
  3. `POST /api/auth/nonce` com o corpo `{"nonce": <n>}` conclui o sign-in e
     semeia `operon_auth`. Nao e `/api/auth/exchange` — essa rota exige um
     bearer que o fluxo de navegador nunca produz.
  4. Dai em diante, toda ESCRITA precisa de dois cabecalhos:
       Origin: http://127.0.0.1:<porta>   sem ele: 403 origin_required
       X-Operon-CSRF: <cookie operon_csrf>  com X-CSRF-Token: 403 csrf_stale

O nome do cabecalho de CSRF e a armadilha central. `X-CSRF-Token`,
`x-csrf-token`, `X-XSRF-TOKEN` e `csrf-token` todos devolvem
`csrf_stale` — a mensagem sugere token velho, mas o token estava correto e o
*nome do cabecalho* estava errado. Com `X-Operon-CSRF` a requisicao passa
direto para a validacao de corpo.

A segunda armadilha e o caso das chaves. O corpo e snake_case
(`project_id`, `input_data`), mas o validador relata o erro em camelCase:
pedir `project_id` enquanto diz `projectId: project_id is required`. E
`input_data` e `record(unknown())` — aceita qualquer objeto sem validar o
interior, entao `{"prompt": "..."}` retorna 200 e **falha depois**, no modelo,
com `messages: at least one message is required`. A chave interna correta e
`request`.

Cobranca, e a regra que este servidor impoe
-------------------------------------------
O daemon autentica com `provider: claude_ai` e reporta
`subscription_type: max`, `rate_limit_tier: default_claude_max_20x`,
`billing_type: stripe_subscription`. Ou seja: o uso sai da assinatura Max, nao
de credito por token. O saldo de credito marca 0 BRL e as conversas rodam
normalmente — se fosse cobranca por token, falhariam.

O campo `total_cost` que as ferramentas devolvem e **contabilidade de tokens**,
nao fatura: mede quanto aquela conversa valeria em API, util para saber o peso
de um pedido. A conversa de verificacao marcou US$ 0,70 com 548 mil tokens de
entrada, e nada disso foi cobrado separado.

Mukhamed pediu explicitamente para nunca gastar por token. Isso esta imposto em
codigo, nao em disciplina: `_exige_assinatura` le /api/auth/status antes de cada
escrita e **recusa** se a conta nao estiver em modo assinatura. Se um dia o
daemon cair para credito pre-pago, as ferramentas param em vez de gastar.
Burlar exige `CLAUDE_SCIENCE_PERMITIR_PAGO=1`, deliberado e explicito.

O recurso escasso de verdade e a janela de rate limit do Max, compartilhada com
esta sessao do Claude Code. Um pedido grande no Claude Science consome parte
dela — por isso as ferramentas reportam tokens, e vale mandar para la o que e
longo, nao o que e trivial.
"""

from __future__ import annotations

import json
import logging
import os
import re
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx

# O httpx loga cada requisicao em INFO. Num servidor MCP stdio isso vai para o
# stderr do cliente e enche o log com a URL de login — que carrega o nonce.
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("httpcore").setLevel(logging.WARNING)

# O SDK Python do MCP renomeou FastMCP -> MCPServer na versao 2.0. Suportar as
# duas evita que a stack quebre dependendo de qual versao o uv resolver.
try:
    from mcp.server.mcpserver import MCPServer as _Server  # mcp >= 2.0
    from mcp.server.mcpserver.exceptions import ToolError as _ToolError
except ImportError:  # pragma: no cover - caminho legado
    from mcp.server.fastmcp import FastMCP as _Server  # mcp 1.x
    from mcp.server.fastmcp.exceptions import ToolError as _ToolError


class ScienceError(_ToolError):
    """Erro de operacao do Claude Science.

    Herda de ToolError porque, a partir do MCP 2.3.0, so a mensagem de um
    ToolError sobrevive ate o cliente. Qualquer outra excecao chega ao modelo
    como "erro interno", sem o texto que diria o que consertar.
    """


# ---------------------------------------------------------------------------
# Configuracao
# ---------------------------------------------------------------------------

RAIZ = Path(__file__).resolve().parents[2]
PROV_LOG = Path(os.environ.get("BIO_STACK_PROVENANCE", RAIZ / "runs" / "provenance.jsonl"))

# Onde a sessao autenticada e guardada entre chamadas. Sem cache, cada
# ferramenta queimaria um nonce de uso unico e um `wsl.exe` de ~2 s.
SESSAO_CACHE = Path(
    os.environ.get(
        "CLAUDE_SCIENCE_SESSION_CACHE", Path.home() / ".claude-science" / ".mcp-sessao.json"
    )
)

DISTRO = os.environ.get("CLAUDE_SCIENCE_WSL_DISTRO", "Ubuntu-24.04")
EXE_WINDOWS = Path(
    os.environ.get(
        "CLAUDE_SCIENCE_EXE",
        Path.home() / "AppData" / "Local" / "Programs" / "ClaudeScience" / "claude-science.exe",
    )
)
# Qual daemon operar quando os dois estao de pe: "wsl", "windows" ou "auto".
# O padrao e wsl porque e la que o conector nvidia-bio esta registrado.
ALVO = os.environ.get("CLAUDE_SCIENCE_TARGET", "wsl").lower()

# Escape hatch da regra de assinatura. So ligue sabendo que vai gastar por token.
PERMITIR_PAGO = os.environ.get("CLAUDE_SCIENCE_PERMITIR_PAGO", "") == "1"

AGENTE_PADRAO = os.environ.get("CLAUDE_SCIENCE_AGENT", "OPERON")
ESPERA_PADRAO_S = int(os.environ.get("CLAUDE_SCIENCE_TIMEOUT", "900"))
INTERVALO_POLL_S = 3.0

# Estados TERMINAIS. Tudo que nao esta aqui conta como "ainda trabalhando".
#
# A primeira versao fazia o contrario — listava os estados em andamento — e
# custou uma conversa paga: o daemon devolveu "processing", que nao estava na
# lista, o poll tratou como terminal e a ferramenta voltou com resposta vazia
# enquanto o agente seguia rodando e cobrando. Enumerar o fim e mais seguro
# que enumerar o meio: um estado desconhecido no meio so faz esperar, e a
# espera tem teto em timeout_s; um estado desconhecido tratado como fim
# descarta trabalho ja pago.
TERMINAIS = frozenset({"completed", "failed", "error", "cancelled", "canceled"})

mcp = _Server("claude-science")


# ---------------------------------------------------------------------------
# Falar com o CLI do daemon
# ---------------------------------------------------------------------------


def _cli(*args: str, timeout: int = 90) -> subprocess.CompletedProcess[str]:
    """Roda o CLI do claude-science no alvo escolhido.

    No WSL o comando vai por `bash -lc`, nao `bash -c`: o PATH de
    ~/.local/bin esta em ~/.profile, lido so em shell de login. Sem o -l, o
    binario simplesmente nao existe.
    """
    if _alvo_resolvido() == "windows":
        if not EXE_WINDOWS.is_file():
            raise ScienceError(
                f"nao achei o executavel do Claude Science em {EXE_WINDOWS}. "
                "Aponte CLAUDE_SCIENCE_EXE para ele ou use CLAUDE_SCIENCE_TARGET=wsl."
            )
        cmd = [str(EXE_WINDOWS), *args]
    else:
        # Sem cifrao no comando: chamado por subprocess, o wsl.exe perde as
        # sequencias $var na traducao da linha de comando, em silencio.
        cmd = ["wsl.exe", "-d", DISTRO, "--", "bash", "-lc", "claude-science " + " ".join(args)]
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except FileNotFoundError as exc:
        raise ScienceError(f"nao consegui executar o CLI do Claude Science: {exc}") from exc
    except subprocess.TimeoutExpired as exc:
        raise ScienceError(f"o CLI do Claude Science nao respondeu em {timeout}s") from exc


_ALVO_CACHE: str | None = None


def _alvo_resolvido() -> str:
    """ "wsl" ou "windows". Com ALVO=auto, escolhe o que estiver rodando."""
    global _ALVO_CACHE
    if ALVO in ("wsl", "windows"):
        return ALVO
    if _ALVO_CACHE:
        return _ALVO_CACHE
    for candidato in ("wsl", "windows"):
        try:
            if candidato == "windows" and not EXE_WINDOWS.is_file():
                continue
            cmd = (
                [str(EXE_WINDOWS), "status"]
                if candidato == "windows"
                else ["wsl.exe", "-d", DISTRO, "--", "bash", "-lc", "claude-science status"]
            )
            saida = subprocess.run(cmd, capture_output=True, text=True, timeout=60, check=False)
            if json.loads(saida.stdout or "{}").get("running"):
                _ALVO_CACHE = candidato
                return candidato
        except (OSError, ValueError, subprocess.SubprocessError):
            continue
    raise ScienceError(
        "nenhum daemon do Claude Science esta rodando. Suba um: "
        f'wsl -d {DISTRO} -- bash -lc "claude-science serve --port 8765 --detached --no-browser"'
    )


def _status_daemon() -> dict[str, Any]:
    r = _cli("status")
    try:
        return json.loads(r.stdout or "{}")
    except json.JSONDecodeError as exc:
        raise ScienceError(
            f"`claude-science status` nao devolveu JSON: {r.stdout[:200]!r}"
        ) from exc


def _base_url() -> str:
    if forcado := os.environ.get("CLAUDE_SCIENCE_URL"):
        return forcado.rstrip("/")
    st = _status_daemon()
    if not st.get("running"):
        raise ScienceError(
            "o daemon do Claude Science nao esta rodando. Suba com "
            "`claude-science serve --detached --no-browser` no alvo escolhido."
        )
    return f"http://127.0.0.1:{st['port']}"


# ---------------------------------------------------------------------------
# Sessao
# ---------------------------------------------------------------------------


def _nonce() -> str:
    """Emite um nonce de login. Uso unico, ~3 min de validade."""
    r = _cli("url")
    m = re.search(r"nonce=([0-9a-f]{32,})", r.stdout + r.stderr)
    if not m:
        raise ScienceError(
            "`claude-science url` nao devolveu um nonce. "
            f"stdout={r.stdout[:160]!r} stderr={r.stderr[:160]!r}"
        )
    return m.group(1)


def _cookies_do_cache(base: str) -> dict[str, str] | None:
    try:
        d = json.loads(SESSAO_CACHE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    # A sessao e por daemon: reusar a do daemon do WSL contra o do Windows
    # autenticaria contra o banco errado.
    if d.get("base") != base:
        return None
    return d.get("cookies") or None


def _grava_cache(base: str, cookies: dict[str, str]) -> None:
    try:
        SESSAO_CACHE.parent.mkdir(parents=True, exist_ok=True)
        SESSAO_CACHE.write_text(
            json.dumps({"base": base, "cookies": cookies, "ts": datetime.now(UTC).isoformat()}),
            encoding="utf-8",
        )
        # O cookie operon_auth e credencial de sessao: vale a mesma coisa que
        # estar logado na UI. 0600 impede outro usuario da maquina de ler.
        SESSAO_CACHE.chmod(0o600)
    except OSError as exc:
        print(f"[claude-science] AVISO: sessao nao cacheada: {exc}", file=sys.stderr, flush=True)


def _autentica(c: httpx.Client, base: str) -> None:
    # Jar limpo antes do sign-in. Os cookies vindos do cache entram sem dominio;
    # os que o daemon semeia agora vem com dominio 127.0.0.1. Os dois convivem no
    # jar com o mesmo nome e cookies.get("operon_csrf") levanta CookieConflict —
    # medido em 2026-10-16, quando a distro foi refeita e o daemon novo devolveu
    # 401 para a sessao em cache. Quem chega aqui tem sessao morta: nada a perder.
    c.cookies.clear()
    n = _nonce()
    c.get("/", params={"nonce": n})
    csrf = c.cookies.get("operon_csrf")
    if not csrf:
        raise ScienceError("o daemon nao semeou o cookie operon_csrf na pagina de login")
    r = c.post(
        "/api/auth/nonce", json={"nonce": n}, headers={"X-Operon-CSRF": csrf, "Origin": base}
    )
    if r.status_code != 200 or not c.cookies.get("operon_auth"):
        raise ScienceError(
            f"sign-in falhou (HTTP {r.status_code}). O nonce vale ~3 min e e de uso unico; "
            "se outra aba consumiu, tente de novo."
        )
    _grava_cache(base, dict(c.cookies))


class Sessao:
    """Cliente HTTP autenticado. Reautentica uma vez em caso de 401."""

    def __init__(self, timeout: float = 120.0) -> None:
        self.base = _base_url()
        self.c = httpx.Client(base_url=self.base, timeout=timeout, follow_redirects=False)
        if cookies := _cookies_do_cache(self.base):
            self.c.cookies.update(cookies)
        else:
            _autentica(self.c, self.base)

    def __enter__(self) -> Sessao:
        return self

    def __exit__(self, *_: object) -> None:
        self.c.close()

    def _cab(self) -> dict[str, str]:
        # O CSRF e lido do cookie a cada chamada: ele pode girar, e um valor
        # guardado de antes do sign-in produz csrf_stale.
        return {"X-Operon-CSRF": self.c.cookies.get("operon_csrf") or "", "Origin": self.base}

    def pedir(self, metodo: str, rota: str, **kw: Any) -> httpx.Response:
        kw.setdefault("headers", {}).update(self._cab())
        try:
            r = self.c.request(metodo, rota, **kw)
        except httpx.HTTPError as exc:
            raise ScienceError(
                f"{metodo} {rota}: falha de rede: {type(exc).__name__}: {exc}"
            ) from exc
        if r.status_code == 401:
            # Cookie expirado: refaz o sign-in e repete uma vez.
            _autentica(self.c, self.base)
            kw["headers"].update(self._cab())
            r = self.c.request(metodo, rota, **kw)
        if r.status_code >= 400:
            raise ScienceError(f"{metodo} {rota} -> HTTP {r.status_code}: {r.text[:400]}")
        return r

    def json(self, metodo: str, rota: str, **kw: Any) -> Any:
        r = self.pedir(metodo, rota, **kw)
        if not r.content:
            return None
        try:
            return r.json()
        except json.JSONDecodeError as exc:
            raise ScienceError(f"{metodo} {rota} nao devolveu JSON: {r.text[:200]!r}") from exc


# ---------------------------------------------------------------------------
# Procedencia
# ---------------------------------------------------------------------------

_PROV_AVISADO = False


def _log_provenance(tool: str, params: dict[str, Any], result: dict[str, Any]) -> None:
    """Registro append-only, mesmo formato dos outros servidores da stack."""
    global _PROV_AVISADO
    try:
        PROV_LOG.parent.mkdir(parents=True, exist_ok=True)
        entry = {
            "ts": datetime.now(UTC).isoformat(),
            "server": "claude-science-mcp",
            "tool": tool,
            "mode": _alvo_resolvido(),
            "params": params,
            "resumo": {
                k: result.get(k)
                for k in ("frame_id", "status", "custo_usd", "tokens_entrada", "tokens_saida")
                if k in result
            },
        }
        with PROV_LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except OSError as exc:
        if not _PROV_AVISADO:
            print(
                f"[claude-science] AVISO: procedencia NAO gravada em {PROV_LOG}: {exc}",
                file=sys.stderr,
                flush=True,
            )
            _PROV_AVISADO = True


# ---------------------------------------------------------------------------
# Helpers de dominio
# ---------------------------------------------------------------------------


def _texto(conteudo: Any) -> str:
    """Achata o conteudo de uma mensagem em texto.

    As mensagens chegam como lista de blocos tipados ([{"type":"text",...}]),
    nao como string. Blocos sem texto (uso de ferramenta, imagem) viram um
    marcador curto em vez de JSON cru, que inundaria o contexto.
    """
    if isinstance(conteudo, str):
        return conteudo
    if isinstance(conteudo, list):
        partes = []
        for b in conteudo:
            if not isinstance(b, dict):
                partes.append(str(b))
            elif "text" in b:
                partes.append(str(b["text"]))
            else:
                partes.append(f"<{b.get('type', 'bloco')}>")
        return "\n".join(partes)
    return json.dumps(conteudo, ensure_ascii=False)


def _exige_assinatura(s: Sessao) -> dict[str, Any]:
    """Recusa qualquer escrita se a conta nao estiver em modo assinatura.

    Mukhamed usa o plano Max e pediu para nunca gastar por token. Checar isso a
    cada escrita custa uma requisicao local e evita a falha silenciosa que
    importa: o daemon mudar de modo (assinatura expirada, troca de org, conta
    deslogada e refeita com chave de API) e as conversas passarem a consumir
    credito sem nada no caminho avisando.
    """
    auth = s.json("GET", "/api/auth/status") or {}
    assinatura = (auth.get("subscription_type") or "").lower()
    cobranca = (auth.get("billing_type") or "").lower()
    ok = assinatura in ("max", "pro", "team", "enterprise") or "subscription" in cobranca
    if not ok and not PERMITIR_PAGO:
        raise ScienceError(
            "recusado: a conta do Claude Science nao esta em modo assinatura "
            f"(subscription_type={auth.get('subscription_type')!r}, "
            f"billing_type={auth.get('billing_type')!r}), entao esta conversa "
            "sairia como cobranca por token. Reentre com a conta do plano Max na UI "
            "do Claude Science. Para gastar por token de proposito, suba o servidor "
            "com CLAUDE_SCIENCE_PERMITIR_PAGO=1."
        )
    return auth


def _resolve_projeto(s: Sessao, projeto: str | None) -> str:
    projs = (s.json("GET", "/api/projects") or {}).get("projects", [])
    if not projs:
        raise ScienceError("o Claude Science nao tem nenhum projeto. Crie um na UI primeiro.")
    if not projeto:
        return projs[0]["project_id"]
    for p in projs:
        if projeto in (p["project_id"], p.get("name")):
            return p["project_id"]
    nomes = ", ".join(f"{p.get('name')} ({p['project_id']})" for p in projs)
    raise ScienceError(f"projeto {projeto!r} nao existe. Disponiveis: {nomes}")


def _frame(s: Sessao, frame_id: str) -> dict[str, Any]:
    return s.json("GET", f"/api/frames/{frame_id}") or {}


def _espera(s: Sessao, frame_id: str, timeout_s: int) -> dict[str, Any]:
    limite = time.monotonic() + timeout_s
    f: dict[str, Any] = {}
    while time.monotonic() < limite:
        f = _frame(s, frame_id)
        if f.get("status") in TERMINAIS:
            return f
        time.sleep(INTERVALO_POLL_S)
    raise ScienceError(
        f"o frame {frame_id} ainda estava em '{f.get('status')}' depois de {timeout_s}s. "
        f"Ele continua rodando no Claude Science: leia com science_read(frame_id='{frame_id}') "
        f"ou interrompa com science_cancel."
    )


def _mensagens(s: Sessao, frame_id: str) -> list[dict[str, Any]]:
    ms = s.json("GET", f"/api/frames/{frame_id}/messages") or []
    return ms if isinstance(ms, list) else ms.get("messages", [])


def _resposta_final(ms: list[dict[str, Any]]) -> str:
    for m in reversed(ms):
        if m.get("role") == "assistant":
            return _texto(m.get("content")).strip()
    return ""


def _custos(f: dict[str, Any]) -> dict[str, Any]:
    return {
        "custo_usd": f.get("total_cost"),
        "tokens_entrada": f.get("input_tokens"),
        "tokens_saida": f.get("output_tokens"),
        "contexto_pct": f.get("context_usage_percent"),
    }


# ---------------------------------------------------------------------------
# Ferramentas
# ---------------------------------------------------------------------------


@mcp.tool()
def science_status() -> dict[str, Any]:
    """Estado do daemon do Claude Science: porta, versao, sandbox, conta e saldo.

    Comece por aqui. Nenhuma outra ferramenta funciona sem o daemon de pe.
    """
    st = _status_daemon()
    if not st.get("running"):
        return {
            "rodando": False,
            "alvo": _alvo_resolvido(),
            "como_subir": (
                f'wsl -d {DISTRO} -- bash -lc "claude-science serve '
                '--port 8765 --detached --no-browser"'
            ),
        }
    d = st.get("daemon", {})
    out: dict[str, Any] = {
        "rodando": True,
        "alvo": _alvo_resolvido(),
        "porta": st.get("port"),
        "pid": st.get("pid"),
        "versao": st.get("version"),
        "sandbox_ativa": d.get("sandbox_active"),
        "exige_token": d.get("require_token"),
        "conversas_ativas": d.get("active_conversations"),
        "data_dir": d.get("data_dir"),
    }
    with Sessao(timeout=30) as s:
        auth = s.json("GET", "/api/auth/status") or {}
        out |= {"conta": auth.get("email"), "autenticado": auth.get("authenticated")}
        try:
            saldo = s.json("GET", "/api/usage/credits-balance") or {}
            dados = saldo.get("data") or {}
            out["saldo"] = f"{dados.get('amount')} {dados.get('currency')}"
        except ScienceError:
            out["saldo"] = "indisponivel"
    return out


@mcp.tool()
def science_projects() -> list[dict[str, Any]]:
    """Projetos do Claude Science, com id para passar em science_ask."""
    with Sessao(timeout=30) as s:
        projs = (s.json("GET", "/api/projects") or {}).get("projects", [])
    return [
        {
            "project_id": p["project_id"],
            "nome": p.get("name"),
            "conversas": p.get("conversation_count"),
            "artefatos": p.get("artifact_count"),
        }
        for p in projs
    ]


@mcp.tool()
def science_agents() -> list[dict[str, Any]]:
    """Agentes disponiveis. OPERON e o orquestrador e o padrao."""
    with Sessao(timeout=30) as s:
        ags = s.json("GET", "/api/agents") or []
    return [
        {
            "nome": a.get("name"),
            "descricao": (a.get("description") or "").split("\n")[0][:200],
            "habilitado": a.get("enabled"),
        }
        for a in ags
    ]


@mcp.tool()
def science_ask(
    pedido: str,
    projeto: str | None = None,
    agente: str = AGENTE_PADRAO,
    nome: str | None = None,
    timeout_s: int = ESPERA_PADRAO_S,
) -> dict[str, Any]:
    """Abre uma conversa no Claude Science, manda o pedido e espera a resposta.

    Consome a janela de rate limit do plano Max — a mesma que esta sessao usa —
    e nao cobranca por token: o servidor recusa se a conta nao estiver em modo
    assinatura. O `custo_usd` devolvido e contabilidade de tokens, para dimen-
    sionar o pedido, nao fatura.

    Vale a pena quando o trabalho (a) e longo e nao deve ocupar o contexto desta
    sessao, (b) precisa dos 24 conectores nativos de biologia do Claude Science,
    ou (c) precisa da sandbox de rede dele. Pergunta trivial se responde aqui.

    Peca a saida enxuta no proprio texto — "devolva so a tabela", "nao imprima
    sequencia". Sem isso o orquestrador despeja o material inteiro na resposta,
    e ela volta para o contexto desta sessao.
    """
    pedido = pedido.strip()
    if not pedido:
        raise ScienceError("pedido vazio")

    with Sessao(timeout=max(120, timeout_s)) as s:
        _exige_assinatura(s)
        pid = _resolve_projeto(s, projeto)
        criado = s.json(
            "POST",
            "/api/frames",
            json={
                # snake_case: o validador relata o erro em camelCase
                # ("projectId: project_id is required"), mas a chave aceita e
                # project_id. Nao reescreva para projectId.
                "project_id": pid,
                "agent_name": agente,
                "name": nome or f"[bio-stack] {pedido[:60]}",
            },
        )
        fid = criado["root_frame_id"]
        # input_data e record(unknown()): aceita qualquer objeto sem validar o
        # interior. A chave errada passa com HTTP 200 e falha depois no modelo,
        # com "messages: at least one message is required". A certa e request.
        s.json("POST", f"/api/frames/{fid}/message", json={"input_data": {"request": pedido}})
        f = _espera(s, fid, timeout_s)
        ms = _mensagens(s, fid)

    out: dict[str, Any] = {
        "frame_id": fid,
        "status": f.get("status"),
        "resposta": _resposta_final(ms),
        "mensagens": len(ms),
        **_custos(f),
    }
    if f.get("status") != "completed":
        erro = (f.get("output_data") or {}).get("error")
        out["erro"] = str(erro)[:500] if erro else f.get("status_description")
    _log_provenance(
        "science_ask",
        {"pedido": pedido[:200], "agente": agente, "projeto": pid, "timeout_s": timeout_s},
        out,
    )
    return out


@mcp.tool()
def science_followup(
    frame_id: str,
    pedido: str,
    timeout_s: int = ESPERA_PADRAO_S,
) -> dict[str, Any]:
    """Continua uma conversa existente, preservando o contexto dela.

    Mais leve que science_ask: o prompt de sistema do orquestrador ja esta em
    cache do lado do daemon, entao gasta menos da janela do Max.
    """
    with Sessao(timeout=max(120, timeout_s)) as s:
        _exige_assinatura(s)
        s.json("POST", f"/api/frames/{frame_id}/message", json={"input_data": {"request": pedido}})
        f = _espera(s, frame_id, timeout_s)
        ms = _mensagens(s, frame_id)
    out = {
        "frame_id": frame_id,
        "status": f.get("status"),
        "resposta": _resposta_final(ms),
        "mensagens": len(ms),
        **_custos(f),
    }
    _log_provenance("science_followup", {"frame_id": frame_id, "pedido": pedido[:200]}, out)
    return out


@mcp.tool()
def science_frames(limite: int = 15, so_ativos: bool = False) -> list[dict[str, Any]]:
    """Conversas recentes, com status e custo. Use para achar um frame_id."""
    with Sessao(timeout=30) as s:
        fr = s.json("GET", "/api/frames") or []
    if so_ativos:
        fr = [f for f in fr if f.get("status") not in TERMINAIS]
    fr = sorted(fr, key=lambda f: f.get("updated_at") or "", reverse=True)[: max(1, limite)]
    return [
        {
            "frame_id": f.get("root_frame_id") or f.get("id"),
            "nome": f.get("name"),
            "agente": f.get("agent_name"),
            "status": f.get("status"),
            "mensagens": f.get("message_count"),
            "custo_usd": f.get("total_cost"),
            "atualizado": f.get("updated_at"),
        }
        for f in fr
    ]


@mcp.tool()
def science_read(
    frame_id: str,
    ultimas: int = 6,
    chars_por_mensagem: int = 2000,
) -> dict[str, Any]:
    """Le uma conversa sem gastar nada.

    Truncado de proposito: uma conversa do orquestrador carrega transcricao de
    ferramentas e passa de 100 mil tokens. Suba `chars_por_mensagem` so quando
    precisar do detalhe de uma mensagem especifica.
    """
    with Sessao(timeout=60) as s:
        f = _frame(s, frame_id)
        ms = _mensagens(s, frame_id)
    recorte = ms[-max(1, ultimas) :]
    return {
        "frame_id": frame_id,
        "nome": f.get("name"),
        "status": f.get("status"),
        "mensagens_total": len(ms),
        **_custos(f),
        "mensagens": [
            {"papel": m.get("role"), "texto": _texto(m.get("content"))[:chars_por_mensagem]}
            for m in recorte
        ],
    }


@mcp.tool()
def science_cancel(frame_id: str) -> dict[str, Any]:
    """Interrompe uma conversa em andamento, parando o gasto."""
    with Sessao(timeout=60) as s:
        s.pedir("POST", f"/api/frames/{frame_id}/cancel")
        f = _espera(s, frame_id, 60)
    out = {"frame_id": frame_id, "status": f.get("status"), **_custos(f)}
    _log_provenance("science_cancel", {"frame_id": frame_id}, out)
    return out


if __name__ == "__main__":
    mcp.run()
