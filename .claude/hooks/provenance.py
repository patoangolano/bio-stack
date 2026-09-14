# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
PostToolUse hook - caderno de laboratorio automatico.

Toda analise que voce rodar fica registrada em runs/provenance.jsonl: qual
ferramenta, quais parametros, qual arquivo saiu, em que horario UTC. Quando
chegar a hora de escrever a secao de Metodos do relatorio ou do TCC, o historico
ja esta la - inclusive as tentativas que nao deram certo, que e justamente o que
ninguem lembra depois.

Isso substitui o habito de anotar comandos a mao, que todo mundo abandona na
terceira semana.

Nao registra: leituras (Read, Grep, Glob) e nada que contenha chave de API.
"""

from __future__ import annotations

import json
import os
import re
import sys
from datetime import UTC, datetime
from pathlib import Path


def _log_path() -> Path:
    """
    Um valor com "${" significa que a expansao nao aconteceu (o bloco env de
    settings.json nao expande). Cair para um caminho ancorado no proprio hook
    e deterministico, em vez de criar um diretorio chamado "${...}".
    """
    raw = os.environ.get("BIO_STACK_PROVENANCE", "")
    if raw and "${" not in raw:
        return Path(raw)
    return Path(__file__).resolve().parents[2] / "runs" / "provenance.jsonl"


LOG = _log_path()

#: Versao do formato de cada linha. O arquivo ja mudou de forma uma vez sem
#: aviso (linhas antigas nao tem `ok` nem `commit`), e um leitor que assuma o
#: formato novo quebra ou, pior, mente sobre as linhas velhas.
SCHEMA = 3


def _commit() -> str | None:
    """
    SHA do HEAD, lido direto do diretorio .git.

    Poderia ser `git rev-parse HEAD`, mas isso e um subprocesso por chamada de
    ferramenta — em 900 operacoes, minutos de overhead. Ler o arquivo custa
    microssegundos e da a mesma resposta.

    Isto e o que amarra resultado a codigo. Sem o SHA, o log responde "o que eu
    rodei" mas nao "com qual versao do script" — e a segunda pergunta e a que
    torna o resultado reproduzivel.
    """
    raiz = Path(__file__).resolve().parents[2] / ".git"
    try:
        head = (raiz / "HEAD").read_text(encoding="utf-8").strip()
    except OSError:
        return None
    if not head.startswith("ref:"):
        return head[:40] or None  # HEAD destacado ja e o proprio SHA

    ref = head[4:].strip()
    try:
        return (raiz / ref).read_text(encoding="utf-8").strip()[:40]
    except OSError:
        pass
    # Ref empacotada: depois de `git gc` o arquivo solto deixa de existir.
    try:
        for linha in (raiz / "packed-refs").read_text(encoding="utf-8").splitlines():
            if linha.endswith(f" {ref}"):
                return linha.split(" ", 1)[0][:40]
    except OSError:
        pass
    return None


COMMIT = _commit()

# Ferramentas de leitura pura nao mudam nada; registra-las so gera ruido.
IGNORED = {"Read", "Grep", "Glob", "TodoWrite", "ListAgents", "ToolSearch", "WebSearch"}

# Os servidores nvidia-bio e nvidia-nemo ja gravam a propria proveniencia, com
# mais detalhe (parametros do modelo, arquivos gerados). Registrar aqui tambem
# duplicaria cada operacao e inflaria as contagens que /metodos usa.
IGNORED_PREFIXES = ("mcp__nvidia-bio__", "mcp__nvidia-nemo__")

# runs/provenance.jsonl e versionado no git de proposito (e o caderno de
# laboratorio), entao um segredo que escape aqui entra no historico para sempre.
# Vale cobrir os formatos amplamente.
SECRET = re.compile(
    r"(nvapi-[A-Za-z0-9_\-]{10,}"
    r"|sk-[A-Za-z0-9_\-]{20,}"
    r"|gh[oprsu]_[A-Za-z0-9]{20,}"
    r"|github_pat_[A-Za-z0-9_]{22,}"
    r"|hf_[A-Za-z0-9]{30,}"
    r"|xox[baprs]-[A-Za-z0-9-]{10,}"
    r"|AKIA[0-9A-Z]{16}"
    r"|-----BEGIN[ A-Z]*PRIVATE KEY-----"
    r"|(?i:api[_-]?key|token|password|passwd|secret|bearer)[\"']?\s*[=:\s]\s*[\"']?[A-Za-z0-9_\-./+]{8,})"
)

MAX_FIELD = 2000


def redact(value):
    """Remove segredos e trunca campos gigantes antes de gravar."""
    if isinstance(value, str):
        cleaned = SECRET.sub("[REDACTED]", value)
        if len(cleaned) > MAX_FIELD:
            return cleaned[:MAX_FIELD] + f"... [+{len(cleaned) - MAX_FIELD} chars]"
        return cleaned
    if isinstance(value, dict):
        return {k: redact(v) for k, v in value.items()}
    if isinstance(value, list):
        return [redact(v) for v in value[:50]]
    return value


def classify(tool_name: str, tool_input: dict) -> str:
    """Rotula a entrada para facilitar filtrar o log depois."""
    if tool_name.startswith("mcp__nvidia-bio"):
        return "nvidia-nim"
    if tool_name.startswith(("mcp__gget", "mcp__biocontext", "mcp__biomcp")):
        return "database-query"
    if tool_name in ("Write", "Edit", "NotebookEdit"):
        return "code-edit"
    if tool_name in ("Bash", "PowerShell"):
        cmd = str(tool_input.get("command", "")).lower()
        bio = ("nextflow", "samtools", "bcftools", "bwa", "blast", "seqkit",
               "fastqc", "salmon", "star", "hisat", "minimap", "scanpy", "python")
        if any(t in cmd for t in bio):
            return "analysis-run"
        return "shell"
    return "other"


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    tool_name = payload.get("tool_name", "")
    if not tool_name or tool_name in IGNORED:
        return 0

    event = payload.get("hook_event_name") or "PostToolUse"
    falhou = event == "PostToolUseFailure"

    # A excecao dos servidores proprios vale so para chamada bem-sucedida.
    # Em `servers/nvidia_bio_mcp/server.py` todo `_log_provenance(...)` roda
    # DEPOIS de obter o resultado, entao quando o POST levanta excecao — o 422
    # que `evo2_score` dava em toda chamada, por exemplo — a linha de log nunca
    # e alcancada. Sem esta ressalva a falha nao aparece em lugar nenhum, que e
    # como tres rotas quebradas passaram meses sem ninguem notar.
    if not falhou and tool_name.startswith(IGNORED_PREFIXES):
        return 0

    tool_input = payload.get("tool_input") or {}
    response = payload.get("tool_response")

    # Qual evento disparou. Sem isto o log so tinha sucesso: `PostToolUse` nao
    # dispara quando a ferramenta falha, entao uma predicao que deu 422 e um
    # script que estourou excecao sumiam do caderno. Uma secao de Metodos
    # escrita a partir de um log so-de-sucesso descreve um experimento que nao
    # aconteceu — e o vies e invisivel, porque nada no arquivo denuncia a falta.
    entry = {
        "v": SCHEMA,
        "ts": datetime.now(UTC).isoformat(timespec="seconds"),
        "session": payload.get("session_id"),
        "commit": COMMIT,
        "kind": classify(tool_name, tool_input),
        "tool": tool_name,
        "ok": not falhou,
        "input": redact(tool_input),
    }

    if falhou:
        motivo = payload.get("tool_error") or payload.get("error") or ""
        if isinstance(motivo, dict):
            motivo = motivo.get("message") or json.dumps(motivo, ensure_ascii=False)
        entry["erro"] = str(motivo)[:400] or "(sem mensagem)"

    # Guarda so o essencial da resposta: caminhos de arquivo e status de erro.
    if isinstance(response, dict):
        slim = {}
        for key in ("files", "filePath", "file_path", "num_models", "confidence"):
            if key in response:
                slim[key] = redact(response[key])
        if response.get("is_error") or response.get("isError"):
            slim["error"] = True
        if slim:
            entry["output"] = slim

    try:
        LOG.parent.mkdir(parents=True, exist_ok=True)
        with LOG.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except OSError:
        pass  # nunca derruba a sessao por causa do log

    return 0


if __name__ == "__main__":
    sys.exit(main())
