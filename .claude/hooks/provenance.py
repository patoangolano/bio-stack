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
    if tool_name.startswith(IGNORED_PREFIXES):
        return 0

    tool_input = payload.get("tool_input") or {}
    response = payload.get("tool_response")

    entry = {
        "ts": datetime.now(UTC).isoformat(timespec="seconds"),
        "session": payload.get("session_id"),
        "kind": classify(tool_name, tool_input),
        "tool": tool_name,
        "input": redact(tool_input),
    }

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
