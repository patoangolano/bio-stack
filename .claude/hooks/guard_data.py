# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""
PreToolUse hook - protege dados de sequenciamento contra comandos destrutivos.

Motivacao pratica, nao burocratica: um FASTQ bruto de um experimento seu pode
ser insubstituivel. Um comando de remocao mal formado apaga semanas de bancada.
Este hook intercepta comandos de shell antes de executarem e bloqueia os que
apagariam dados brutos ou reescreveriam historico do git.

Cobre POSIX e PowerShell. Isso importa: a maquina alvo e Windows, e
`Remove-Item -Recurse -Force data\\raw` e a forma que um agente escreve
naturalmente aqui - nenhuma regra escrita so para `rm` a pegaria.

Protocolo de hook do Claude Code:
  - entrada: JSON no stdin com tool_name e tool_input
  - saida:   exit 0 libera; exit 2 bloqueia e devolve stderr ao modelo
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

# Sem isto, acentos na mensagem de bloqueio saem corrompidos no console do
# Windows (cp1252). A saida tambem e mantida em ASCII puro por seguranca.
try:
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except (AttributeError, OSError):
    pass

# Diretorios e extensoes tratados como dados brutos irrecuperaveis. Ancorados
# em limite de caminho para nao casarem dentro de nomes maiores - "bam" nao
# pode bloquear "bambu_resultados/".
PROTECTED = (
    r"(?:^|[\s/\\'\"=])"
    r"(?:data[/\\]raw|raw|fastq|fq|bam|cram|sra)"
    r"(?:[/\\]|[\s'\"]|$)"
)

# Verbos de remocao do PowerShell. `rm` ja e alias de Remove-Item no PS, entao
# a regra POSIX cobre esse caso - mas so ele.
PS_REMOVE = r"(?:Remove-Item|Remove-ItemProperty|Clear-Content|\bri\b)"

# (regex, explicacao). A explicacao volta para o modelo, entao precisa dizer
# o que fazer em vez de apenas o que nao fazer.
RULES: list[tuple[re.Pattern[str], str]] = [
    (
        re.compile(r"\brm\s+(?:-[a-zA-Z]+\s+)*\S*" + PROTECTED, re.I),
        "remocao dentro de diretorio de dados brutos. Se a intencao e liberar "
        "espaco, mova para uma pasta de descarte e confirme com o usuario antes.",
    ),
    (
        re.compile(PS_REMOVE + r"[^|;]*" + PROTECTED, re.I),
        "remocao de dados brutos via PowerShell. Mesma regra: mova em vez de "
        "apagar, e confirme com o usuario.",
    ),
    (
        re.compile(PROTECTED + r"[^|;]*\|\s*(?:" + PS_REMOVE + r"|rm)\b", re.I),
        "pipeline encanando dados brutos para um comando de remocao.",
    ),
    (
        re.compile(r"\bfind\b[^|;]*" + PROTECTED + r"[^|;]*(?:-delete|-exec\s+rm)", re.I),
        "find com -delete ou -exec rm sobre dados brutos.",
    ),
    (
        re.compile(
            r"\brm\s+-[a-zA-Z]*[rf][a-zA-Z]*\s+(?:-[a-zA-Z-]+\s+)*"
            r"(?:/|~|\*|\$HOME)(?:\s|/|\*|$)"
        ),
        "remocao recursiva na raiz, no home ou com curinga solto. "
        "Quase sempre um erro de digitacao.",
    ),
    (
        re.compile(r"--no-preserve-root"),
        "--no-preserve-root desarma a unica protecao que o rm tem. "
        "Nunca e necessario em analise de dados.",
    ),
    (
        re.compile(r"\bgit\s+push\b.*(--force(?!-with-lease)|-f\b)"),
        "push forcado reescreve historico compartilhado. Use --force-with-lease, "
        "e so depois de confirmar com quem mais usa o repositorio.",
    ),
    (
        re.compile(r"\bgit\s+(reset\s+--hard|clean\s+-[a-zA-Z]*f)"),
        "descarta alteracoes nao commitadas de forma irreversivel. "
        "Faca `git stash` primeiro se houver qualquer duvida.",
    ),
    (
        re.compile(r"\b(mkfs|dd\s+if=.*\s+of=/dev/|fdisk|diskpart|Format-Volume)\b", re.I),
        "operacao de disco de baixo nivel. Fora do escopo de analise de dados.",
    ),
    (
        re.compile(
            r"\b(curl|wget|iwr|Invoke-WebRequest)\b.*\|\s*"
            r"(bash|sh|pwsh|powershell|iex|Invoke-Expression)",
            re.I,
        ),
        "baixar e executar script direto do pipe. Baixe o arquivo, leia, e so entao rode.",
    ),
]

DATA_EXT = (
    ".fastq", ".fq", ".bam", ".cram", ".vcf", ".bcf",
    ".fastq.gz", ".fq.gz", ".vcf.gz", ".bam.bai", ".sam",
)


def _overwrites_existing(command: str) -> str | None:
    """
    Redirecionamento `>` so destroi se o alvo JA existe.

    A versao anterior desta regra bloqueava qualquer `comando > saida.bam`,
    inclusive a criacao de arquivo novo - que e o idioma canonico de
    bioinformatica e acontece em praticamente toda sessao. Checar a existencia
    e o que separa "sobrescrever resultado" de "gerar resultado".

    `>>` (append) nunca destroi, entao e ignorado.
    """
    for m in re.finditer(r"(?<!>)>(?!>)\s*([\"']?)([^\s\"'|;&]+)\1", command):
        alvo = m.group(2)
        if not alvo.lower().endswith(DATA_EXT):
            continue
        try:
            if Path(alvo).exists():
                return alvo
        except OSError:
            continue
    return None


#: Ferramentas que escrevem arquivo sem passar por shell. O guard original so
#: olhava Bash/PowerShell, entao `Write` sobre um FASTQ em data/raw passava sem
#: nenhuma checagem — o caminho mais curto para destruir dado bruto era tambem
#: o unico desprotegido.
WRITE_TOOLS = ("Write", "Edit", "MultiEdit", "NotebookEdit")

#: Ancorado em limite de caminho: `data/raw` casa, `metadata/rawstats` nao.
RAW_DIR = re.compile(r"(?:^|[/\\])data[/\\]raw(?:[/\\]|$)", re.I)


def _guard_write(tool_name: str, tool_input: dict) -> str | None:
    """
    Bloqueia escrita dentro de data/raw e sobrescrita de arquivo de dado bruto.

    Criar arquivo novo em qualquer outro lugar e livre — o objetivo e proteger
    o insubstituivel, nao atrapalhar a analise.
    """
    alvo = str(
        tool_input.get("file_path")
        or tool_input.get("notebook_path")
        or tool_input.get("path")
        or ""
    )
    if not alvo:
        return None

    if RAW_DIR.search(alvo):
        return (
            f"{alvo} esta dentro de data/raw, que guarda dado bruto irrecuperavel.\n"
            "  Escreva em data/proc/ ou runs/ e deixe o bruto intacto. Se o arquivo "
            "precisa mesmo entrar no acervo bruto, o usuario copia manualmente."
        )

    if alvo.lower().endswith(DATA_EXT):
        try:
            existe = Path(alvo).exists()
        except OSError:
            existe = False
        if existe:
            return (
                f"{alvo} ja existe e e um arquivo de dado bruto.\n"
                "  Sobrescrever apaga o original sem historico. Escreva num nome "
                "novo e compare antes de substituir."
            )
    return None


def _finish(is_agy: bool, blocked: bool, reason: str = "") -> int:
    if is_agy:
        if blocked:
            print(json.dumps({"decision": "deny", "reason": reason}))
        else:
            print(json.dumps({"decision": "allow"}))
        return 0
    if blocked:
        print(reason, file=sys.stderr)
        return 2
    return 0


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # entrada malformada nunca deve travar a sessao

    is_agy = "hookEventName" in payload or "toolCall" in payload
    if is_agy:
        tc = payload.get("toolCall") or {}
        raw_name = tc.get("name", "")
        raw_args = tc.get("args") or {}
        if raw_name in ("write_file", "edit_file"):
            tool_name = "Write"
            tool_input = {"file_path": raw_args.get("Path") or raw_args.get("file_path") or ""}
        elif raw_name == "run_command":
            tool_name = "PowerShell"
            tool_input = {"command": raw_args.get("CommandLine") or raw_args.get("command") or ""}
        else:
            tool_name = raw_name
            tool_input = raw_args
    else:
        tool_name = payload.get("tool_name")
        tool_input = payload.get("tool_input") or {}

    if tool_name in WRITE_TOOLS:
        motivo = _guard_write(tool_name, tool_input)
        if motivo:
            return _finish(
                is_agy,
                True,
                "[guard_data] escrita bloqueada pela bio-stack.\n"
                f"  motivo: {motivo}",
            )
        return _finish(is_agy, False)

    if tool_name not in ("Bash", "PowerShell"):
        return _finish(is_agy, False)

    command = str(tool_input.get("command", ""))
    if not command:
        return _finish(is_agy, False)

    alvo = _overwrites_existing(command)
    if alvo:
        return _finish(
            is_agy,
            True,
            "[guard_data] comando bloqueado pela bio-stack.\n"
            f"  motivo: o redirecionamento sobrescreveria {alvo}, que ja existe.\n"
            "  Escreva num nome novo e compare antes de substituir.",
        )

    for pattern, reason in RULES:
        if pattern.search(command):
            return _finish(
                is_agy,
                True,
                "[guard_data] comando bloqueado pela bio-stack.\n"
                f"  motivo: {reason}\n"
                f"  comando: {command[:300]}\n"
                "Se voce tem certeza e o usuario confirmou explicitamente, "
                "peca para ele rodar manualmente no terminal dele.",
            )

    return _finish(is_agy, False)


if __name__ == "__main__":
    sys.exit(main())
