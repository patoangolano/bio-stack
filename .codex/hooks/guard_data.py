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


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0  # entrada malformada nunca deve travar a sessao

    if payload.get("tool_name") not in ("Bash", "PowerShell"):
        return 0

    tool_input = payload.get("tool_input") or {}
    command = str(tool_input.get("command", ""))
    if not command:
        return 0

    alvo = _overwrites_existing(command)
    if alvo:
        print(
            "[guard_data] comando bloqueado pela bio-stack.\n"
            f"  motivo: o redirecionamento sobrescreveria {alvo}, que ja existe.\n"
            "  Escreva num nome novo e compare antes de substituir.",
            file=sys.stderr,
        )
        return 2

    for pattern, reason in RULES:
        if pattern.search(command):
            print(
                "[guard_data] comando bloqueado pela bio-stack.\n"
                f"  motivo: {reason}\n"
                f"  comando: {command[:300]}\n"
                "Se voce tem certeza e o usuario confirmou explicitamente, "
                "peca para ele rodar manualmente no terminal dele.",
                file=sys.stderr,
            )
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
