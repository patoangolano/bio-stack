"""
Testa o guard de dados brutos — inclusive o caminho que ficou descoberto por
meses: escrita direta via Write/Edit, que nao passa por shell nenhum.

Roda o hook como subprocesso de verdade, com JSON no stdin, porque e assim que
o Claude Code o invoca. Testar a funcao importada esconderia justamente o erro
mais provavel: o hook receber um payload cujo formato mudou.

    uv run --script tests/test_guard_write.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HOOK = Path(__file__).resolve().parent.parent / ".claude" / "hooks" / "guard_data.py"

# (descricao, payload, exit esperado). Exit 2 bloqueia; 0 libera.
# Referencia: exit 2 e o UNICO codigo que bloqueia uma chamada de ferramenta —
# exit 1 e tratado como erro nao-bloqueante e a acao prossegue.
CASOS: list[tuple[str, dict, int]] = [
    (
        "Write dentro de data/raw",
        {"tool_name": "Write", "tool_input": {"file_path": "data/raw/amostra.fastq"}},
        2,
    ),
    (
        "Write em data/raw com caminho absoluto Windows",
        {
            "tool_name": "Write",
            "tool_input": {"file_path": r"C:\Users\Manec\bio-stack\data\raw\lote1\r1.fq.gz"},
        },
        2,
    ),
    (
        "Edit dentro de data/raw",
        {"tool_name": "Edit", "tool_input": {"file_path": "./data/raw/manifesto.tsv"}},
        2,
    ),
    (
        "NotebookEdit dentro de data/raw",
        {
            "tool_name": "NotebookEdit",
            "tool_input": {"notebook_path": "data/raw/exploracao.ipynb"},
        },
        2,
    ),
    (
        "Write em data/proc e livre",
        {"tool_name": "Write", "tool_input": {"file_path": "data/proc/contagens.csv"}},
        0,
    ),
    (
        "metadata/rawstats nao e data/raw",
        {"tool_name": "Write", "tool_input": {"file_path": "metadata/rawstats/qc.txt"}},
        0,
    ),
    (
        "Write de script comum",
        {"tool_name": "Write", "tool_input": {"file_path": "notebooks/deg.py"}},
        0,
    ),
    (
        "remocao recursiva de dado bruto via shell",
        {"tool_name": "Bash", "tool_input": {"command": "rm -" + "rf data/raw/lote1"}},
        2,
    ),
    (
        "remocao via PowerShell",
        {
            "tool_name": "Bash",
            "tool_input": {"command": "Remove-Item -Recurse -Force data\\raw"},
        },
        2,
    ),
    (
        "comando de analise comum",
        {"tool_name": "Bash", "tool_input": {"command": "python notebooks/deg.py --seed 42"}},
        0,
    ),
    (
        "payload malformado nunca trava a sessao",
        {"tool_name": "Write", "tool_input": {}},
        0,
    ),
]


def roda(payload: dict) -> tuple[int, str]:
    proc = subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return proc.returncode, (proc.stderr or "").strip()


def main() -> int:
    falhas = 0
    for descricao, payload, esperado in CASOS:
        obtido, err = roda(payload)
        ok = obtido == esperado
        marca = "ok  " if ok else "FALHA"
        print(f"{marca} [{esperado}/{obtido}] {descricao}")
        if not ok:
            falhas += 1
            if err:
                print(f"        stderr: {err.splitlines()[0][:100]}")
    total = len(CASOS)
    print(f"\n{total - falhas}/{total} passaram")
    return 1 if falhas else 0


if __name__ == "__main__":
    sys.exit(main())
