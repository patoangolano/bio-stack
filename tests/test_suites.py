"""Expoe ao pytest as suites autonomas `tests/check_*.py`.

    uv run pytest -q

As suites nao sao testes de pytest: cada `check_*.py` e um script PEP 723 com
metadados de dependencia inline, roda por `uv run --script`, imprime o proprio
relatorio e devolve o veredito no codigo de saida. Esse formato existe por um
motivo — `check_determinism.py` precisa de scanpy 1.11.5 e `check_cache.py` do
pacote mcp, versoes que nao sao as do .venv do projeto, e o uv monta um
ambiente por script.

O que faltava era a ponte: sem ela, `uv run pytest` coletava os arquivos pelo
prefixo `test_`, importava cada um, e o `sys.exit()` no nivel do modulo subia
como SystemExit e matava a sessao inteira com INTERNALERROR. O prefixo
`check_` tira os scripts da coleta; este arquivo os devolve como um teste cada,
executado como subprocesso, que e exatamente como eles foram escritos para
rodar.

`smoke_design.py` fica de fora de proposito: faz chamadas de rede e gasta
credito da NVIDIA.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
ROOT = TESTS.parent

#: Teto generoso: `check_determinism.py` resolve e instala scanpy num ambiente
#: proprio antes de rodar PCA, UMAP e Leiden. Na primeira execucao, com o cache
#: do uv frio, isso passa de dois minutos.
TIMEOUT_S = 900

SUITES = sorted(p.name for p in TESTS.glob("check_*.py"))

#: Vazio quando o uv nao esta no PATH; o skipif abaixo cuida do caso.
UV = shutil.which("uv") or ""


@pytest.mark.skipif(not UV, reason="uv fora do PATH; cada suite resolve as proprias deps")
@pytest.mark.parametrize("suite", SUITES, ids=[s[:-3] for s in SUITES])
def test_suite_autonoma(suite: str) -> None:
    # VIRTUAL_ENV aponta para o .venv do projeto quando o pytest roda sob
    # `uv run`. Deixa-lo aqui faz o uv avisar que o ambiente do script nao bate
    # com ele a cada chamada; o script tem o proprio ambiente, por design.
    env = {k: v for k, v in os.environ.items() if k != "VIRTUAL_ENV"}

    proc = subprocess.run(
        [UV, "run", "--script", str(TESTS / suite)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        cwd=str(ROOT),
        env=env,
        timeout=TIMEOUT_S,
        check=False,
    )
    # O relatorio da suite e a informacao util quando ela falha: cada linha diz
    # qual asercao caiu. Vai para o stdout capturado do pytest, que so e
    # impresso em caso de falha.
    print(proc.stdout)
    if proc.stderr:
        print(proc.stderr)
    assert proc.returncode == 0, f"{suite} saiu com {proc.returncode}"
