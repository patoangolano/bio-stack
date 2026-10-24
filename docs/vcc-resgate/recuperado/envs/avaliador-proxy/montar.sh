#!/usr/bin/env bash
# envs/avaliador-proxy/montar.sh — monta o ambiente isolado do avaliador local
# do Virtual Cell Challenge 2026 (etapa 5 de docs/VCC-2026-ROTA.md).
#
# Rodar do Windows, com `bash -l` (nao `bash`): o PATH do `uv` esta em
# /usr/local/bin, mas o do pipx e dos shims do pixi vem de ~/.profile, lido so
# em shell de login, e este script chama `python3` do sistema para a conferencia.
#
#   wsl -d Ubuntu-24.04 -- bash -l /mnt/c/Users/Manec/bio-stack/envs/avaliador-proxy/montar.sh
#
# Idempotente: reexecutar reinstala as versoes exatas de requirements.in sobre o
# venv existente. Passe --recriar para apagar o venv antes.
set -euo pipefail

PROJETO=/mnt/c/Users/Manec/bio-stack
RAIZ_ENV="$HOME/envs/avaliador-proxy"
VENV="$RAIZ_ENV/.venv"
REQ="$PROJETO/envs/avaliador-proxy/requirements.in"
PY=3.12

# O venv NAO mora em /mnt/c de proposito: sao ~4 500 arquivos pequenos entre
# numba, llvmlite e pyarrow, e cada import atravessaria a fronteira de sistema
# de arquivos. Medido no projeto: ordens de grandeza mais lento.
if [[ "${1:-}" == "--recriar" ]]; then
    rm -rf "$VENV"
fi

if [[ ! -r "$REQ" ]]; then
    echo "erro: $REQ nao encontrado — o projeto esta montado em $PROJETO?" >&2
    exit 1
fi

command -v uv >/dev/null || { echo "erro: uv nao esta no PATH (esperado /usr/local/bin/uv)" >&2; exit 1; }

mkdir -p "$RAIZ_ENV"
[[ -d "$VENV" ]] || uv venv "$VENV" --python "$PY"

# --exact remove o que nao esta no requirements.in: sem isso, uma dependencia
# que saiu do arquivo continua importavel e o ambiente deixa de ser descrito
# pelo que esta versionado.
uv pip install --python "$VENV/bin/python" --exact -r "$REQ"

# O freeze completo fica ao lado, como em envs/bonej-headless: o requirements.in
# declara a intencao, o requirements.txt registra o que entrou de fato,
# incluindo as transitivas que nenhum pin nosso fixa.
uv pip freeze --python "$VENV/bin/python" > "$PROJETO/envs/avaliador-proxy/requirements.txt"

echo
echo "=== conferencia: import real, nao --help ==="
# `--help` do cell-eval2 nao importa polars, nao resolve backend de DE e nao
# toca o gate de CPM. Estes tres toques exercitam as tres coisas.
"$VENV/bin/python" - <<'PY'
import importlib.metadata as md

import numpy as np
import polars as pl
import scipy.sparse as sp

from cell_eval2.config import EvalConfig
from cell_eval2.de_compute import _ref_cpm_from_cells
from cell_eval2.run import _cache_backend, _cache_device

for nome in ("cell-eval2", "pdex", "numba-mwu", "numpy", "scipy", "polars", "pyarrow", "numba"):
    print(f"  {nome:12s} {md.version(nome)}")

cfg = EvalConfig.from_preset("vcc2026")
# Backend fixado: `_cache_backend` devolve a string verbatim, sem sondar driver.
from dataclasses import replace
cfg_pdex = replace(cfg, de=replace(cfg.de, backend="pdex"), device="cpu")
print(f"  preset vcc2026: metrics={cfg.metrics} pert_col={cfg.pert_col!r} "
      f"cache_strict={cfg.cache_strict} filtro_cpm={cfg.filter.filter_gene_min_cpm_cell}")
print(f"  resolvido com backend fixo: device={_cache_device(cfg_pdex)!r} "
      f"backend={_cache_backend(cfg_pdex)!r}")

# O gate de DE do preset: media, por gene, do CPM POR CELULA das celulas de
# referencia. Duas celulas, dois genes, matriz [[1,1],[0,1]]: o CPM por celula
# do gene 0 e 0,5e6 e 0, media 2,5e5; do gene 1 e 0,5e6 e 1e6, media 7,5e5. O
# CPM do PSEUDOBULK daria 3,33e5 e 6,67e5 — numeros diferentes, e o gate usa o
# primeiro par. E a definicao que reproduz o universo de DE oficial.
X = sp.csr_matrix(np.array([[1.0, 1.0], [0.0, 1.0]]))
cpm = _ref_cpm_from_cells(X, n_genes=2)
assert np.allclose(cpm, [2.5e5, 7.5e5]), cpm
print(f"  gate de CPM por celula (referencia): {cpm.tolist()}")

# pdex e numba_mwu importaveis. O import do numba nao compila nada — a primeira
# chamada compila —, entao isto prova o wheel, nao o JIT; o JIT e exercitado de
# verdade na corrida de verificacao, que e o motivo de ela existir.
import numba_mwu  # noqa: F401
import pdex  # noqa: F401
print(f"  numba_mwu importado, pdex importado, polars {pl.__version__}")
print("OK")
PY
echo
echo "venv: $VENV/bin/python"
echo "proximo passo: bash -l $PROJETO/scripts/vcc/avaliador_corrida.sh"
