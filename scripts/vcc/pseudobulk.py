"""Pseudobulk no espaco em que o placar do VCC 2026 compara.

`bulk_lognorm` e porte linha a linha de `cell_eval2.prep.bulk_lognorm_means`
(cell-eval2 0.19.0): soma as contagens do grupo, normaliza a soma para
`ALVO_SOMA` e aplica `log1p` uma vez. NAO e a media de log1p por celula -- essa
e outra coisa, e era o que o `ctrl_means` do PIE carregava. O teste
`test_bulk_lognorm_identico_ao_do_cell_eval2_instalado` confere bit a bit
contra o pacote instalado em `envs/cell-eval2`.

`somas_por_grupo` le o CSR de um `.h5ad` em blocos de linhas, sem materializar a
matriz: a submissao completa tem ~2,1e9 nao nulos.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

# moments.DEFAULT_BULK_TARGET_SUM do cell-eval2 0.19.0.
ALVO_SOMA = 50_000.0
LINHAS_POR_BLOCO = 4096
_SEPARADOR = "\x1f"


def bulk_lognorm(somas: np.ndarray, alvo_soma: float = ALVO_SOMA) -> np.ndarray:
    """`log1p(alvo_soma * P / sum_g P)` por linha de uma matriz [grupos, genes] de SOMAS.

    Grupo todo nulo devolve zeros, como no cell-eval2. A ordem das operacoes e a
    mesma do original, para o resultado sair identico bit a bit.
    """
    somas = np.asarray(somas, dtype=np.float64)
    if somas.ndim != 2:
        raise ValueError(f"somas tem de ser [grupos, genes]; veio {somas.shape}")
    totais = somas.sum(axis=1)
    escala = np.zeros_like(totais)
    ok = totais > 0
    escala[ok] = alvo_soma / totais[ok]
    return np.log1p(somas * escala[:, None])


def somas_por_grupo(
    caminho: Path,
    colunas: Sequence[str],
    linhas_por_bloco: int = LINHAS_POR_BLOCO,
) -> tuple[list[tuple[str, ...]], np.ndarray, list[str]]:
    """Soma das contagens por grupo de celulas, definido pelos valores de `colunas` em obs.

    Devolve (chaves dos grupos em ordem de primeira aparicao, somas [grupos, genes]
    em float64, eixo de genes do arquivo).
    """
    with h5py.File(caminho, "r") as f:
        g_obs, g_var, x = f["obs"], f["var"], f["X"]
        assert isinstance(g_obs, h5py.Group) and isinstance(g_var, h5py.Group)
        if not isinstance(x, h5py.Group):
            raise ValueError(f"{caminho}: X nao e esparso; rode validar_submissao.py primeiro")
        obs, var = ad.io.read_elem(g_obs), ad.io.read_elem(g_var)
        assert isinstance(obs, pd.DataFrame) and isinstance(var, pd.DataFrame)
        genes = var.index.astype(str).tolist()

        rotulo = obs[colunas[0]].astype(str)
        for c in colunas[1:]:
            rotulo = rotulo + _SEPARADOR + obs[c].astype(str)
        codigos, unicos = pd.factorize(rotulo.to_numpy())
        chaves = [tuple(str(u).split(_SEPARADOR)) for u in unicos]

        dados, indices = x["data"], x["indices"]
        assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
        indptr = np.asarray(x["indptr"]).astype(np.int64)
        ng = len(genes)
        somas = np.zeros((len(chaves), ng), dtype=np.float64)
        n_linhas = indptr.size - 1
        for r0 in range(0, n_linhas, linhas_por_bloco):
            r1 = min(r0 + linhas_por_bloco, n_linhas)
            a, b = int(indptr[r0]), int(indptr[r1])
            bloco = sp.csr_matrix(
                (np.asarray(dados[a:b], dtype=np.float64), np.asarray(indices[a:b]), indptr[r0 : r1 + 1] - a),
                shape=(r1 - r0, ng),
            )
            # So os grupos presentes no bloco: a submissao vem ordenada por grupo,
            # entao um bloco de 4.096 linhas toca ~11 dos 900, nao todos.
            presentes, local = np.unique(codigos[r0:r1], return_inverse=True)
            um = sp.csr_matrix(
                (np.ones(r1 - r0), (local, np.arange(r1 - r0))), shape=(presentes.size, r1 - r0)
            )
            somas[presentes] += sp.csr_matrix(um @ bloco).toarray()
    return chaves, somas, genes
