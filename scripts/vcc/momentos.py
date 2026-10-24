"""Momentos do pool de controle de um contexto: o denominador de todo o gerador.

Uma varredura por contexto, em blocos de linhas, sem materializar a matriz
(104,5e6 nao nulos em D). Produz, em float64:

    P[g]  soma de contagens por gene
    L_i   biblioteca de cada celula
    A[g]  media ARITMETICA do CPM por celula -- `(1/n) sum_i 1e6*x_ig/L_i`
    pi[g] P[g]/S           composicao do pool
    b[g]  log1p(50_000*pi) == `bulk_lognorm` do pool inteiro, bit a bit
    gate  A > 5.0, ESTRITAMENTE

O gate e a parte que mais parece detalhe e nao e. O cell-eval2 filtra genes em
`de_compute.py:1257` por `genes[ref_cpm > threshold]`, com `ref_cpm` vindo da
media ARITMETICA de CPM por celula e `cpm_factor` exatamente 1,0 na competicao
(`de_compute.py:487`). Nao e o CPM do pseudobulk: medido, os dois conjuntos
diferem em 82 genes em D, 29 em E e 48 em F. Sao os genes do gate -- e so eles
-- que `de_sig_jaccard`, `de_lfc_nmae` e as duas de direcao leem.

O recast para float64 antes de qualquer soma e obrigatorio, nao higiene: o
`.data` do h5ad e float32, e somar 104e6 parcelas em float32 acumula erro. `A`
e o denominador do gate e do log2FC.

`sum_g A[g] == 1e6` e identidade, nao medida: cada celula contribui exatamente
1e6 de CPM, logo a media entre celulas tambem. Nao existe "fator global" a
corrigir a jusante.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pseudobulk import ALVO_SOMA

PISO_CPM_DO_GATE: float = 5.0  # de_compute.py:1257, estritamente maior
LINHAS_POR_BLOCO: int = 2048
CPM: float = 1e6


@dataclass(frozen=True)
class MomentosControle:
    contexto: str
    genes: list[str]
    n_celulas: int
    soma_por_gene: np.ndarray  # P [n_genes] float64
    soma_total: float  # S
    bibliotecas: np.ndarray  # L [n_celulas] float64
    lote: np.ndarray  # codigos de ntc_id [n_celulas] int16
    cpm_medio_por_celula: np.ndarray  # A [n_genes] float64
    pi: np.ndarray  # [n_genes] float64
    w: np.ndarray  # ALVO_SOMA * pi
    bulk: np.ndarray  # b = log1p(w)
    gate: np.ndarray  # [n_genes] bool
    alcancavel: np.ndarray  # pi > 0
    zerofrac: np.ndarray  # [n_genes] float64
    nnz: int
    contagem_maxima: float


def _abrir(caminho: Path) -> tuple[h5py.File, h5py.Group, pd.DataFrame, list[str]]:
    f = h5py.File(caminho, "r")
    try:
        x, g_obs, g_var = f["X"], f["obs"], f["var"]
        if not isinstance(x, h5py.Group):
            raise ValueError(f"{caminho}: X nao e esparso")
        assert isinstance(g_obs, h5py.Group) and isinstance(g_var, h5py.Group)
        obs, var = ad.io.read_elem(g_obs), ad.io.read_elem(g_var)
        assert isinstance(obs, pd.DataFrame) and isinstance(var, pd.DataFrame)
        return f, x, obs, var.index.astype(str).tolist()
    except BaseException:
        f.close()
        raise


def _blocos(x: h5py.Group, n_linhas: int, linhas_por_bloco: int, n_genes: int):
    """Produz (r0, bloco CSR float64) percorrendo as linhas sem materializar nada."""
    dados, indices = x["data"], x["indices"]
    assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
    indptr = np.asarray(x["indptr"]).astype(np.int64)
    for r0 in range(0, n_linhas, linhas_por_bloco):
        r1 = min(r0 + linhas_por_bloco, n_linhas)
        a, b = int(indptr[r0]), int(indptr[r1])
        yield (
            r0,
            sp.csr_matrix(
                (
                    np.asarray(dados[a:b], dtype=np.float64),
                    np.asarray(indices[a:b]),
                    indptr[r0 : r1 + 1] - a,
                ),
                shape=(r1 - r0, n_genes),
            ),
        )


def momentos_do_controle(
    caminho: Path,
    col_contexto: str = "context",
    col_lote: str = "ntc_id",
    linhas_por_bloco: int = LINHAS_POR_BLOCO,
) -> MomentosControle:
    f, x, obs, genes = _abrir(Path(caminho))
    with f:
        n_genes = len(genes)
        n_celulas = len(obs)
        contexto = str(obs[col_contexto].astype(str).iloc[0]) if col_contexto in obs else ""
        if col_lote in obs:
            codigos, _ = pd.factorize(obs[col_lote].astype(str).to_numpy(), sort=True)
            lote = codigos.astype(np.int16)
        else:
            lote = np.zeros(n_celulas, dtype=np.int16)

        soma_por_gene = np.zeros(n_genes, dtype=np.float64)
        cpm_acumulado = np.zeros(n_genes, dtype=np.float64)
        nao_nulos_por_gene = np.zeros(n_genes, dtype=np.int64)
        bibliotecas = np.zeros(n_celulas, dtype=np.float64)
        nnz = 0
        contagem_maxima = 0.0

        for r0, bloco in _blocos(x, n_celulas, linhas_por_bloco, n_genes):
            lib = np.asarray(bloco.sum(axis=1), dtype=np.float64).ravel()
            bibliotecas[r0 : r0 + lib.size] = lib
            soma_por_gene += np.asarray(bloco.sum(axis=0), dtype=np.float64).ravel()
            nao_nulos_por_gene += bloco.getnnz(axis=0)
            nnz += int(bloco.nnz)
            if bloco.nnz:
                contagem_maxima = max(contagem_maxima, float(bloco.data.max()))
            # `A` precisa da soma de x_ig/L_i: escala cada LINHA por 1/L_i antes
            # de somar nas celulas. Celula de biblioteca nula nao contribui.
            escala = np.zeros_like(lib)
            np.divide(CPM, lib, out=escala, where=lib > 0)
            escalado = sp.csr_matrix(sp.diags(escala) @ bloco)
            cpm_acumulado += np.asarray(escalado.sum(axis=0), dtype=np.float64).ravel()

    soma_total = float(soma_por_gene.sum())
    pi = soma_por_gene / soma_total if soma_total > 0 else np.zeros_like(soma_por_gene)
    # `w` e `P * (ALVO_SOMA/S)`, nao `ALVO_SOMA * pi`: e a ORDEM de operacoes do
    # `bulk_lognorm` (uma escala por linha, depois a multiplicacao), e so com ela
    # `b` sai bit a bit igual ao que o cell-eval2 calcula do mesmo pool. As duas
    # formas diferem por um ulp, o que ja basta para `np.array_equal` recusar --
    # e esse arredondamento e justamente o que faz o delta nulo nao ser identidade.
    escala = ALVO_SOMA / soma_total if soma_total > 0 else 0.0
    w = soma_por_gene * escala
    cpm_medio = cpm_acumulado / n_celulas
    return MomentosControle(
        contexto=contexto,
        genes=genes,
        n_celulas=n_celulas,
        soma_por_gene=soma_por_gene,
        soma_total=soma_total,
        bibliotecas=bibliotecas,
        lote=lote,
        cpm_medio_por_celula=cpm_medio,
        pi=pi,
        w=w,
        bulk=np.log1p(w),
        gate=cpm_medio > PISO_CPM_DO_GATE,
        alcancavel=pi > 0,
        zerofrac=1.0 - nao_nulos_por_gene / n_celulas,
        nnz=nnz,
        contagem_maxima=contagem_maxima,
    )


def carregar_pool(caminho: Path) -> sp.csr_matrix:
    """CSR [n_celulas, n_genes] com `data` uint16 e `indices` int32.

    uint16 porque a contagem maxima medida nos tres contextos e 4.104, e porque
    o pool inteiro de D em float32 custaria 1,2 GB contra 598 MB assim. int32
    nos indices porque o indice de gene vai no maximo a 18.532 -- mas o
    `indptr` fica em int64, que e o dtype que a submissao inteira exige.
    """
    f, x, obs, genes = _abrir(Path(caminho))
    with f:
        dados, indices, indptr = x["data"], x["indices"], x["indptr"]
        assert isinstance(dados, h5py.Dataset) and isinstance(indices, h5py.Dataset)
        assert isinstance(indptr, h5py.Dataset)
        valores = np.asarray(dados[:], dtype=np.float64)
        if not np.array_equal(valores, np.floor(valores)):
            raise ValueError(f"{caminho}: contagens fracionarias no pool de controle")
        if valores.size and (valores.min() < 0 or valores.max() > np.iinfo(np.uint16).max):
            raise ValueError(f"{caminho}: contagem fora da faixa de uint16")
        m = sp.csr_matrix(
            (
                valores.astype(np.uint16),
                np.asarray(indices[:], dtype=np.int32),
                np.asarray(indptr[:], dtype=np.int64),
            ),
            shape=(len(obs), len(genes)),
        )
    m.sort_indices()
    return m
