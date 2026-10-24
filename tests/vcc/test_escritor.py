"""O escritor incremental: 900 grupos de 400 celulas num `.h5ad` sem `vstack`.

A submissao tem 360.000 celulas e cerca de 2,0e9 nao nulos -- 15 GB. Nada dela
pode ser materializado, e e justamente por isso que o `indptr` e o lugar mais
facil de errar: anexar o `indptr` de cada grupo CRU da 900 blocos de 401
entradas (360.900 em vez de 360.001) que reiniciam em zero a cada grupo. Isso
reprova os itens 13, 15 e 17 do contrato de uma vez, e o diagnostico aponta para
"matriz corrompida" em vez de "um off-by-one no escritor".

O segundo lugar: gravar `obs`/`var` como datasets h5py a mao. O
`validar_submissao._ler_tabela` devolve DataFrame VAZIO, e caem os itens 2, 3, 5
e 9 a 12 juntos.
"""

from __future__ import annotations

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from escritor import EscritorCSR, conferir_teto_por_celula

NC = 400


def grupo(rng: np.random.Generator, n_genes: int, n_celulas: int = NC) -> sp.csr_matrix:
    denso = rng.poisson(0.4, size=(n_celulas, n_genes)).astype(np.float32)
    denso[:, 0] = rng.integers(1, 9, size=n_celulas)  # garante linha nao nula
    m = sp.csr_matrix(denso)
    m.eliminate_zeros()
    return m


def escrever(caminho, genes, grupos):
    e = EscritorCSR(caminho, genes)
    for matriz, contexto, alvo in grupos:
        e.acrescentar(matriz, contexto, alvo)
    return e.fechar()


def test_escritor_monta_indptr_e_atributos(tmp_path):
    """T14 -- tres grupos de 400 linhas dao UM `indptr` de 1.201 entradas."""
    rng = np.random.default_rng(2)
    genes = [f"G{i:04d}" for i in range(500)]
    partes = [(grupo(rng, 500), "D", a) for a in ("ABCC1", "ACLY", "ADNP")]
    caminho = tmp_path / "sub.h5ad"
    nnz = escrever(caminho, genes, partes)

    with h5py.File(caminho, "r") as f:
        x = f["X"]
        assert isinstance(x, h5py.Group)
        indptr = np.asarray(x["indptr"])
        assert indptr.size == 3 * NC + 1
        assert int(indptr[0]) == 0
        assert np.all(np.diff(indptr) >= 0)
        assert int(indptr[-1]) == nnz
        assert indptr.dtype == np.int64
        indices = x["indices"]
        dados = x["data"]
        assert isinstance(indices, h5py.Dataset) and isinstance(dados, h5py.Dataset)
        assert indices.dtype == np.int32
        assert dados.dtype == np.float32
        assert dados.size == indices.size == nnz
        assert x.attrs["encoding-type"] == "csr_matrix"
        assert x.attrs["encoding-version"] == "0.1.0"
        assert list(np.asarray(x.attrs["shape"])) == [3 * NC, 500]
        obs_g, var_g = f["obs"], f["var"]
        assert isinstance(obs_g, h5py.Group) and isinstance(var_g, h5py.Group)
        obs = ad.io.read_elem(obs_g)
        var = ad.io.read_elem(var_g)
    assert isinstance(obs, pd.DataFrame) and isinstance(var, pd.DataFrame)
    assert len(obs) == 3 * NC
    assert list(obs.columns) == ["target_gene", "context"]
    assert obs["target_gene"].astype(str).tolist()[:2] == ["ABCC1", "ABCC1"]
    assert sorted(set(obs["target_gene"].astype(str))) == ["ABCC1", "ACLY", "ADNP"]
    assert set(obs["context"].astype(str)) == {"D"}
    assert var.index.astype(str).tolist() == genes
    assert obs.index.astype(str).tolist()[0] == "D_ABCC1_000"


def test_o_arquivo_escrito_e_identico_ao_vstack(tmp_path):
    """T14 -- o teste de corretude do escritor: `vstack` dos grupos, bit a bit.

    `vstack` nao cabe em producao (15 GB), mas cabe em tres grupos, e e o unico
    oraculo honesto para o `indptr` montado a mao.
    """
    rng = np.random.default_rng(4)
    genes = [f"G{i:04d}" for i in range(300)]
    matrizes = [grupo(rng, 300) for _ in range(3)]
    caminho = tmp_path / "sub.h5ad"
    escrever(caminho, genes, [(m, "E", a) for m, a in zip(matrizes, "abc", strict=True)])

    esperado = sp.csr_matrix(sp.vstack(matrizes, format="csr"))
    lido = ad.read_h5ad(caminho).X
    assert isinstance(lido, sp.csr_matrix)
    assert np.array_equal(lido.indptr, esperado.indptr)
    assert np.array_equal(lido.indices, esperado.indices)
    assert np.array_equal(lido.data, esperado.data)


def test_indptr_cru_anexado_seria_pego(tmp_path):
    """T14, controle -- o defeito que este escritor existe para nao ter.

    Reproduz o erro a mao e mostra que ele NAO e silencioso: 1.203 entradas em
    vez de 1.201, reiniciando em zero, o que quebra a monotonicidade.
    """
    rng = np.random.default_rng(6)
    matrizes = [grupo(rng, 200, n_celulas=NC) for _ in range(3)]
    cru = np.concatenate([m.indptr for m in matrizes])
    assert cru.size == 3 * (NC + 1)
    assert not np.all(np.diff(cru) >= 0)
    certo = np.concatenate(
        [[0], *[m.indptr[1:] + sum(k.nnz for k in matrizes[:i]) for i, m in enumerate(matrizes)]]
    )
    assert certo.size == 3 * NC + 1
    assert np.all(np.diff(certo) >= 0)
    assert int(certo[-1]) == sum(m.nnz for m in matrizes)


def test_chunk_pequeno_nao_muda_o_arquivo(tmp_path):
    """O tamanho de chunk do h5py e desempenho, nao conteudo."""
    rng = np.random.default_rng(8)
    genes = [f"G{i:04d}" for i in range(150)]
    matrizes = [grupo(rng, 150, n_celulas=50) for _ in range(2)]
    saidas = []
    for chunk in (262_144, 97):
        caminho = tmp_path / f"sub_{chunk}.h5ad"
        e = EscritorCSR(caminho, genes, elementos_por_chunk=chunk)
        for m, a in zip(matrizes, "xy", strict=True):
            e.acrescentar(m, "F", a)
        e.fechar()
        saidas.append(ad.read_h5ad(caminho).X)
    a, b = saidas
    assert isinstance(a, sp.csr_matrix) and isinstance(b, sp.csr_matrix)
    assert np.array_equal(a.indptr, b.indptr)
    assert np.array_equal(a.data, b.data)


def test_escritor_recusa_matriz_fora_do_eixo(tmp_path):
    rng = np.random.default_rng(9)
    genes = [f"G{i:04d}" for i in range(100)]
    e = EscritorCSR(tmp_path / "sub.h5ad", genes)
    with pytest.raises(ValueError, match="eixo de genes"):
        e.acrescentar(grupo(rng, 99, n_celulas=10), "D", "ABCC1")


def test_escritor_recusa_fracionario_e_negativo(tmp_path):
    """O contrato pede contagem inteira nao negativa; o escritor nao deixa passar."""
    genes = [f"G{i:04d}" for i in range(10)]
    e = EscritorCSR(tmp_path / "sub.h5ad", genes)
    fracionario = sp.csr_matrix(np.array([[1.5, 0, 2.0]] + [[1.0, 0, 0]] * 2, dtype=np.float32))
    fracionario.resize(3, 10)
    with pytest.raises(ValueError, match="inteir"):
        e.acrescentar(fracionario, "D", "ABCC1")
    negativo = sp.csr_matrix(np.array([[-1.0, 0, 2.0]] + [[1.0, 0, 0]] * 2, dtype=np.float32))
    negativo.resize(3, 10)
    with pytest.raises(ValueError, match="negativ"):
        e.acrescentar(negativo, "D", "ABCC1")


def test_conferir_teto_por_celula(tmp_path):
    """A conferencia do teto e feita sobre o ARQUIVO, lido em blocos, nao em memoria."""
    rng = np.random.default_rng(10)
    genes = [f"G{i:04d}" for i in range(80)]
    m = grupo(rng, 80, n_celulas=120)
    caminho = tmp_path / "sub.h5ad"
    escrever(caminho, genes, [(m, "D", "ABCC1")])
    maior, acima = conferir_teto_por_celula(caminho, linhas_por_bloco=17)
    esperado = float(np.asarray(m.sum(axis=1)).ravel().max())
    assert maior == pytest.approx(esperado)
    assert acima == 0
    maior_bloco_grande, _ = conferir_teto_por_celula(caminho, linhas_por_bloco=4096)
    assert maior_bloco_grande == maior


def test_eixo_real_de_18533_genes(tmp_path, ref):
    """Um grupo no eixo oficial: o `var` tem de sair com os 18.533 simbolos em ordem."""
    genes = pd.read_csv(ref / "gene_names.csv")["gene_name"].astype(str).tolist()
    assert len(genes) == 18_533
    rng = np.random.default_rng(12)
    linhas, colunas = 40, len(genes)
    denso = np.zeros((linhas, colunas), dtype=np.float32)
    for i in range(linhas):
        onde = rng.choice(colunas, size=500, replace=False)
        denso[i, onde] = rng.integers(1, 20, size=500)
    caminho = tmp_path / "sub.h5ad"
    nnz = escrever(caminho, genes, [(sp.csr_matrix(denso), "D", "ABCC1")])
    assert nnz == linhas * 500
    a = ad.read_h5ad(caminho)
    assert a.shape == (linhas, 18_533)
    assert a.var.index.astype(str).tolist() == genes
    assert a.var.index[0] == "TSPAN6" and a.var.index[-1] == "PRH1"


def test_conferir_teto_ignora_linha_vazia(tmp_path):
    """`np.add.reduceat` REPETE o valor anterior numa linha vazia, e isso mente.

    O contrato proibe celula de soma zero, entao a submissao nunca tem uma -- mas
    `conferir_teto_por_celula` le um arquivo qualquer, e sem zerar as linhas
    vazias ela reportaria a soma da linha anterior como se fosse delas. O teste
    de mutacao pegou esta lacuna: nenhuma fixture tinha linha vazia, e remover o
    `np.where` passava.
    """
    genes = [f"G{i:03d}" for i in range(6)]
    denso = np.array(
        [
            [10.0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],  # vazia: o reduceat repetiria 10
            [0, 0, 3.0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
        ],
        dtype=np.float32,
    )
    caminho = tmp_path / "vazia.h5ad"
    escrever(caminho, genes, [(sp.csr_matrix(denso), "D", "ABCC1")])
    maior, acima = conferir_teto_por_celula(caminho, linhas_por_bloco=2)
    assert maior == pytest.approx(10.0)
    assert acima == 0
    # Com teto baixo, so a linha de 10 conta -- nao as vazias herdando o 10.
    _, acima_baixo = conferir_teto_por_celula(caminho, linhas_por_bloco=2, teto=5.0)
    assert acima_baixo == 1
