"""Escritor incremental do `.h5ad` da submissao: 900 grupos, nada materializado.

A submissao final tem 360.000 celulas, 18.533 genes e cerca de 2,0e9 nao nulos,
que em `float32` + `int32` sao 15 GB. Nenhum `sp.vstack` sobre mais de um grupo
cabe, logo `X/data` e `X/indices` sao estendidos por `resize` uma vez por grupo
e o `indptr` e acumulado em memoria -- 360.001 entradas em int64, 2,88 MB.

Duas armadilhas, as duas com diagnostico enganoso:

**O `indptr`.** O de cada grupo comeca em 0 e termina no nnz do grupo. Anexar
cru da 900 blocos de 401 entradas -- 360.900 em vez de 360.001 -- que reiniciam
em zero a cada grupo. O `validar_submissao.py` reprova os itens
`estrutura_indptr_tamanho`, `estrutura_indptr_monotonico` e `estrutura_fim_dentro`
de uma vez, e isso se le como "matriz corrompida" em vez de um off-by-one aqui.
O certo e `indptr_local[1:] + deslocamento`.

**`obs` e `var`.** Gravados como dataset h5py a mao, o `_ler_tabela` do validador
devolve DataFrame VAZIO e caem os itens de rotulo juntos -- `colunas_obs`,
`dimensao_de_genes`, `eixo_de_genes`, `contextos`, `alvos_por_contexto`,
`celulas_por_alvo`. Vao por `ad.io.write_elem`, que grava o `encoding-type`
`dataframe` e as categorias do jeito que o leitor espera.

`indptr` em int64 nao e preferencia: 2,0e9 esta a 93,3% de 2^31-1. O
`validar_submissao.py` le com `.astype(np.int64)` e passaria, mas qualquer
consumidor que monte o CSR com `indptr` int32 estoura em silencio.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import cast

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ELEMENTOS_POR_CHUNK: int = 262_144
LINHAS_POR_BLOCO: int = 4096
TETO_DE_CONTAGEM: float = 1_000_000.0


class EscritorCSR:
    """Acumula grupos de 400 celulas num `.h5ad` CSR, um `resize` por grupo."""

    def __init__(
        self,
        caminho: Path,
        genes: Sequence[str],
        col_alvo: str = "target_gene",
        col_contexto: str = "context",
        elementos_por_chunk: int = ELEMENTOS_POR_CHUNK,
    ) -> None:
        self.caminho = Path(caminho)
        self.genes = [str(g) for g in genes]
        self.col_alvo = col_alvo
        self.col_contexto = col_contexto
        self._indptr: list[np.ndarray] = [np.zeros(1, dtype=np.int64)]
        self._alvos: list[str] = []
        self._contextos: list[str] = []
        self._indice: list[str] = []
        self._n_obs = 0
        self._nnz = 0
        self._fechado = False
        self.caminho.parent.mkdir(parents=True, exist_ok=True)
        self._arquivo = h5py.File(self.caminho, "w")
        grupo = self._arquivo.create_group("X")
        grupo.attrs["encoding-type"] = "csr_matrix"
        grupo.attrs["encoding-version"] = "0.1.0"
        grupo.create_dataset(
            "data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(elementos_por_chunk,)
        )
        grupo.create_dataset(
            "indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(elementos_por_chunk,)
        )

    def acrescentar(self, matriz: sp.csr_matrix, contexto: str, alvo: str) -> None:
        if self._fechado:
            raise ValueError("escritor ja fechado")
        local = sp.csr_matrix(matriz)
        forma = cast("tuple[int, int]", local.shape)
        if forma[1] != len(self.genes):
            raise ValueError(
                f"{contexto}/{alvo}: matriz com {forma[1]} colunas, "
                f"eixo de genes tem {len(self.genes)}"
            )
        dados = np.asarray(matriz.data, dtype=np.float64)
        if dados.size:
            if not np.all(np.isfinite(dados)):
                raise ValueError(f"{contexto}/{alvo}: valor nao finito")
            if dados.min() < 0:
                raise ValueError(f"{contexto}/{alvo}: contagem negativa")
            if not np.array_equal(dados, np.floor(dados)):
                raise ValueError(f"{contexto}/{alvo}: contagem nao inteira")
        local.sort_indices()

        x = self._arquivo["X"]
        assert isinstance(x, h5py.Group)
        d, i = x["data"], x["indices"]
        assert isinstance(d, h5py.Dataset) and isinstance(i, h5py.Dataset)
        d.resize((self._nnz + local.nnz,))
        i.resize((self._nnz + local.nnz,))
        d[self._nnz :] = np.asarray(local.data, dtype=np.float32)
        i[self._nnz :] = np.asarray(local.indices, dtype=np.int32)

        # A linha que o teste existe para fixar: o `indptr` local SEM o zero
        # inicial, deslocado pelo nnz ja escrito.
        self._indptr.append(np.asarray(local.indptr[1:], dtype=np.int64) + self._nnz)
        self._nnz += int(local.nnz)
        self._alvos.extend([alvo] * forma[0])
        self._contextos.extend([contexto] * forma[0])
        self._indice.extend(f"{contexto}_{alvo}_{k:03d}" for k in range(forma[0]))
        self._n_obs += forma[0]

    def fechar(self) -> int:
        if self._fechado:
            return self._nnz
        indptr = np.concatenate(self._indptr)
        if indptr.size != self._n_obs + 1:
            raise ValueError(f"indptr com {indptr.size} entradas, esperado {self._n_obs + 1}")
        obs = pd.DataFrame(
            {
                self.col_alvo: pd.Categorical(self._alvos),
                self.col_contexto: pd.Categorical(self._contextos),
            },
            index=pd.Index(self._indice, dtype=object, name=None),
        )
        var = pd.DataFrame(index=pd.Index(self.genes, dtype=object, name=None))
        x = self._arquivo["X"]
        assert isinstance(x, h5py.Group)
        x.create_dataset("indptr", data=indptr, dtype=np.int64)
        x.attrs["shape"] = np.array([self._n_obs, len(self.genes)], dtype=np.int64)
        self._arquivo.attrs["encoding-type"] = "anndata"
        self._arquivo.attrs["encoding-version"] = "0.1.0"
        ad.io.write_elem(self._arquivo, "obs", obs)
        ad.io.write_elem(self._arquivo, "var", var)
        for vazio in ("obsm", "varm", "obsp", "varp", "layers", "uns"):
            ad.io.write_elem(self._arquivo, vazio, {})
        self._arquivo.close()
        self._fechado = True
        return self._nnz

    def __enter__(self) -> EscritorCSR:
        return self

    def __exit__(self, *_: object) -> None:
        if not self._fechado:
            self.fechar()


def conferir_teto_por_celula(
    caminho: Path, linhas_por_bloco: int = LINHAS_POR_BLOCO, teto: float = TETO_DE_CONTAGEM
) -> tuple[float, int]:
    """`(maior soma por celula, quantas celulas acima do teto)`, varrido no ARQUIVO.

    A guarda e medida sobre o que foi escrito, nunca reconstruida de um
    `max(kernel, controle)` em memoria -- no codigo de referencia do organizador
    essa reconstrucao divergiu do real em 2 de 4 sementes.
    """
    with h5py.File(Path(caminho), "r") as f:
        x = f["X"]
        assert isinstance(x, h5py.Group)
        dados = x["data"]
        assert isinstance(dados, h5py.Dataset)
        indptr = np.asarray(x["indptr"]).astype(np.int64)
        n_linhas = indptr.size - 1
        maior = 0.0
        acima = 0
        for r0 in range(0, n_linhas, linhas_por_bloco):
            r1 = min(r0 + linhas_por_bloco, n_linhas)
            a, b = int(indptr[r0]), int(indptr[r1])
            if b == a:
                continue
            valores = np.asarray(dados[a:b], dtype=np.float64)
            # `np.add.reduceat` nao serve aqui, e o teste de linha vazia mostrou
            # por que: ele REPETE o valor anterior numa linha vazia e levanta
            # `IndexError: index out-of-bounds` quando as linhas vazias estao no
            # fim do bloco (o indice de inicio delas e igual ao tamanho do
            # pedaco). `repeat` + `bincount` nao tem nenhum dos dois casos --
            # contagem zero simplesmente nao gera indice.
            por_linha = np.diff(indptr[r0 : r1 + 1])
            linhas = np.repeat(np.arange(r1 - r0), por_linha)
            somas = np.bincount(linhas, weights=valores, minlength=r1 - r0)
            maior = max(maior, float(somas.max()))
            acima += int((somas > teto).sum())
    return maior, acima
