"""Falsificacoes de submissao do VCC 2026 para testar `validar_submissao.py`.

Porte de `docs/vcc-resgate/scripts/t4_forja.py` e `t8_rotulos.py`, com duas
mudancas deliberadas:

* Cada variante tem o PROPRIO gerador (`default_rng(7)`). No original o `rng`
  era global ao modulo, entao o conteudo de uma falsificacao dependia de quais
  outras tinham sido geradas antes no mesmo processo -- e o `pytest-randomly`
  do projeto reordena os testes.
* As colunas de cada celula sao sorteadas por estrato (um gene por faixa de
  `NG // K`), o que as torna estritamente crescentes por construcao e dispensa
  o laco de 120.000 linhas do original. Os defeitos injetados continuam nas
  mesmas posicoes, entao os laudos antigos valem como oraculo.

Os laudos esperados sao os `laudo_*.json` gerados pelo validador antigo sobre
as falsificacoes originais, resgatados em `docs/vcc-resgate/arrays/`.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

RAIZ = Path(__file__).resolve().parents[2]
REF = RAIZ / "data" / "raw" / "vcc" / "controls-final"
LAUDOS = RAIZ / "docs" / "vcc-resgate" / "arrays"
sys.path.insert(0, str(RAIZ / "scripts" / "vcc"))

# As 15 variantes que o validador antigo julgou e cujos laudos sobreviveram.
VARIANTES = [
    "limpo",
    "celula_nula",
    "contagem",
    "desordem",
    "eixo_invertido",
    "folga",
    "fracionario",
    "indice_fora",
    "indptr_nao_mono",
    "nan",
    "negativo",
    "ntc",
    "pert_nula",
    "repetido",
    "zero",
]

NC = 400  # celulas por (contexto, alvo)
K = 20  # nao nulos por celula


def _referencia() -> tuple[list[str], list[str]]:
    if not (REF / "gene_names.csv").exists():
        pytest.skip(f"pacote controls-final ausente em {REF} (vcc datasets download controls-final)")
    genes = pd.read_csv(REF / "gene_names.csv")["gene_name"].astype(str).tolist()
    alvos = pd.read_csv(REF / "pert_counts.csv")["target_gene"].astype(str).tolist()
    return genes, alvos


def laudo_antigo(variante: str) -> dict:
    return json.loads((LAUDOS / f"laudo_{variante}.json").read_text(encoding="utf-8"))


def _x_base(rng: np.random.Generator, nobs: int, ng: int):
    largura = ng // K
    cols = np.arange(K) * largura + rng.integers(0, largura, size=(nobs, K))
    dat = rng.integers(1, 50, size=(nobs, K)).astype(np.float32).ravel()
    ind = cols.astype(np.int32).ravel()
    indptr = np.arange(0, nobs * K + 1, K, dtype=np.int64)
    return dat, ind, indptr


def _defeito_em_x(variante: str, dat, ind, indptr, ng: int):
    if variante == "repetido":  # indice duplicado na linha 5
        ind[5 * K + 1] = ind[5 * K]
    elif variante == "desordem":  # troca duas colunas na linha 9
        ind[9 * K], ind[9 * K + 1] = ind[9 * K + 1], ind[9 * K]
    elif variante == "fracionario":
        dat[100] = 3.5
    elif variante == "negativo":
        dat[200] = -1.0
    elif variante == "nan":
        dat[300] = np.nan
    elif variante == "zero":
        dat[400] = 0.0
    elif variante == "celula_nula":  # linha 7 inteira em zero
        dat[7 * K : 8 * K] = 0.0
    elif variante == "pert_nula":  # o primeiro alvo inteiro em zero
        dat[0 : NC * K] = 0.0
    elif variante == "folga":
        ind = np.concatenate([ind, np.zeros(3, np.int32)])
        dat = np.concatenate([dat, np.zeros(3, np.float32)])
    elif variante == "indice_fora":
        ind[500] = ng + 5
    elif variante == "indptr_nao_mono":
        indptr[10] = indptr[11] + 5
    return dat, ind, indptr


def _obs(variante: str, alvos: list[str]) -> pd.DataFrame:
    ca: list[str] = []
    for a in alvos:
        ca.extend([a] * NC)
    if variante == "ntc":  # 5 linhas de controle, tiradas do primeiro alvo
        ca[0:5] = ["non-targeting"] * 5
    if variante == "contagem":  # alvos[0] com 395, alvos[1] com 405
        ca[0:5] = [alvos[1]] * 5
    # As variantes de rotulo usam categorias em ordem alfabetica e as demais na
    # ordem oficial, como os originais: o validador nao pode depender disso.
    categorias = sorted(set(ca)) if variante in ("ntc", "contagem") else alvos
    colunas: dict[str, pd.Categorical] = {
        "target_gene": pd.Categorical(ca, categories=categorias),
        "context": pd.Categorical(["D"] * len(ca), categories=["D"]),
    }
    if variante == "sem_contexto":
        del colunas["context"]
    return pd.DataFrame(colunas, index=pd.Index([f"c{i}" for i in range(len(ca))], dtype=object))


def escrever_forja(variante: str, caminho: Path) -> Path:
    genes, alvos = _referencia()
    ng = len(genes)
    obs = _obs(variante, alvos)
    rng = np.random.default_rng(7)
    dat, ind, indptr = _x_base(rng, len(obs), ng)
    dat, ind, indptr = _defeito_em_x(variante, dat, ind, indptr, ng)
    eixo = list(reversed(genes)) if variante == "eixo_invertido" else genes
    var = pd.DataFrame(index=pd.Index(eixo, dtype=object))
    ad.AnnData(obs=obs, var=var).write_h5ad(caminho)
    with h5py.File(caminho, "r+") as f:
        if "X" in f:
            del f["X"]
        g = f.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.array([len(obs), ng], dtype=np.int64)
        g.create_dataset("data", data=dat.astype(np.float32))
        g.create_dataset("indices", data=ind.astype(np.int32))
        g.create_dataset("indptr", data=indptr.astype(np.int64))
    return caminho


@pytest.fixture(scope="session")
def ref() -> Path:
    _referencia()
    return REF


@pytest.fixture(scope="session")
def forja(tmp_path_factory: pytest.TempPathFactory) -> Callable[[str], Path]:
    """Devolve o caminho da falsificacao pedida, construindo-a uma vez por sessao."""
    pasta = tmp_path_factory.mktemp("forjas")
    feitas: dict[str, Path] = {}

    def obter(variante: str) -> Path:
        if variante not in feitas:
            feitas[variante] = escrever_forja(variante, pasta / f"forja_{variante}.h5ad")
        return feitas[variante]

    return obter
