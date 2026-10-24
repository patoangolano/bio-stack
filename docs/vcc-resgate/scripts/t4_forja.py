"""Forja um h5ad que satisfaz o contrato e injeta defeitos, um por arquivo."""
import sys, numpy as np, pandas as pd, anndata as ad, h5py, scipy.sparse as sp
from pathlib import Path
sys.path.insert(0, "scripts/vcc")
OUT = Path(sys.argv[1])
genes = pd.read_csv("data/raw/vcc/controls-final/gene_names.csv")["gene_name"].astype(str).tolist()
alvos = pd.read_csv("data/raw/vcc/controls-final/pert_counts.csv")["target_gene"].astype(str).tolist()
NG, NA, NC = len(genes), len(alvos), 400
nobs = NA * NC
rng = np.random.default_rng(7)
K = 20  # nao nulos por celula

def obs_df():
    nomes, ca, cc = [], [], []
    for a in alvos:
        nomes.extend(f"D_{a}_{i:03d}" for i in range(NC)); ca.extend([a]*NC); cc.extend(["D"]*NC)
    return pd.DataFrame({"target_gene": pd.Categorical(ca, categories=alvos),
                         "context": pd.Categorical(cc, categories=["D"])},
                        index=pd.Index(nomes, dtype=object))

def construir(variante):
    cols = np.sort(rng.choice(NG, size=(nobs, K), replace=True), axis=1)
    # garante indices estritamente crescentes por linha (sem repetido)
    for r in range(nobs):
        u = np.unique(cols[r])
        if u.size < K: u = np.sort(rng.permutation(NG)[:K])
        cols[r] = u
    data = rng.integers(1, 50, size=(nobs, K)).astype(np.float32)
    indptr = np.arange(0, nobs*K+1, K, dtype=np.int64)
    ind = cols.ravel().astype(np.int32); dat = data.ravel()
    if variante == "repetido":       # indice duplicado na linha 5
        ind[5*K+1] = ind[5*K]
    elif variante == "desordem":     # troca duas colunas na linha 9
        ind[9*K], ind[9*K+1] = ind[9*K+1], ind[9*K]
    elif variante == "fracionario":
        dat[100] = 3.5
    elif variante == "negativo":
        dat[200] = -1.0
    elif variante == "nan":
        dat[300] = np.nan
    elif variante == "zero":
        dat[400] = 0.0
    elif variante == "celula_nula":  # linha 7 inteira em zero -> soma zero
        dat[7*K:8*K] = 0.0
    elif variante == "pert_nula":    # alvo 0 inteiro em zero
        dat[0:NC*K] = 0.0
    elif variante == "folga":
        ind = np.concatenate([ind, np.zeros(3, np.int32)]); dat = np.concatenate([dat, np.zeros(3, np.float32)])
    elif variante == "indice_fora":
        ind[500] = NG + 5
    elif variante == "indptr_nao_mono":
        indptr[10] = indptr[11] + 5
    elif variante == "eixo_invertido":  # var na ordem reversa
        pass
    return dat, ind, indptr

def escrever(variante, caminho):
    dat, ind, indptr = construir(variante)
    gs = list(reversed(genes)) if variante == "eixo_invertido" else genes
    vazio = ad.AnnData(obs=obs_df(), var=pd.DataFrame(index=pd.Index(gs, dtype=object)))
    vazio.write_h5ad(caminho, compression=None)
    with h5py.File(caminho, "r+") as f:
        if "X" in f: del f["X"]
        g = f.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"; g.attrs["encoding-version"] = "0.1.0"
        g.attrs["shape"] = np.array([nobs, NG], dtype=np.int64)
        g.create_dataset("data", data=dat.astype(np.float32))
        g.create_dataset("indices", data=ind.astype(np.int32))
        g.create_dataset("indptr", data=indptr.astype(np.int64))

for v in sys.argv[2:]:
    p = OUT / f"forja_{v}.h5ad"
    escrever(v, p)
    print("escrito", p, p.stat().st_size // 1024, "KiB")
