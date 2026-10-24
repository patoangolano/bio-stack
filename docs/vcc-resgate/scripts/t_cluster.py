import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path
VCC = Path("data/raw/vcc/controls-final")
a = ad.read_h5ad(VCC / "context_D.h5ad")
x = sp.csr_matrix(a.X); x.eliminate_zeros()
g = a.obs["ntc_id"].astype(str).to_numpy()
t = np.asarray(x.sum(axis=1)).ravel().astype(np.float64)
P = np.asarray(x.sum(axis=0)).ravel().astype(np.float64); T = P.sum()
grupos = np.unique(g); G = grupos.size
print(f"grupos (ntc_id): {G} | celulas/grupo: min {np.bincount(pd.factorize(g)[0]).min()} "
      f"max {np.bincount(pd.factorize(g)[0]).max()}")

loo = np.empty((G, x.shape[1]))
for k, grp in enumerate(grupos):
    sel = g == grp
    Pk = P - np.asarray(x[sel].sum(axis=0)).ravel()
    Tk = T - t[sel].sum()
    loo[k] = np.log1p(5e4 * Pk / Tk)
m = loo.mean(axis=0)
se_grupo = np.sqrt((G - 1) / G * np.square(loo - m).sum(axis=0))

se_cel = pd.read_csv("runs/vcc/ctrl_bulk/por_gene.tsv", sep="\t")["se_D"].to_numpy()
ok = (se_cel > 0) & (se_grupo > 0)
r = se_grupo[ok] / se_cel[ok]
print(f"\nSE delete-1-celula  : mediana {np.median(se_cel[ok]):.4e}")
print(f"SE delete-1-guia(46): mediana {np.median(se_grupo[ok]):.4e}")
print(f"razao guia/celula   : mediana {np.median(r):.3f}  p10 {np.percentile(r,10):.3f}  "
      f"p90 {np.percentile(r,90):.3f}")

print("\n--- casos degenerados de jackknife_delete1 ---")
import sys; sys.path.insert(0,"scripts/vcc"); import ctrl_bulk as cb
for n in (1, 2):
    xd = np.array([[3.0,0.0],[5.0,2.0]])[:n]
    c = sp.csc_matrix(sp.csr_matrix(xd))
    Pp = xd.sum(axis=0); tt = xd.sum(axis=1)
    out = cb.jackknife_delete1(np.asarray(c.indptr), np.asarray(c.indices),
          np.asarray(c.data,np.float64), Pp, tt, float(Pp.sum()))
    print(f"  n_celulas={n} -> se={out}")
