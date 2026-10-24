import numpy as np, anndata as ad, scipy.sparse as sp
from pathlib import Path
VCC = Path("data/raw/vcc/controls-final")
for c in ("E","F"):
    a = ad.read_h5ad(VCC / f"context_{c}.h5ad")
    x = sp.csr_matrix(a.X); x.eliminate_zeros()
    pb32 = np.asarray(x.sum(axis=0), dtype=np.float64).ravel()
    pb64 = np.asarray(x.astype(np.float64).sum(axis=0)).ravel()
    d = np.abs(pb32-pb64)
    print(f"{c}: dtype {x.dtype} max gene {pb64.max():.6e} | max|dif| {d.max():.3e} "
          f"| genes afetados {int((d>0).sum())}")
    t32 = np.asarray(x.sum(axis=1), dtype=np.float64).ravel()
    t64 = np.asarray(x.astype(np.float64).sum(axis=1)).ravel()
    print(f"   umi/celula max|dif| {np.abs(t32-t64).max():.3e}")
    del a, x

# teste mais afiado: soma float32 de contagens inteiras com total acima de 2^24
rng = np.random.default_rng(3)
for alvo in (1.2e7, 1.6e7, 2.0e7, 4.0e7, 8.0e7):
    v = rng.poisson(alvo/18400, 18400).astype(np.float32)
    s32 = float(np.asarray(sp.csr_matrix(v.reshape(-1,1)).sum(axis=0)).ravel()[0])
    s64 = float(v.astype(np.float64).sum())
    print(f"total~{alvo:.1e}: float32 {s32:.1f} vs float64 {s64:.1f}  erro {s32-s64:+.1f} "
          f"({abs(s32-s64)/s64:.2e} rel)")
