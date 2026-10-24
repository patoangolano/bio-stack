import numpy as np, anndata as ad, scipy.sparse as sp, pandas as pd
from pathlib import Path
VCC = Path("data/raw/vcc/controls-final")
a = ad.read_h5ad(VCC / "context_D.h5ad")
print("tipo X:", type(a.X).__name__, "dtype:", a.X.dtype)
x = sp.csr_matrix(a.X); x.eliminate_zeros()
print("dtype apos csr_matrix:", x.dtype, " nnz:", x.nnz)

# caminho do codigo: soma na dtype da matriz, cast depois
pb_codigo = np.asarray(x.sum(axis=0), dtype=np.float64).ravel()
print("dtype de x.sum(axis=0):", np.asarray(x.sum(axis=0)).dtype)
# referencia: soma de verdade em float64
x64 = x.astype(np.float64)
pb_ref = np.asarray(x64.sum(axis=0)).ravel()

d = np.abs(pb_codigo - pb_ref)
print(f"\npseudobulk max gene: {pb_ref.max():.6e}   (limite inteiro exato float32 = 1.6777e7)")
print(f"genes acima de 2^24: {int((pb_ref > 2**24).sum())}")
print(f"max|dif| pseudobulk codigo vs float64: {d.max():.6e}  genes diferentes: {int((d>0).sum())}")
if d.max() > 0:
    i = int(np.argmax(d))
    print(f"  pior gene idx={i}: codigo={pb_codigo[i]:.1f} ref={pb_ref[i]:.1f} dif={d[i]:.1f}")
    b_c = np.log1p(5e4*pb_codigo/pb_codigo.sum()); b_r = np.log1p(5e4*pb_ref/pb_ref.sum())
    print(f"  impacto em b_ctrl: max|dif|={np.abs(b_c-b_r).max():.3e}")
    cc = 1e6*pb_codigo/pb_codigo.sum(); cr = 1e6*pb_ref/pb_ref.sum()
    print(f"  impacto no universo (cpm_pb>5): {int((cc>5).sum())} vs {int((cr>5).sum())}")

# a outra soma: umi_por_celula
t_c = np.asarray(x.sum(axis=1), dtype=np.float64).ravel()
t_r = np.asarray(x64.sum(axis=1)).ravel()
print(f"\numi_por_celula max: {t_r.max():.1f}  max|dif|: {np.abs(t_c-t_r).max():.3e}")

# e a media do CPM por celula, que usa diags @ x em float32?
esc = sp.csr_matrix(sp.diags(1e6/t_c) @ x)
print(f"dtype de diags(1e6/umi) @ x: {esc.dtype}")
cm_c = np.asarray(esc.mean(axis=0)).ravel()
esc64 = sp.csr_matrix(sp.diags(1e6/t_r) @ x64)
cm_r = np.asarray(esc64.mean(axis=0)).ravel()
print(f"cpm_medio_celula max|dif|: {np.abs(cm_c-cm_r).max():.3e}  "
      f"universo: {int((cm_c>5).sum())} vs {int((cm_r>5).sum())}")
