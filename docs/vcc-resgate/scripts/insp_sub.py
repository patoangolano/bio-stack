import anndata as ad
import numpy as np

p = "/mnt/c/Users/Manec/bio-stack/runs/vcc/submissao/submissao_v1.h5ad"
a = ad.read_h5ad(p, backed="r")
print("forma:", a.shape)
print("obs colunas:", list(a.obs.columns))
print(a.obs.head(3))
for c in a.obs.columns:
    v = a.obs[c].astype(str).to_numpy()
    u, n = np.unique(v, return_counts=True)
    print(f"  {c}: {u.size} valores; contagens distintas {np.unique(n)[:6]}")
    print(f"     primeiros: {list(u[:6])}")
print("var_names[:4]:", list(a.var_names[:4]), "n_var:", a.n_vars)
x = a.X
print("tipo de X:", type(x).__name__)
bloco = a[:3].to_memory().X
print("dtype:", bloco.dtype, "amostra:", np.asarray(bloco[0].todense()).ravel()[:6])
