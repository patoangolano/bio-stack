import anndata as ad, pandas as pd, numpy as np
from pathlib import Path
VCC = Path("data/raw/vcc/controls-final")
genes = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
print(f"gene_names.csv: {len(genes)} genes, primeiros {genes[:4]}")
ref = None
for c in ("D", "E", "F"):
    a = ad.read_h5ad(VCC / f"context_{c}.h5ad", backed="r")
    vn = [str(v) for v in a.var_names]
    igual_ordem = vn == genes
    igual_conj = set(vn) == set(genes)
    print(f"context_{c}: n_vars={a.n_vars} n_obs={a.n_obs} | "
          f"ordem==csv: {igual_ordem} | conjunto==csv: {igual_conj} | primeiros {vn[:4]}")
    if ref is None:
        ref = vn
    else:
        print(f"           ordem igual a context_D: {vn == ref}")
    a.file.close()
