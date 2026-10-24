import anndata as ad, numpy as np, json
from pathlib import Path
VCC = Path("data/raw/vcc/controls-final")
print(json.dumps(json.loads((VCC/"manifest.json").read_text()), indent=1)[:900])
for c in ("D","E","F"):
    a = ad.read_h5ad(VCC / f"context_{c}.h5ad", backed="r")
    print(f"\ncontext_{c}: n_obs={a.n_obs} obs cols={list(a.obs.columns)}")
    for col in a.obs.columns:
        v = a.obs[col]
        u = v.unique()
        print(f"   {col}: {len(u)} valores unicos -> {list(u[:6])}")
    a.file.close()
