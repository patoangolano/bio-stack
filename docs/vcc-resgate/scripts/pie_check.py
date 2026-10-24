import json, sys, csv
from pathlib import Path
import numpy as np
root = Path(sys.argv[1])
vcc_genes_csv = Path(sys.argv[2])
gt = json.loads((root / "PIE_sources/_cb1aaa4/gene_text/meta.json").read_text(encoding="utf-8"))
print("gene_text top keys:", list(gt.keys()))
keys = gt["keys"]
print("gene_text n_keys:", len(keys), "dim:", gt.get("dim"), "layout:", gt.get("layout"), "index:", gt.get("index"))
with vcc_genes_csv.open(encoding="utf-8") as fh:
    vcc = [r["gene_name"] for r in csv.DictReader(fh)]
ks = set(keys)
missing = [g for g in vcc if g not in ks]
print("VCC genes:", len(vcc), "unique:", len(set(vcc)), "missing from gene_text:", len(missing))
print("first missing:", missing[:25])
# dataset metas
for name in ["PIE_replogle_nadig_essential", "PIE_tahoe100m", "PIE_jiang", "PIE_arc_vcc_25", "PIE_x_atlas_orion"]:
    m = json.loads((root / name / "meta.json").read_text(encoding="utf-8"))
    info = {k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in m.items()}
    print(name, info)
