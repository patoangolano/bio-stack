import json, sys, csv
from pathlib import Path
root = Path(sys.argv[1])
vcc_csv = Path(sys.argv[2])
gt = set(json.loads((root / "PIE_sources/_cb1aaa4/gene_text/meta.json").read_text(encoding="utf-8"))["keys"])
with vcc_csv.open(encoding="utf-8") as fh:
    vcc = [r["gene_name"] for r in csv.DictReader(fh)]
union = {}
for name in ["PIE_replogle_nadig_essential", "PIE_tahoe100m", "PIE_jiang", "PIE_arc_vcc_25", "PIE_x_atlas_orion"]:
    m = json.loads((root / name / "meta.json").read_text(encoding="utf-8"))
    g = m["genes"]
    miss = [x for x in g if x not in gt]
    print(name, "n_genes", len(g), "missing_from_gene_text", len(miss), miss[:5])
    for x in g:
        union.setdefault(x, len(union))
print("union genes (5 dirs):", len(union))
print("VCC genes in training union:", sum(1 for x in vcc if x in union), "of", len(vcc))
tr = json.loads((root / "PIE_splits/replogle_xdataset/train.json").read_text(encoding="utf-8"))
perts = set()
nrows = 0
per_ds = {}
for key, lst in tr.items():
    ds = key.split(".", 1)[0]
    per_ds.setdefault(ds, set()).update(lst)
    perts.update(lst)
    nrows += len(lst)
print("train.json keys:", len(tr), "rows:", nrows, "unique perts:", len(perts))
print({k: len(v) for k, v in per_ds.items()})
n_p, n_g = len(perts), len(union)
print("one (n_p,n_g) float64 array GiB:", round(n_p * n_g * 8 / 2**30, 2))
