import json, sys, csv
from pathlib import Path
root = Path(sys.argv[1]); pert_csv = Path(sys.argv[2])
with pert_csv.open(encoding="utf-8") as fh:
    targets = sorted({r["target_gene"] for r in csv.DictReader(fh)})
tr = json.loads((root / "PIE_splits/replogle_xdataset/train.json").read_text(encoding="utf-8"))
by_ds = {}
for key, lst in tr.items():
    by_ds.setdefault(key.split(".", 1)[0], set()).update(lst)
allp = set().union(*by_ds.values())
hit = [t for t in targets if t in allp]
print("VCC targets:", len(targets), "with donor evidence in train.json:", len(hit), "without:", len(targets) - len(hit))
for ds, s in by_ds.items():
    print(ds, sum(1 for t in targets if t in s))
src_keys = {}
for s in ["esm2", "ncbi_text", "string_space", "depmap_gene_effect", "perturbation_text"]:
    p = root / "PIE_sources" / "_adv1_cb1aaa4" / f"{s}_meta.json"
    m = json.loads(p.read_text(encoding="utf-8"))
    ks = set(m["keys"])
    print(s, "targets covered (exact key):", sum(1 for t in targets if t in ks))
gt = set(json.loads((root / "PIE_sources/_cb1aaa4/gene_text/meta.json").read_text(encoding="utf-8"))["keys"])
print("union genes of 4 evidence dirs:")
u = set()
for name in ["PIE_tahoe100m", "PIE_jiang", "PIE_arc_vcc_25", "PIE_x_atlas_orion"]:
    u.update(json.loads((root / name / "meta.json").read_text(encoding="utf-8"))["genes"])
n_p, n_g = len(allp), len(u)
gib = n_p * n_g * 8 / 2**30
print("n_p", n_p, "n_g", n_g, "per (P,G) 8-byte array GiB", round(gib, 2), "build peak est (9 global + 2x4 per-ds int64) GiB", round(gib * (9 + 8), 1))
print("on-disk est GiB: 4 f32 + 5 u16 + 8 u8 =", round(n_p * n_g * (4*4 + 5*2 + 8*1) / 2**30, 1))
