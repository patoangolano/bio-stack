import json, sys
import numpy as np
base = "C:/Users/Manec/bio-stack/data/raw/pie/arrays/PIE_arc_vcc_25/preprocessed"
m = json.load(open(base + "/meta.json"))
print("meta keys:", list(m.keys()))
print("n genes in meta:", len(m.get("genes", [])))
for k, v in m.items():
    if k != "genes":
        print(" ", k, "=", (str(v)[:200] if not isinstance(v, list) else f"list[{len(v)}] {v[:6]}"))
for f in ("fold_changes", "delta_p", "fdr", "tested", "ctrl_means", "pert_ids", "ctx_ids"):
    a = np.load(f"{base}/{f}.npy", allow_pickle=True)
    print(f"{f}: shape={a.shape} dtype={a.dtype}", a.ravel()[:4] if a.dtype.kind in "fiub" else a.ravel()[:4])
