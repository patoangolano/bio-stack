import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "scripts/vcc")
METAS = Path("data/raw/pie"); ARRAYS = Path("data/raw/pie/arrays")
VCC = Path("data/raw/vcc/controls-final")
REPOS = {"orion":"PIE_x_atlas_orion","replogle":"PIE_replogle_nadig_essential",
         "arc_vcc_25":"PIE_arc_vcc_25","jiang":"PIE_jiang"}
prof = {"arc_vcc_25":2016.0,"replogle":46.0,"jiang":13.0,"orion":2.0}
alvos = set(pd.read_csv(VCC/"pert_counts.csv")["target_gene"].astype(str))
split = json.loads((METAS/"PIE_splits"/"replogle_xdataset"/"train.json").read_text(encoding="utf-8"))
treino = {}
for k, v in split.items():
    treino.setdefault(str(k).split(".",1)[0], set()).update(str(p) for p in v)
print("chaves do split:", sorted(split)[:8], "...", len(split))
print("datasets no split:", {k: len(v) for k,v in treino.items()})

def perts_por_contexto(repo):
    rd = REPOS[repo]
    meta = json.loads((METAS/rd/"meta.json").read_text(encoding="utf-8"))
    ip = {int(i):str(p) for p,i in meta["pert_to_id"].items()}
    ic = {int(i):str(c) for c,i in meta["context_to_id"].items()}
    cx = np.load(ARRAYS/rd/"preprocessed"/"ctx_ids.npy", allow_pickle=False)
    px = np.load(ARRAYS/rd/"preprocessed"/"pert_ids.npy", allow_pickle=False)
    out = {}
    for c,p in zip(cx.tolist(), px.tolist(), strict=True):
        out.setdefault(ic[int(c)], set()).add(ip[int(p)])
    return out

peso_media = {}
total_ctx = 0
cobertura = {}
for repo in REPOS:
    meta = json.loads((METAS/REPOS[repo]/"meta.json").read_text(encoding="utf-8"))
    ds = str(meta["dataset"])
    cat = perts_por_contexto(repo)
    for ctx in sorted(cat):
        usaveis = sorted(alvos & cat[ctx] & treino.get(ds, set()))
        if not usaveis: continue
        total_ctx += 1
        peso_media[f"{repo}.{ctx}"] = (prof[repo], len(usaveis))
        for a in usaveis: cobertura.setdefault(a, []).append(f"{repo}.{ctx}")
print(f"\ncontextos que contribuem: {total_ctx}")
for k,(w,n) in peso_media.items(): print(f"  {k:24s} peso={w:7.1f} alvos={n}")
import collections
soma = collections.Counter()
for k,(w,n) in peso_media.items(): soma[k.split(".")[0]] += w
tot = sum(soma.values())
print("\npeso TOTAL acumulado na media agrupada, por repo (so onde media!=0):")
for r,s in sorted(soma.items(), key=lambda x:-x[1]):
    print(f"  {r:12s} {s:9.1f}  ({100*s/tot:.1f}% se todos os genes tiverem media em todos)")
print(f"\nalvos distintos cobertos: {len(cobertura)}")
mult = {a:v for a,v in cobertura.items() if len(v)>1}
print(f"alvos com mais de um contexto doador: {len(mult)}")
for a,v in list(mult.items())[:12]: print("   ", a, v)
