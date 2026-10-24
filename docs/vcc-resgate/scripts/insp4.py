import json
import numpy as np
import pandas as pd
R = "C:/Users/Manec/bio-stack"
base = f"{R}/data/raw/pie/arrays/PIE_arc_vcc_25/preprocessed"
m = json.load(open(base + "/meta.json"))
eixo = pd.read_csv(f"{R}/data/raw/vcc/controls-final/gene_names.csv")["gene_name"].astype(str).to_numpy()
alvos = set(pd.read_csv(f"{R}/data/raw/vcc/controls-final/pert_counts.csv")["target_gene"].astype(str))
fc=np.load(base+"/fold_changes.npy"); te=np.load(base+"/tested.npy")
pid=np.load(base+"/pert_ids.npy"); cid=np.load(base+"/ctx_ids.npy")
gp=np.asarray(m["genes"],dtype=str); pos={g:i for i,g in enumerate(eixo)}
col=np.array([i for i,g in enumerate(gp) if g in pos]); inv={v:k for k,v in m["pert_to_id"].items()}
for ctx_nome, ctx in m["context_to_id"].items():
    res=[]
    for lin in range(fc.shape[0]):
        if int(cid[lin])!=ctx: continue
        n=inv[int(pid[lin])]
        v=fc[lin,col]; t=te[lin,col]
        mv=t&np.isfinite(v)&(v>0)&(np.abs(np.log2(np.where(v>0,v,1.0)))>0.1)
        res.append((int(mv.sum()), n, n in alvos))
    res.sort(reverse=True)
    q=np.array([r[0] for r in res]); vcc=sum(1 for r in res if r[2])
    print(f"{ctx_nome}: {len(res)} perts, {vcc} entre os 300 VCC")
    for th in (100,200,400,800):
        print(f"   >= {th:4d} moventes: {(q>=th).sum():4d}  (dos VCC: {sum(1 for r in res if r[0]>=th and r[2]):3d})")
