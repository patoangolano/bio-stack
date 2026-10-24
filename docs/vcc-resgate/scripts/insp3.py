import json
import numpy as np
import pandas as pd
R = "C:/Users/Manec/bio-stack"
base = f"{R}/data/raw/pie/arrays/PIE_arc_vcc_25/preprocessed"
m = json.load(open(base + "/meta.json"))
eixo = pd.read_csv(f"{R}/data/raw/vcc/controls-final/gene_names.csv")["gene_name"].astype(str).to_numpy()
alvos = set(pd.read_csv(f"{R}/data/raw/vcc/controls-final/pert_counts.csv")["target_gene"].astype(str))
fc = np.load(base+"/fold_changes.npy"); tested = np.load(base+"/tested.npy")
pid = np.load(base+"/pert_ids.npy"); cid = np.load(base+"/ctx_ids.npy")
genes_pie = np.asarray(m["genes"], dtype=str)
pos = {g:i for i,g in enumerate(eixo)}
col = np.array([i for i,g in enumerate(genes_pie) if g in pos])
inv = {v:k for k,v in m["pert_to_id"].items()}
nomes = [inv[int(x)] for x in pid]
print("perts arc_vcc_25 total:", len(set(nomes)))
print("por contexto:", {c: int((cid==i).sum()) for c,i in m["context_to_id"].items()})
inter = sorted({n for n in nomes if n in alvos})
print("intersecao com os 300 alvos VCC:", len(inter))
cont = []
for lin in range(fc.shape[0]):
    n = inv[int(pid[lin])]
    if n not in alvos: continue
    v = fc[lin, col]; t = tested[lin, col]
    mv = t & np.isfinite(v) & (v>0) & (np.abs(np.log2(np.where(v>0,v,1.0)))>0.263)
    mv2 = t & np.isfinite(v) & (v>0) & (np.abs(np.log2(np.where(v>0,v,1.0)))>0.1)
    cont.append((int(mv.sum()), int(mv2.sum()), n, int(cid[lin])))
cont.sort(reverse=True)
print(f"linhas na intersecao: {len(cont)}")
q = np.array([c[0] for c in cont]); q2 = np.array([c[1] for c in cont])
print("genes moventes (|log2fc|>0.263): min/med/max =", q.min(), int(np.median(q)), q.max())
print("genes moventes (|log2fc|>0.1)  : min/med/max =", q2.min(), int(np.median(q2)), q2.max())
for th in (10, 50, 100, 200, 500):
    print(f"  >= {th:4d} moventes(0.263): {(q>=th).sum():4d}   moventes(0.1): {(q2>=th).sum():4d}")
print("top 10:", cont[:10])
