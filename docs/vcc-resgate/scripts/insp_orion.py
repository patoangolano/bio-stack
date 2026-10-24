import json

import numpy as np
import pandas as pd

R = "C:/Users/Manec/bio-stack"
base = f"{R}/data/raw/pie/arrays/PIE_x_atlas_orion/preprocessed"
with open(f"{base}/meta.json", encoding="utf-8") as fh:
    m = json.load(fh)
print("dataset:", m["dataset"], "num_rows:", m["num_rows"], "num_genes:", m["num_genes"])
print("contextos:", m["context_to_id"])
eixo = pd.read_csv(f"{R}/data/raw/vcc/controls-final/gene_names.csv")["gene_name"]
eixo = eixo.astype(str).to_numpy()
alvos = set(pd.read_csv(f"{R}/data/raw/vcc/controls-final/pert_counts.csv")
            ["target_gene"].astype(str))
gp = np.asarray(m["genes"], dtype=str)
pos = {g: i for i, g in enumerate(eixo)}
col = np.array([i for i, g in enumerate(gp) if g in pos])
print("genes orion:", gp.size, "cruzam com o nosso eixo:", col.size)

pid = np.load(f"{base}/pert_ids.npy")
cid = np.load(f"{base}/ctx_ids.npy")
inv = {v: k for k, v in m["pert_to_id"].items()}
nomes = np.array([inv[int(x)] for x in pid])
for ctx_nome, ctx in m["context_to_id"].items():
    sel = cid == ctx
    n_vcc = sum(1 for n in nomes[sel] if n in alvos)
    print(f"  {ctx_nome}: {int(sel.sum())} linhas, {n_vcc} alvos do painel VCC")
