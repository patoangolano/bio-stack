import numpy as np, pandas as pd
VCC="data/raw/vcc/controls-final"
genes=pd.read_csv(f"{VCC}/gene_names.csv")["gene_name"].astype(str).tolist()
with np.load("runs/vcc/ctrl_bulk/ctrl_bulk.npz",allow_pickle=True) as z:
    print("chaves:",z.files)
    for c in ("D","E","F"):
        b=np.asarray(z[f"b_ctrl_{c}"],dtype=np.float64).ravel()
        p=np.expm1(b)
        print(f"{c}: n={b.size} soma_expm1={p.sum():.6f} desvio_de_5e4={p.sum()-5e4:+.3e} zeros={(p<=0).sum()}")
    if "genes" in z.files:
        g=[str(x) for x in np.asarray(z["genes"]).ravel()]
        print("ctrl_bulk genes identico ao painel na ordem:", g==genes)
