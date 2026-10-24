import numpy as np, pandas as pd
pg = pd.read_csv("runs/vcc/ctrl_bulk/por_gene.tsv", sep="\t")
z = np.load("runs/vcc/ctrl_bulk/ctrl_bulk.npz", allow_pickle=True)
print("gene order npz == tsv:", [str(g) for g in z["genes"]] == pg["gene"].astype(str).tolist())
for c in ("D","E","F"):
    b_tsv, b_npz = pg[f"b_ctrl_{c}"].to_numpy(), z[f"b_ctrl_{c}"]
    s = np.expm1(b_tsv).sum()
    print(f"{c}: sum(expm1(b_ctrl)) = {s:.6f}  (esperado 5e4, desvio rel {abs(s-5e4)/5e4:.2e}) | "
          f"tsv vs npz max|dif| {np.abs(b_tsv-b_npz).max():.3e} | "
          f"genes com b_ctrl==0: {int((b_tsv==0).sum())}")
