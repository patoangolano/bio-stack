import json
from pathlib import Path
import numpy as np, pandas as pd
VCC=Path("data/raw/vcc/controls-final"); T=Path("runs/vcc/transfer"); CTRL=Path("runs/vcc/ctrl_bulk")
genes_painel=pd.read_csv(VCC/"gene_names.csv")["gene_name"].astype(str).tolist()
alvos=pd.read_csv(VCC/"pert_counts.csv")["target_gene"].astype(str).tolist()
pg=pd.read_csv(CTRL/"por_gene.tsv",sep="\t")
print("colunas por_gene:",list(pg.columns))
print("linhas por_gene:",len(pg))
print("ordem bate com gene_names.csv:", list(pg[pg.columns[0]].astype(str))==genes_painel if pg.columns[0] in ("gene","gene_name") else "verificar")
for c in ("D","E","F"):
    b=np.asarray(pg[f"b_ctrl_{c}"],dtype=np.float64)
    print(f"{c}: b_ctrl zeros={int((b==0).sum())} min={b.min():.6f} max={b.max():.4f} soma expm1={np.expm1(b).sum():.1f}")
