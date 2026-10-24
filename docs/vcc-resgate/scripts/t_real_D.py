import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path
VCC = Path("data/raw/vcc/controls-final")
pg = pd.read_csv("runs/vcc/ctrl_bulk/por_gene.tsv", sep="\t")
genes = pd.read_csv(VCC / "gene_names.csv")["gene_name"].astype(str).tolist()
print("gene col == csv order:", pg["gene"].astype(str).tolist() == genes)

a = ad.read_h5ad(VCC / "context_D.h5ad")
x = sp.csr_matrix(a.X); x.eliminate_zeros()
t = np.asarray(x.sum(axis=1)).ravel().astype(np.float64)
P = np.asarray(x.sum(axis=0)).ravel().astype(np.float64)
T = P.sum()
n = x.shape[0]

b = np.log1p(5e4 * P / T)
cpm_pb = 1e6 * P / T
# media do CPM por celula, feita do jeito mais burro possivel (laco por blocos)
acc = np.zeros(x.shape[1])
for i0 in range(0, n, 2000):
    blo = x[i0:i0+2000].toarray().astype(np.float64)
    acc += (blo * (1e6 / t[i0:i0+2000])[:, None]).sum(axis=0)
cpm_cel = acc / n

print(f"\nb_ctrl      max|dif| vs tsv: {np.abs(b - pg['b_ctrl_D'].to_numpy()).max():.3e}")
print(f"cpm_pb      max|dif| vs tsv: {np.abs(cpm_pb - pg['cpm_pb_D'].to_numpy()).max():.3e}")
print(f"cpm_cel     max|dif| vs tsv: {np.abs(cpm_cel - pg['cpm_cel_D'].to_numpy()).max():.3e}")
print(f"\nchecagem de unidade  b_ctrl == log1p(cpm_pb/20): "
      f"{np.abs(b - np.log1p(cpm_pb/20)).max():.3e}   (5e4/1e6 = 1/20)")

u_gt  = int((cpm_cel >  5.0).sum()); u_ge = int((cpm_cel >= 5.0).sum())
print(f"\nuniverso D: >5.0 -> {u_gt}   >=5.0 -> {u_ge}   oficial 9633")
print(f"genes exatamente em 5.0: {int((cpm_cel == 5.0).sum())}")
d = np.abs(cpm_cel - 5.0)
print(f"genes a menos de 0,01 CPM do corte: {int((d < 0.01).sum())}  "
      f"a menos de 0,1: {int((d < 0.1).sum())}")
print(f"universo no tsv (no_universo_D): {int(pg['no_universo_D'].sum())}")
print(f"\nse_D  tsv: mediana {np.median(pg['se_D']):.4e}  max {pg['se_D'].max():.4e}  "
      f"min {pg['se_D'].min():.4e}  zeros {(pg['se_D']==0).sum()}")
print(f"genes com pseudobulk zero: {int((P==0).sum())}")
