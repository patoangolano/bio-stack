import numpy as np
z = np.load("runs/vcc/ctrl_bulk/ctrl_bulk.npz", allow_pickle=True)
for c in ("D","E","F"):
    pb = z[f"pseudobulk_{c}"]
    print(f"  {c}: pseudobulk max gene {pb.max():.4e}  margem ate 2^24: {2**24/pb.max():.2f}x")
# demonstracao do limite: soma float32 de inteiros acima de 2^24 perde contagem
n = 40000
v = np.full(n, 500.0, dtype=np.float32)   # soma exata = 2.0e7 > 2^24
print(f"\n  soma float32 de {n}x500 = {np.float32(v).sum():.1f} | exato = {n*500}")
import scipy.sparse as sp
m = sp.csr_matrix(v.reshape(-1,1))
print(f"  via scipy csr .sum(axis=0) float32: {float(np.asarray(m.sum(axis=0)).ravel()[0]):.1f}")
print(f"  via float64:                        {float(np.asarray(m.astype(np.float64).sum(axis=0)).ravel()[0]):.1f}")
