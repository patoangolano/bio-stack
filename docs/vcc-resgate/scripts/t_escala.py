import sys
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, "scripts/vcc")
import ctrl_bulk as cb
ESCALA = cb.ESCALA_PSEUDOBULK

rng = np.random.default_rng(11)
n, g = 18400, 300
# tamanho de biblioteca realista (mediana ~5.000 UMI/celula)
t_alvo = rng.lognormal(np.log(5000), 0.5, n)
# expressao por gene: alguns muito altos (mitocondriais), muitos baixos
peso = rng.gamma(0.3, 1.0, g) + 1e-4
peso[:5] = [0.05, 0.03, 0.02, 0.015, 0.01]   # genes dominantes
peso /= peso.sum()
lam = t_alvo[:, None] * peso[None, :]
xd = rng.poisson(lam).astype(np.float64)
xd[:, -1] = 0.0

P = xd.sum(axis=0); T = P.sum(); tc = xd.sum(axis=1)
b = np.log1p(ESCALA * P / T)

# referencia estavel: dois passes
loo = np.log1p(ESCALA * (P[None, :] - xd) / (T - tc)[:, None])
m = loo.mean(axis=0)
se_ref = np.sqrt((n - 1) / n * np.square(loo - m).sum(axis=0))

x = sp.csr_matrix(xd); x.eliminate_zeros(); c = sp.csc_matrix(x)
se_cod = cb.jackknife_delete1(np.asarray(c.indptr), np.asarray(c.indices),
                              np.asarray(c.data, np.float64), P, tc, float(T))

ok = se_ref > 0
rel = np.abs(se_cod - se_ref)[ok] / se_ref[ok]
ordem = np.argsort(-rel)
print(f"b_ctrl: min {b.min():.3f} max {b.max():.3f}")
print(f"SE ref: mediana {np.median(se_ref):.3e}  min {se_ref[ok].min():.3e}  max {se_ref.max():.3e}")
print(f"erro relativo: mediana {np.median(rel):.3e}  p99 {np.percentile(rel,99):.3e}  max {rel.max():.3e}")
print("\npiores 8 genes (b_ctrl, se_ref, se_codigo, rel):")
idx = np.flatnonzero(ok)[ordem[:8]]
for i in idx:
    print(f"  g{i:4d} b={b[i]:7.4f}  ref={se_ref[i]:.6e}  cod={se_cod[i]:.6e}  rel={np.abs(se_cod[i]-se_ref[i])/se_ref[i]:.3e}")
neg = int((se_cod == 0).sum() - (se_ref == 0).sum())
print(f"\ngenes zerados por np.maximum(desvio,0) a mais que a referencia: {neg}")
