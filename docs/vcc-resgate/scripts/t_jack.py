import sys
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "scripts/vcc")
import ctrl_bulk as cb

ESCALA = cb.ESCALA_PSEUDOBULK
rng = np.random.default_rng(7)

def naive(xd):
    n, g = xd.shape
    P = xd.sum(axis=0)
    T = P.sum()
    loo = np.empty((n, g))
    for i in range(n):
        Pi = P - xd[i]
        Ti = T - xd[i].sum()
        loo[i] = np.log1p(ESCALA * Pi / Ti)
    media = loo.mean(axis=0)
    var = (n - 1) / n * np.square(loo - media).sum(axis=0)
    return np.sqrt(var), loo

def fechada(xd):
    x = sp.csr_matrix(xd)
    x.eliminate_zeros()
    P = np.asarray(x.sum(axis=0), dtype=np.float64).ravel()
    t = np.asarray(x.sum(axis=1), dtype=np.float64).ravel()
    T = float(P.sum())
    c = sp.csc_matrix(x)
    return cb.jackknife_delete1(
        np.asarray(c.indptr), np.asarray(c.indices),
        np.asarray(c.data, dtype=np.float64), P, t, T)

for nome, xd in [
    ("esparso 200x50", (rng.poisson(0.4, (200, 50))).astype(np.float64)),
    ("denso 200x50",   (rng.poisson(6.0, (200, 50))).astype(np.float64)),
    ("misto 300x40",   (rng.poisson(rng.gamma(0.5, 4.0, 40)[None, :], (300, 40))).astype(np.float64)),
]:
    xd[:, 0] = 0.0                    # gene totalmente nulo
    xd[xd.sum(axis=1) == 0, 1] = 1.0  # nenhuma celula vazia
    a, loo = naive(xd)
    b = fechada(xd)
    # bloco menor que BLOCO_GENES e tambem um caso com varios blocos
    err = np.abs(a - b)
    rel = err / np.maximum(a, 1e-300)
    print(f"{nome}: max|dif|={err.max():.3e}  max rel={np.nanmax(rel[a>0]):.3e}  "
          f"se_naive[0]={a[0]:.6e} se_fech[0]={b[0]:.6e}")

# forca varios blocos de genes
cb.BLOCO_GENES = 7
xd = rng.poisson(1.5, (150, 33)).astype(np.float64)
a, _ = naive(xd); b = fechada(xd)
print(f"multiblocos (BLOCO_GENES=7): max rel={np.nanmax(np.abs(a-b)/np.maximum(a,1e-300)):.3e}")
