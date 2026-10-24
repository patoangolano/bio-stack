"""delta_kd1.npz reproduz delta_previsto.npz? Se nao, k=1 nao e a v4 e a
varredura mede o erro da decomposicao em vez de medir k."""

import numpy as np

R = "C:/Users/Manec/bio-stack/runs/vcc/transfer"
with np.load(f"{R}/delta_previsto.npz", allow_pickle=True) as a, \
        np.load(f"{R}/varredura_k/delta_kd1.npz", allow_pickle=True) as b:
    alvos_a = np.asarray(a["alvos"]).astype(str)
    alvos_b = np.asarray(b["alvos"]).astype(str)
    print("ordem de alvos identica:", bool(np.array_equal(alvos_a, alvos_b)))
    pos = {t: i for i, t in enumerate(alvos_b)}
    reindex = np.array([pos[t] for t in alvos_a])
    for ctx in ("D", "E", "F"):
        x = np.asarray(a[f"delta_{ctx}"], dtype=np.float64)
        y = np.asarray(b[f"delta_{ctx}"], dtype=np.float64)[reindex]
        d = np.abs(x - y)
        rel = d / np.maximum(np.abs(x), 1e-12)
        print(f"{ctx}: max|dif|={d.max():.3e}  p99.9|dif|={np.quantile(d, 0.999):.3e}  "
              f"rms={np.sqrt((d ** 2).mean()):.3e}  "
              f"max|x|={np.abs(x).max():.4f}  rms|x|={np.sqrt((x ** 2).mean()):.4f}")
        print(f"   fracao de entradas com |dif| > 1e-6: {(d > 1e-6).mean():.6f}"
              f"   > 1e-4: {(d > 1e-4).mean():.6f}")
        print(f"   corr(x,y) = {np.corrcoef(x.ravel(), y.ravel())[0, 1]:.9f}"
              f"   max rel (onde |x|>1e-3) = "
              f"{rel[np.abs(x) > 1e-3].max() if (np.abs(x) > 1e-3).any() else float('nan'):.3e}")
