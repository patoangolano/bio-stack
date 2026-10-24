#!/usr/bin/env bash
# Monta a tabela dos seis membros por braco de k, lendo os scored_*.csv do WSL.
set -uo pipefail
O="$HOME/vcc-avaliador/out-vcc-D-orion-k40-s0"
python3 - "$O" <<'PY'
import csv, pathlib, sys
O = pathlib.Path(sys.argv[1])
MEMBROS = [
    ("de_wilcoxon_direction_fidelity_yield_raw", "fidelity"),
    ("de_wilcoxon_sig_jaccard", "jaccard"),
    ("de_wilcoxon_direction_reach_raw", "reach"),
    ("de_wilcoxon_lfc_nmae", "lfc_nmae"),
    ("expr_mse_unbiased_capped_norm", "expr_mse"),
    ("pds_cosine", "pds_cosine"),
    ("avg_score", "AVG_SCORE"),
]
REGISTRADO = {1: 0.275, 2: 0.460, 3: 0.590, 4: 0.647, 5: 0.627, 6: 0.587}
arms = {}
for f in sorted(O.glob("scored_k*.csv")):
    k = int(f.stem.removeprefix("scored_k"))
    with f.open(encoding="utf-8") as fh:
        arms[k] = {r["metric"]: (float(r["from_replicate"]) if r["from_replicate"] not in ("", "nan") else None)
                   for r in csv.DictReader(fh)}
ks = sorted(arms)
print("braco             " + "".join(f"{'k=' + str(k):>10s}" for k in ks))
for chave, rot in MEMBROS:
    print(f"{rot:18s}" + "".join(
        (f"{arms[k][chave]:10.4f}" if arms[k].get(chave) is not None else f"{'—':>10s}") for k in ks))
print(f"{'registrado 22/10':18s}" + "".join(f"{REGISTRADO.get(k, float('nan')):10.3f}" for k in ks))
melhor = max(ks, key=lambda k: arms[k]["avg_score"] or -9)
print(f"\notimo medido: k={melhor}  avg_score={arms[melhor]['avg_score']:+.4f}")
for k in ks:
    v, r = arms[k]["avg_score"], REGISTRADO.get(k)
    if v is not None and r is not None:
        print(f"  k={k}: medido {v:+.4f}  registrado {r:+.3f}  diferenca {v - r:+.4f}")
PY
