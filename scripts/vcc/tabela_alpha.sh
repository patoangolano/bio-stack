#!/usr/bin/env bash
# Monta a tabela dos seis membros por braco de alpha, lendo os scored_a*.csv do WSL.
#
# Espelha o `tabela_k.sh`, com uma diferenca que importa na leitura: a referencia
# aqui nao e um registro antigo, e o proprio braco `alpha = 0`. Ele reproduz o
# escalar `k=4` bit a bit por construcao (c=1, fator=1), entao a coluna dele e o
# zero da comparacao e qualquer movimento nas outras e efeito de alpha e de nada
# mais. Por isso a ultima linha e a DIFERENCA contra alpha=0, e nao um valor.
set -uo pipefail
O="${O:-$HOME/vcc-avaliador/out-vcc-D-orion-k40-s0}"
python3 - "$O" <<'PY'
import csv
import pathlib
import sys

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

bracos = {}
for f in sorted(O.glob("scored_a*.csv")):
    # `scored_a075` -> 0,75. Os rotulos vao sem ponto no nome do arquivo porque
    # um ponto no meio do rotulo confundiria quem le o diretorio.
    alpha = int(f.stem.removeprefix("scored_a")) / 100.0
    with f.open(encoding="utf-8") as fh:
        bracos[alpha] = {
            r["metric"]: (
                float(r["from_replicate"]) if r["from_replicate"] not in ("", "nan") else None
            )
            for r in csv.DictReader(fh)
        }

if not bracos:
    sys.exit(f"nenhum scored_a*.csv em {O}")
alphas = sorted(bracos)
base = bracos.get(0.0)

print("membro            " + "".join(f"{'a=' + f'{a:g}':>10s}" for a in alphas))
for chave, rotulo in MEMBROS:
    print(
        f"{rotulo:18s}"
        + "".join(
            (f"{bracos[a][chave]:10.4f}" if bracos[a].get(chave) is not None else f"{'—':>10s}")
            for a in alphas
        )
    )
if base is not None:
    print(
        f"{'delta vs a=0':18s}"
        + "".join(
            (
                f"{bracos[a]['avg_score'] - base['avg_score']:+10.4f}"
                if bracos[a].get("avg_score") is not None and base.get("avg_score") is not None
                else f"{'—':>10s}"
            )
            for a in alphas
        )
    )

melhor = max(alphas, key=lambda a: bracos[a]["avg_score"] or -9)
print(f"\notimo medido: alpha={melhor:g}  avg_score={bracos[melhor]['avg_score']:+.4f}")
if base is not None and base.get("avg_score") is not None:
    ganho = bracos[melhor]["avg_score"] - base["avg_score"]
    print(f"ganho sobre o escalar (alpha=0): {ganho:+.4f}")
    # Por membro: onde o ganho veio, e a quem custou. A curva de k mostrou que os
    # membros discordam, entao um ganho medio pode esconder uma troca.
    print("\npor membro, alpha otimo contra alpha=0:")
    for chave, rotulo in MEMBROS[:-1]:
        a, b = base.get(chave), bracos[melhor].get(chave)
        if a is None or b is None:
            continue
        print(f"  {rotulo:12s} {a:7.4f} -> {b:7.4f}  {b - a:+.4f}")
PY
