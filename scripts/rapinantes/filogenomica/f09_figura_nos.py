"""
f09 — figura dos nos-chave: uma linha por no, uma coluna por conjunto/sensibilidade.

Cada celula: gCF / sCF (supermatriz) e, embaixo, a probabilidade local do wASTRAL.
Cor = gCF (0-100). Celula vazia com "ausente": o no nao esta na arvore da
supermatriz daquele conjunto (o PP do wASTRAL aparece se o no estiver la).

    uv run python scripts/rapinantes/filogenomica/f09_figura_nos.py <nos_chave.tsv> <saida_sem_extensao>
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def main() -> None:
    linhas = list(csv.DictReader(Path(sys.argv[1]).open(encoding="utf-8"), delimiter="\t"))
    saida = Path(sys.argv[2])
    conjuntos = list(dict.fromkeys(r["conjunto"] for r in linhas))
    nos = list(dict.fromkeys(r["no"] for r in linhas))
    idx = {(r["conjunto"], r["no"]): r for r in linhas}
    m = np.full((len(nos), len(conjuntos)), np.nan)
    fig, ax = plt.subplots(figsize=(2.1 * len(conjuntos) + 4.5, 0.42 * len(nos) + 1.6))
    for i, no in enumerate(nos):
        for j, cj in enumerate(conjuntos):
            r = idx.get((cj, no))
            if not r:
                continue
            pp = r.get("localPP") or ""
            pp_txt = f"PP {float(pp):.2f}" if pp else ("PP —" if r.get("no_wASTRAL") == "nao" else "")
            if r.get("na_supermatriz") == "sim" and r.get("gCF"):
                m[i, j] = float(r["gCF"])
                txt = f"{float(r['gCF']):.0f} / {float(r['sCF']):.0f}" if r.get("sCF") else f"{float(r['gCF']):.0f}"
            else:
                txt = "ausente"
            ax.text(j, i, f"{txt}\n{pp_txt}", ha="center", va="center", fontsize=7)
    im = ax.imshow(m, cmap="viridis", vmin=0, vmax=100, aspect="auto", alpha=0.75)
    ax.set_xticks(range(len(conjuntos)), conjuntos, rotation=20, ha="right", fontsize=8)
    ax.set_yticks(range(len(nos)), nos, fontsize=8)
    n_genes = {r["conjunto"]: r.get("genes_informativos", "") for r in linhas if r["no"] == nos[0]}
    ax.set_title("Nós-chave: gCF / sCFL na supermatriz (%) e probabilidade local do wASTRAL\n"
                 f"(gCF sobre as árvores de genes de cada conjunto; n de locos informativos no 1º nó: {n_genes})",
                 fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("gCF (%)", fontsize=8)
    fig.tight_layout()
    saida.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(saida.with_suffix(f".{ext}"), dpi=200)
    print(saida.with_suffix(".png"))


if __name__ == "__main__":
    main()
