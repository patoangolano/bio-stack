"""
Figuras: frequencias de quarteto das tres resolucoes de cada ramo contestado,
por tipo de dado e estimador do B10K, e a fracao de quartetos de cada
arvore-fonte da Open Tree para a base de Afroaves.

    uv run python scripts/rapinantes/reanalise/figura_quartetos.py

Le b10k/quartetos_ramos.tsv e opentree/votos_por_arvore.tsv; grava
figuras/quartetos_ramos.png|svg e figuras/opentree_afroaves.png|svg.
A ordem e as cores das resolucoes sao fixas por ramo (a de referencia, da arvore
63K, sempre primeiro e em verde), na mesma ordem de tabelas_md.py.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tabelas_md import PERG, TIPOS  # noqa: E402

RAIZ = Path(__file__).resolve().parents[3]
R = RAIZ / "runs" / "rapinantes" / "reanalise"
FIG = R / "figuras"
CORES = ["#1b9e77", "#d95f02", "#7570b3"]
PAINEIS = ["afroaves_base", "afroaves_mono_se_AS", "australaves_base", "psittacopasserae",
           "accipitrimorphae", "accipitriformes_ss"]
# Stiller et al. 2024, dados-fonte da Fig. 2a, no 26 (63K, 363 taxa, ASTRAL-MP)
ARTIGO = {"afroaves_base": (0.3495, 0.3456, 0.3049)}


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    linhas = list(csv.DictReader(open(R / "b10k" / "quartetos_ramos.tsv", encoding="utf-8"), delimiter="\t"))
    fig, eixos = plt.subplots(3, 2, figsize=(15, 13))
    for ax, p in zip(eixos.flat, PAINEIS, strict=True):
        titulo, chaves, curtos = PERG[p]
        barras: list[tuple[str, list[float]]] = []
        if p in ARTIGO:
            barras.append(("artigo: 63K, 363 táxons, ASTRAL-MP (Fig. 2a)", list(ARTIGO[p])))
        for tipo, nome in TIPOS.items():
            for x in linhas:
                if x["pergunta"] != p or not x.get("EN") or not x["analise"].startswith(tipo + " |"):
                    continue
                desc = x["analise"].split(" | ")[1]
                if "busca livre" in desc and "wASTRAL" not in desc:
                    continue
                est = "wASTRAL aLRT" if "wASTRAL" in desc else "ASTRAL-IV"
                tax = "363" if desc.startswith("363") else "42"
                qs = [float(x.get(f"q[{k}]") or "nan") for k in chaves]
                barras.append((f"{nome}, {tax} táxons, {est} (EN {float(x['EN']):,.0f})".replace(",", " "), qs))
        for i, (_rot, qs) in enumerate(barras):
            esq = 0.0
            for j, v in enumerate(qs):
                ax.barh(i, v, left=esq, color=CORES[j], edgecolor="white", label=curtos[j] if i == 0 else None)
                if v >= 0.06:  # segmento estreito demais para o numero
                    ax.text(esq + v / 2, i, f"{v:.3f}".replace(".", ","), ha="center", va="center", fontsize=6.5,
                            color="white")
                esq += v
        ax.axvline(1 / 3, color="k", ls=":", lw=1)
        ax.set_yticks(range(len(barras)))
        ax.set_yticklabels([b[0] for b in barras], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("frequência de quarteto (as três somam 1)")
        ax.set_title(titulo, fontsize=10)
        ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False)
    fig.suptitle("B10K (Stiller et al. 2024): quartetos em volta dos ramos das aves de rapina\n"
                 "árvores de genes colapsadas em aLRT < 0,95 (ASTRAL-IV, pontuando a topologia 63K) ou ponderadas "
                 "pelo aLRT (wASTRAL, busca livre); linha pontilhada = 1/3; EN = número efetivo de árvores de genes",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"quartetos_ramos.{ext}", dpi=180)
    plt.close(fig)

    votos = list(csv.DictReader(open(R / "opentree" / "votos_por_arvore.tsv", encoding="utf-8"), delimiter="\t"))
    sel = [v for v in votos if "=" in v["afroaves_base_quartetos"]]
    sel.sort(key=lambda v: (int(v["ano"] or 0), v["arvore"]))
    fig, ax = plt.subplots(figsize=(9, 0.32 * len(sel) + 1.8))
    nomes = ["Accipitrimorphae+Strigiformes", "Strigiformes+Coraciimorphae", "Accipitrimorphae+Coraciimorphae"]
    curtos = ["A+S", "S+C", "A+C"]
    for i, v in enumerate(sel):
        d = dict(kv.split("=") for kv in v["afroaves_base_quartetos"].split(";"))
        esq = 0.0
        for j, a in enumerate(nomes):
            val = float(d.get(a, 0))
            ax.barh(i, val, left=esq, color=CORES[j], edgecolor="white", label=curtos[j] if i == 0 else None)
            esq += val
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels([f"{v['primeiro_autor'][:18]} {v['ano']} ({v['arvore']}, {v['folhas']} folhas)"
                        for v in sel], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fração de 2 000 quartetos amostrados (A, S, C, externo) por árvore-fonte")
    ax.set_title(f"Open Tree 16.1: base de Afroaves nas árvores-fonte (n = {len(sel)})\n"
                 "A = Accipitrimorphae, S = Strigiformes, C = Coraciimorphae", fontsize=10)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=3, frameon=False)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIG / f"opentree_afroaves.{ext}", dpi=180)
    plt.close(fig)
    print("figuras ok:", FIG)


if __name__ == "__main__":
    main()
