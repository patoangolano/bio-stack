"""Figuras do paper das reliquias de esmalte (docs/paper_reliquias.md).

Todo numero vem dos TSV das execucoes, nao de digitacao:
  Fig. 1  sintenia NRG1..PSD3 (coordenadas medidas no s24/A8)
  Fig. 2  indice de retencao corrigido pelo nulo, duas consultas (s33_indice.tsv)
  Fig. 3  lesoes por mil aa, com dente x sem dente (s29_lesoes.tsv)
  Fig. 4  lesoes do ENAM compartilhadas, mapeadas na arvore das aves
          (s31_ENAM_compartilhadas.tsv)

Cor (skill dataviz, paleta de referencia validada com validate_palette.js,
--pairs all, modo claro): azul #2a78d6 = com dente, laranja #eb6834 = sem dente.
Codificacao secundaria por forma do marcador (cheio x vazado) para impressao em
escala de cinza. Texto sempre em tinta de texto, nunca na cor da serie.

    uv run python scripts/paper_figuras.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
RUNS = RAIZ / "runs" / "atavismo"
OUT = RAIZ / "runs" / "relatorios" / "paper_fig"

COM, SEM = "#2a78d6", "#eb6834"
TINTA, TINTA2, GRADE, SUPERF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
NEUTRO = "#9a9890"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.edgecolor": TINTA2,
        "axes.labelcolor": TINTA,
        "xtick.color": TINTA2,
        "ytick.color": TINTA2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": SUPERF,
        "axes.facecolor": SUPERF,
        "savefig.facecolor": SUPERF,
        "axes.linewidth": 0.8,
    }
)

GRUPO = {
    "Amis": "Crocodylia",
    "Ggan": "Crocodylia",
    "Cpor": "Crocodylia",
    "Acar": "Squamata",
    "Pbiv": "Squamata",
    "Pmur": "Squamata",
    "Cpic": "Testudines",
    "Cmyd": "Testudines",
    "Gevg": "Testudines",
    "Psin": "Testudines",
}
NOME = {
    "Amis": "Alligator",
    "Ggan": "Gavialis",
    "Cpor": "Crocodylus",
    "Acar": "Anolis",
    "Pbiv": "Python",
    "Pmur": "Podarcis",
    "Cpic": "Chrysemys",
    "Cmyd": "Chelonia",
    "Gevg": "Gopherus",
    "Psin": "Pelodiscus",
    "Ggal": "Gallus",
    "Shab": "Strigops",
    "Apla": "Anas",
    "Tgut": "Taeniopygia",
    "Dnov": "Dromaius",
    "Scam": "Struthio",
    "Arow": "Apteryx",
    "Tgua": "Tinamus",
}
ORDEM_GRUPOS = ["Crocodylia", "Squamata", "Testudines", "Aves"]


def le(caminho: Path) -> list[dict]:
    with caminho.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def salva(fig: Figure, nome: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / nome, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  {OUT / nome}")


def fig1_sintenia() -> None:
    """Dois tracos, escala propria cada um, ancoras ligadas."""
    jac = {
        "NRG1": (156_331_250, 156_595_276),
        "ENAM": (157_293_890, 157_307_957),
        "AMBN": (157_326_619, 157_341_756),
        "AMTN": (157_359_292, 157_374_521),
        "LPL": (157_488_689, 157_512_736),
        "PSD3": (157_601_092, 157_780_011),
    }
    fra = {
        "NRG1": (54_642_186, 55_109_784),
        "LPL": (55_214_993, 55_231_804),
        "PSD3": (55_231_641, 55_383_432),
    }
    rel = {"ENAM": (55_155_480, 55_157_241), "AMBN": (55_168_855, 55_168_992)}
    fig, ax = plt.subplots(figsize=(7.0, 2.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.25, 1.25)
    ax.axis("off")

    def x(pos: float, lo: float, hi: float) -> float:
        return 0.04 + 0.92 * (pos - lo) / (hi - lo)

    tracos = [
        (1.0, "Alligator mississippiensis  chr3 (NC_081826.1)", jac, 156_300_000, 157_800_000),
        (0.15, "Gallus gallus  chrZ (NC_052572.1, GRCg7b)", fra, 54_620_000, 55_400_000),
    ]
    pos_x: dict[tuple[float, str], float] = {}
    for y, rotulo, genes, lo, hi in tracos:
        ax.plot([0.04, 0.96], [y, y], color=NEUTRO, lw=1.2, zorder=1)
        ax.text(0.04, y + 0.17, rotulo, fontsize=8, color=TINTA, style="italic")
        cluster = []
        for g, (a, b) in genes.items():
            esmalte = g in ("ENAM", "AMBN", "AMTN")
            xa, xb = x(a, lo, hi), x(b, lo, hi)
            larg = max(xb - xa, 0.012)
            ax.add_patch(
                FancyBboxPatch(
                    (xa, y - 0.06),
                    larg,
                    0.12,
                    boxstyle="round,pad=0,rounding_size=0.01",
                    fc=COM if esmalte else "#d9d7d0",
                    ec="none",
                    zorder=2,
                )
            )
            if esmalte:  # os tres estao a ~1 mm um do outro: um rotulo so para o cluster
                cluster.append(xa + larg / 2)
            else:
                ax.text(
                    xa + larg / 2, y - 0.16, g, ha="center", va="top", fontsize=7.5, color=TINTA
                )
            pos_x[(y, g)] = xa + larg / 2
        if cluster:
            ax.text(
                sum(cluster) / len(cluster),
                y + 0.10,
                "ENAM · AMBN · AMTN",
                ha="center",
                va="bottom",
                fontsize=7.5,
                color=TINTA,
                fontweight="bold",
            )
    xs = []
    for _g, (a, b) in rel.items():
        xa = x((a + b) / 2, 54_620_000, 55_400_000)
        xs.append(xa)
        ax.plot([xa], [0.15], marker="v", ms=8, mfc="none", mec=SEM, mew=1.6, zorder=3)
    ax.text(
        sum(xs) / len(xs),
        0.30,
        "ΨENAM · ΨAMBN",
        ha="center",
        fontsize=7.5,
        color=TINTA,
        fontweight="bold",
    )
    for g in ("NRG1", "LPL", "PSD3"):
        ax.plot(
            [pos_x[(1.0, g)], pos_x[(0.15, g)]],
            [0.92, 0.23],
            color=NEUTRO,
            lw=0.7,
            ls=(0, (2, 2)),
            zorder=0,
        )
    legenda = [
        Line2D(
            [],
            [],
            marker="s",
            ls="",
            mfc=COM,
            mec="none",
            ms=7,
            label="gene do esmalte (com dente)",
        ),
        Line2D(
            [],
            [],
            marker="v",
            ls="",
            mfc="none",
            mec=SEM,
            mew=1.6,
            ms=7,
            label="relíquia detectada (sem dente)",
        ),
    ]
    ax.legend(
        handles=legenda, loc="lower right", frameon=False, fontsize=7.5, bbox_to_anchor=(1.0, -0.22)
    )
    salva(fig, "fig1_sintenia.png")


def fig2_indice() -> None:
    """Pirulito horizontal: uma linha por especie, agrupada por clado; dois paineis
    (uma consulta cada), cada um com o seu eixo — sao consultas diferentes."""
    linhas = le(RUNS / "indice_nulo" / "s33_indice.tsv")
    ordem = [
        "Amis",
        "Ggan",
        "Cpor",
        "Acar",
        "Pbiv",
        "Pmur",
        "Cpic",
        "Cmyd",
        "Gevg",
        "Psin",
        "Ggal",
        "Shab",
        "Apla",
        "Tgut",
        "Dnov",
        "Scam",
        "Arow",
        "Tgua",
    ]
    ypos = {e: len(ordem) - 1 - i for i, e in enumerate(ordem)}
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.6), sharey=True)
    for ax, (fonte, titulo) in zip(
        axes,
        [
            ("Amis", "a  consulta: proteínas do jacaré"),
            ("Acar", "b  consulta: proteínas do Anolis"),
        ],
        strict=True,
    ):
        dados = {r["especie"]: r for r in linhas if r["consulta"] == fonte}
        xmax = max(float(r["indice_nulo"]) for e, r in dados.items() if e != fonte)
        for e in ordem:
            y = ypos[e]
            if e == fonte:
                ax.text(0, y, "  (a própria consulta)", va="center", fontsize=6.5, color=TINTA2)
                continue
            r = dados[e]
            v = float(r["indice_nulo"])
            com = r["dente"] == "COM"
            cor = COM if com else SEM
            ax.plot([0, v], [y, y], color=cor, lw=1.6, solid_capstyle="round", zorder=2)
            ax.plot(
                v,
                y,
                marker="o",
                ms=4.8,
                ls="",
                mfc=cor if com else SUPERF,
                mec=cor,
                mew=1.3,
                zorder=3,
            )
            ax.text(
                v + xmax * 0.025,
                y,
                f"{v:.2f}".replace(".", ","),
                va="center",
                fontsize=6.3,
                color=TINTA2,
            )
        for corte in (2.5, 5.5, 9.5):  # separa Crocodylia | Squamata | Testudines | Aves
            ax.axhline(len(ordem) - 1 - corte, color=GRADE, lw=0.8)
        ax.set_xlim(0, xmax * 1.22)
        ax.set_title(titulo, loc="left", fontsize=8.5, color=TINTA)
        ax.set_xlabel("índice de retenção (excesso sobre o nulo)")
        ax.grid(axis="x", color=GRADE, lw=0.5)
        ax.set_axisbelow(True)
    axes[0].set_yticks(
        [ypos[e] for e in ordem], [NOME[e] for e in ordem], fontsize=7, style="italic"
    )
    for nome, meio in (("Crocodylia", 1), ("Squamata", 4), ("Testudines", 7.5), ("Aves", 13.5)):
        axes[1].text(
            1.02,
            (len(ordem) - 1 - meio) / (len(ordem) - 1) * 0.96 + 0.02,
            nome,
            transform=axes[1].transAxes,
            fontsize=7,
            color=TINTA2,
            va="center",
            rotation=90,
        )
    legenda = [
        Line2D([], [], marker="o", ls="", mfc=COM, mec=COM, ms=5, label="com dente"),
        Line2D([], [], marker="o", ls="", mfc=SUPERF, mec=SEM, mew=1.3, ms=5, label="sem dente"),
    ]
    axes[1].legend(handles=legenda, loc="lower right", frameon=False, fontsize=7.5)
    fig.tight_layout()
    salva(fig, "fig2_indice.png")


def fig3_lesoes() -> None:
    # Com dente: TODO alinhamento, como no total do texto (12.160 aa). Sem dente: so
    # as reliquias com lesao — as demais nao tem o que mostrar numa taxa.
    linhas = [r for r in le(RUNS / "lesoes" / "s29_lesoes.tsv") if r["alinhou"] == "True"]
    for r in linhas:
        n = int(r["q_fim"]) - int(r["q_ini"]) + 1
        r["n_aa"] = n
        r["les"] = int(r["frameshifts"]) + int(r["stops"])
        r["taxa"] = 1000 * r["les"] / n
    com = sorted(
        (r for r in linhas if r["dente"] == "COM"), key=lambda r: (r["gene"], r["especie"])
    )
    sem = sorted(
        (r for r in linhas if r["dente"] == "sem" and r["les"] > 0), key=lambda r: -r["taxa"]
    )
    filas = sem + com
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for i, r in enumerate(filas):
        cor = COM if r["dente"] == "COM" else SEM
        y = len(filas) - 1 - i
        ax.plot([0, r["taxa"]], [y, y], color=cor, lw=2, solid_capstyle="round", zorder=2)
        ax.plot(
            r["taxa"],
            y,
            marker="o",
            ms=4.5,
            ls="",
            mfc=cor if r["dente"] == "COM" else "none",
            mec=cor,
            mew=1.3,
            zorder=3,
        )
        ax.text(
            r["taxa"] + 0.5,
            y,
            f"{r['les']} / {r['n_aa']} aa",
            fontsize=6.5,
            va="center",
            color=TINTA2,
        )
    ax.set_yticks(
        range(len(filas)),
        [f"{r['gene']} · {NOME[r['especie']]}" for r in reversed(filas)],
        fontsize=7,
    )
    ax.axhline(len(com) - 0.5, color=GRADE, lw=0.8)
    ax.text(
        ax.get_xlim()[1] if ax.get_xlim()[1] > 1 else 25,
        len(com) - 0.2,
        "sem dente ↑   com dente ↓",
        ha="right",
        fontsize=7,
        color=TINTA2,
    )
    ax.set_xlabel("lesões (frameshift + stop) por mil aminoácidos alinhados")
    ax.grid(axis="x", color=GRADE, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(r["taxa"] for r in filas) * 1.25)
    tot_com = sum(r["les"] for r in com), sum(r["n_aa"] for r in com)
    aa = f"{tot_com[1]:,}".replace(",", ".")
    taxa = f"{1000 * tot_com[0] / tot_com[1]:.2f}".replace(".", ",")
    ax.set_title(
        f"com dente: {tot_com[0]} lesão em {aa} aa ({taxa} por mil)",
        loc="left",
        fontsize=8,
        color=TINTA,
    )
    salva(fig, "fig3_lesoes.png")


def fig4_arvore() -> None:
    linhas = le(RUNS / "nivel4" / "s31_ENAM_compartilhadas.tsv")
    cont: Counter = Counter()
    for r in linhas:
        cont[(r["leitura"].replace("data ", ""), r["tipo"])] += 1
    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    ax.axis("off")
    folhas = ["Apteryx", "Dromaius", "Struthio", "Strigops", "Taeniopygia", "Anas"]
    ys: dict[str, float] = {f: 5.0 - i for i, f in enumerate(folhas)}
    ys["Crocodylia"] = -1.0
    xf = 0.92
    no = {"Novaeratitae": (0.62, (ys["Apteryx"] + ys["Dromaius"]) / 2)}
    no["Palaeognathae"] = (0.45, (no["Novaeratitae"][1] + ys["Struthio"]) / 2)
    no["Psittacopasserae"] = (0.62, (ys["Strigops"] + ys["Taeniopygia"]) / 2)
    no["Neognathae"] = (0.45, (no["Psittacopasserae"][1] + ys["Anas"]) / 2)
    no["Neornithes"] = (0.28, (no["Palaeognathae"][1] + no["Neognathae"][1]) / 2)
    no["raiz"] = (0.10, (no["Neornithes"][1] + ys["Crocodylia"]) / 2)
    filhos = {
        "Novaeratitae": ["Apteryx", "Dromaius"],
        "Palaeognathae": ["Novaeratitae", "Struthio"],
        "Psittacopasserae": ["Strigops", "Taeniopygia"],
        "Neognathae": ["Psittacopasserae", "Anas"],
        "Neornithes": ["Palaeognathae", "Neognathae"],
        "raiz": ["Neornithes", "Crocodylia"],
    }

    def pos(n: str) -> tuple[float, float]:
        return no[n] if n in no else (xf, ys[n])

    for pai, fs in filhos.items():
        xp, yp = pos(pai)
        for f in fs:
            xc, yc = pos(f)
            cor = COM if f == "Crocodylia" else TINTA2
            ax.plot([xp, xp, xc], [yp, yc, yc], color=cor, lw=1.1)
    for f in folhas + ["Crocodylia"]:
        ax.text(
            xf + 0.01,
            ys[f],
            f + ("  (com dente: 0 lesões)" if f == "Crocodylia" else ""),
            va="center",
            fontsize=8,
            style="italic",
            color=TINTA,
        )
    rotulos = {
        "Palaeognathae": "Palaeognathae",
        "Novaeratitae (Apteryx+Dromaius)": "Novaeratitae",
        "Psittacopasserae (Strigops+Taeniopygia)": "Psittacopasserae",
    }
    for chave, n in rotulos.items():
        stops, dels = cont.get((chave, "stop"), 0), cont.get((chave, "del"), 0)
        partes = [f"{stops} stop{'s' if stops != 1 else ''}"] + (
            [f"{dels} del. 10 nt"] if dels else []
        )
        xp, yp = pos(n)
        ax.plot(xp, yp, marker="o", ms=7, mfc="none", mec=SEM, mew=1.6, zorder=3)
        ax.text(
            xp - 0.01,
            yp + 0.25,
            f"{n}\n" + " + ".join(partes),
            ha="right",
            fontsize=6.8,
            color=TINTA,
        )
    xp, yp = pos("Neornithes")
    inc = sum(1 for r in linhas if r["leitura"].startswith("INCONGRUENTE"))
    ax.text(
        xp - 0.01,
        yp + 0.25,
        f"Neornithes\n0 congruentes ({inc} homoplasias)",
        ha="right",
        fontsize=6.8,
        color=TINTA2,
    )
    ax.set_xlim(-0.05, 1.25)
    ax.set_ylim(-1.6, 5.6)
    salva(fig, "fig4_arvore.png")


def main() -> int:
    print("figuras:")
    fig1_sintenia()
    fig2_indice()
    fig3_lesoes()
    fig4_arvore()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
