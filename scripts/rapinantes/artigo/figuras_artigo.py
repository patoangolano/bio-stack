# figuras_artigo.py
"""Figuras do artigo ABNT (reports/artigo/figuras/).

Figura 1: cladogramas da supermatriz (ML) e do wASTRAL dos mesmos 300 locos de
aminoácidos, lado a lado, com as linhagens rapinantes destacadas e os nós em disputa
anotados com o suporte. Figura 2: fator de concordância de sítios (sCFL) dos nós
profundos, decomposto em nó e duas alternativas NNI, lido da árvore ML.

Uso: uv run python scripts/rapinantes/artigo/figuras_artigo.py
"""

from __future__ import annotations

from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Phylo import NewickIO
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

RAIZ = Path(__file__).resolve().parents[3]
EXEC = RAIZ / "runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit"
ARVORE_ML = EXEC / "05_concatenado/concatenado.treefile"
ARVORE_WASTRAL = EXEC / "07_arvore_especies/especies.tre"
NOS_CHAVE = RAIZ / "runs/rapinantes/filogenomica/nos_chave_AA.tsv"
SAIDA = RAIZ / "reports/artigo/figuras"

# Paleta categórica de referência (skill dataviz), três primeiras posições: as únicas
# que passam no validador com todos os pares em jogo. Cariamiformes fica em tinta
# primária e é identificado pelo nome da ponta.
COR = {
    "Accipitrimorphae": "#2a78d6",
    "Falconiformes": "#eb6834",
    "Strigiformes": "#1baf7a",
    "Cariamiformes": "#0b0b0b",
}
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
CINZA = "#898781"
CINZA_CLARO = "#c3c2b7"

GRUPOS = {
    "Accipitrimorphae": [
        "Gymnogyps_californianus",
        "Sarcoramphus_papa",
        "Sagittarius_serpentarius",
        "Pandion_haliaetus",
        "Elanus_axillaris",
        "Gypaetus_barbatus",
        "Leptodon_cayanensis",
        "Aegypius_monachus",
        "Pithecophaga_jefferyi",
        "Aquila_chrysaetos",
        "Harpia_harpyja",
        "Astur_gentilis",
        "Buteo_buteo",
        "Haliaeetus_albicilla",
    ],
    "Falconiformes": [
        "Falco_naumanni",
        "Falco_peregrinus",
        "Ibycter_americanus",
        "Herpetotheres_cachinnans",
    ],
    "Strigiformes": ["Tyto_alba", "Athene_noctua", "Strix_aluco", "Ninox_strenua"],
    "Cariamiformes": ["Cariama_cristata"],
}
PSITTACOPASSERAE = [
    "Nestor_notabilis",
    "Strigops_habroptila",
    "Melopsittacus_undulatus",
    "Corvus_cornix",
    "Taeniopygia_guttata",
    "Chiroxiphia_lanceolata",
    "Acanthisitta_chloris",
]
CORACIIMORPHAE = [
    "Colius_striatus",
    "Leptosomus_discolor",
    "Apaloderma_vittatum",
    "Buceros_rhinoceros",
    "Merops_nubicus",
    "Dryobates_pubescens",
    "Indicator_indicator",
    "Pogoniulus_pusillus",
]
ACC, FAL, STR, CAR = (set(GRUPOS[k]) for k in GRUPOS)
PSI, COR_ = set(PSITTACOPASSERAE), set(CORACIIMORPHAE)
TELLURAVES = ACC | FAL | STR | CAR | PSI | COR_

# Nós anotados: (rótulo, conjunto de pontas). Só aparecem se existirem na árvore.
NOS_ANOTADOS = [
    ("Telluraves", TELLURAVES),
    ("Australaves", CAR | FAL | PSI),
    ("Eufalconimorphae", FAL | PSI),
    ("Afroaves", ACC | STR | COR_),
    ("Strigiformes + Coraciimorphae", STR | COR_),
    ("Accipitrimorphae", ACC),
    ("Cariamiformes + Falconiformes", CAR | FAL),
    ("Hieraves", ACC | STR),
    ("rapinantes + $\\it{Leptosomus}$", ACC | STR | FAL | CAR | {"Leptosomus_discolor"}),
]


def carregar(caminho: Path):
    with caminho.open(encoding="utf-8") as fh:
        arvore = next(NewickIO.parse(fh))
    arvore.root_with_outgroup("Gallus_gallus")
    arvore.ladderize(reverse=True)
    return arvore


def grupo_de(nomes: set[str]) -> str | None:
    for grupo, membros in GRUPOS.items():
        if nomes <= set(membros):
            return grupo
    return None


def suporte(clado, tipo: str) -> str:
    """UFBoot da supermatriz (rótulo SH-aLRT/aBayes/UFBoot) ou localPP do wASTRAL."""
    if tipo == "ml":
        rotulo = clado.name or ""
        return rotulo.split("/")[-1] if "/" in rotulo else ""
    if clado.confidence is None:
        return ""
    return f"{clado.confidence:.2f}".replace(".", ",")


def desenhar(ax, arvore, tipo: str, titulo: str) -> None:
    pontas = arvore.get_terminals()
    y: dict = {p: float(i) for i, p in enumerate(pontas)}
    altura: dict = {}

    def calc_altura(c):
        if c.is_terminal():
            altura[c] = 0
        else:
            altura[c] = 1 + max(calc_altura(f) for f in c.clades)
        return altura[c]

    calc_altura(arvore.root)
    h_raiz = altura[arvore.root]

    def calc_y(c):
        if c.is_terminal():
            return y[c]
        ys = [calc_y(f) for f in c.clades]
        y[c] = (min(ys) + max(ys)) / 2
        return y[c]

    calc_y(arvore.root)
    x = {c: h_raiz - altura[c] for c in y}
    nomes = {c: {t.name for t in c.get_terminals()} for c in y}

    def cor(c) -> tuple[str, float]:
        g = grupo_de(nomes[c])
        if g:
            return COR[g], 1.3
        if nomes[c] <= TELLURAVES:
            return TINTA_2, 0.8
        return CINZA_CLARO, 0.8

    for c in y:
        if c.is_terminal():
            continue
        cc, lw = cor(c)
        ys = [y[f] for f in c.clades]
        ax.plot([x[c], x[c]], [min(ys), max(ys)], color=cc, lw=lw, solid_capstyle="butt")
        for f in c.clades:
            cf, lwf = cor(f)
            ax.plot([x[c], x[f]], [y[f], y[f]], color=cf, lw=lwf, solid_capstyle="butt")

    for p in pontas:
        nome = p.name.replace("_", " ")
        g = grupo_de({p.name})
        ax.text(
            h_raiz + 0.35,
            y[p],
            nome,
            fontsize=5.3,
            style="italic",
            va="center",
            color=TINTA if (g or p.name in TELLURAVES) else CINZA,
        )

    # Rótulos de nó: nome + suporte, logo acima do ramo que chega ao nó.
    alvos = {frozenset(s): rot for rot, s in NOS_ANOTADOS}
    for c in y:
        if c.is_terminal():
            continue
        rot = alvos.get(frozenset(nomes[c]))
        if rot is None:
            continue
        s = suporte(c, tipo)
        ax.plot(x[c], y[c], "o", ms=3.2, color=TINTA, zorder=5)
        ax.text(
            x[c] - 0.15,
            y[c] - 0.35,
            f"{rot} {s}".strip(),
            fontsize=5.4,
            ha="right",
            va="bottom",
            color=TINTA,
            zorder=6,
            bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85),
        )

    # largura do painel em unidades da árvore: 20 % para rótulos de nó, 45 % para a árvore, 35 % para os nomes
    largura = h_raiz / 0.45
    ax.set_xlim(-0.20 * largura, h_raiz + 0.35 * largura)
    ax.set_ylim(len(pontas) - 0.3, -1.2)
    ax.axis("off")
    ax.set_title(titulo, fontsize=7.5, loc="left", color=TINTA)


def figura1() -> Path:
    fig, eixos = plt.subplots(1, 2, figsize=(6.3, 8.6), dpi=300)
    desenhar(eixos[0], carregar(ARVORE_ML), "ml", "A  Supermatriz particionada (ML; UFBoot)")
    desenhar(
        eixos[1],
        carregar(ARVORE_WASTRAL),
        "wastral",
        "B  wASTRAL das 300 árvores de genes (localPP)",
    )
    marcas = [
        Line2D([], [], color=COR[g], lw=2)
        for g in ("Accipitrimorphae", "Falconiformes", "Strigiformes")
    ]
    marcas += [Line2D([], [], color=TINTA_2, lw=1), Line2D([], [], color=CINZA_CLARO, lw=1)]
    fig.legend(
        marcas,
        [
            "Accipitrimorphae",
            "Falconiformes",
            "Strigiformes",
            "outros Telluraves (inclui $\\it{Cariama}$)",
            "fora de Telluraves",
        ],
        loc="lower center",
        ncol=3,
        fontsize=6,
        frameon=False,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.05, wspace=0.04)
    destino = SAIDA / "figura1_arvores.png"
    fig.savefig(destino, dpi=300)
    plt.close(fig)
    return destino


# Nós da Figura 2, na ordem da árvore; os nomes das alternativas só onde foram
# identificadas pelas árvores restritas (H2 e H7b, seção 4.10 dos Métodos).
NOS_FIG2 = [
    ("Telluraves", "Telluraves", None),
    ("Australaves", "Australaves", None),
    ("Eufalconimorphae", "Eufalconimorphae", "Cariamiformes + Falconiformes"),
    ("Psittacopasserae", "Psittacopasserae", None),
    ("Falconiformes", "Falconiformes", None),
    ("Afroaves", "Afroaves", None),
    ("Strigiformes + Coraciimorphae (H1/H3)", "Strigiformes + Coraciimorphae", "Hieraves"),
    ("Accipitrimorphae", "Accipitrimorphae", None),
    ("Pandion + Accipitridae", "$\\it{Pandion}$ + Accipitridae", None),
    ("Strigiformes", "Strigiformes", None),
]


def formatar_p(p: float) -> str:
    """p com vírgula decimal; abaixo de 0,001 em notação de potência de 10."""
    if p < 0.001:
        mant, exp = f"{p:.1e}".split("e")
        return f"{mant.replace('.', ',')} × $10^{{{int(exp)}}}$"
    return f"{p:.2g}".replace(".", ",")


def um_decimal(v: float) -> str:
    """Arredonda meio para cima, como o texto do artigo (51,65 -> 51,7; o .1f daria 51,6)."""
    return str(Decimal(str(v)).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)).replace(".", ",")


def figura2() -> Path:
    tab = pd.read_csv(NOS_CHAVE, sep="\t")
    tab = tab[tab["conjunto"] == "AA_MAFFT"].set_index("no")
    fig, ax = plt.subplots(figsize=(6.3, 3.6), dpi=300)
    for i, (no, rotulo, alt_nome) in enumerate(NOS_FIG2):
        lin = tab.loc[no]
        alts = sorted([lin["sDF1"], lin["sDF2"]], reverse=True)
        segs = [(lin["sCF"], COR["Accipitrimorphae"]), (alts[0], CINZA), (alts[1], CINZA_CLARO)]
        inicio = 0.0
        for j, (v, c) in enumerate(segs):
            # 1 pt de folga entre segmentos, no lugar do contorno
            ax.barh(i, v - 0.4, left=inicio + 0.2, height=0.62, color=c, linewidth=0)
            cor_txt = "white" if j < 2 else TINTA
            ax.text(
                inicio + v / 2,
                i,
                um_decimal(v),
                ha="center",
                va="center",
                fontsize=6,
                color=cor_txt,
            )
            inicio += v
        if alt_nome:
            ax.text(
                100.8,
                i,
                f"maior alternativa:\n{alt_nome}",
                fontsize=5.4,
                va="center",
                color=TINTA_2,
                linespacing=1.3,
            )
        p = lin["p_sDF1_vs_sDF2"]
        ax.text(
            -1.2,
            i + 0.38,
            f"assimetria: p = {formatar_p(p)}",
            fontsize=4.8,
            ha="right",
            va="center",
            color=CINZA,
        )
        ax.text(-1.2, i - 0.08, rotulo, fontsize=6.4, ha="right", va="center", color=TINTA)
    ax.axvline(100 / 3, color=TINTA, lw=0.8, ls=(0, (3, 2)), zorder=0)
    ax.text(
        100 / 3 + 0.8, -0.62, "1/3 (politomia)", fontsize=5.6, ha="left", va="bottom", color=TINTA
    )
    ax.set_ylim(len(NOS_FIG2) - 0.5, -0.9)
    ax.set_xlim(0, 100)
    ax.set_xlabel(
        "sítios (%) que favorecem cada resolução do quarteto (sCFL)", fontsize=6.6, color=TINTA
    )
    ax.tick_params(axis="x", labelsize=6, colors=TINTA_2, length=2)
    ax.set_yticks([])
    for lado in ("top", "right", "left"):
        ax.spines[lado].set_visible(False)
    ax.spines["bottom"].set_color(CINZA_CLARO)
    marcas = [
        Rectangle((0, 0), 1, 1, color=c) for c in (COR["Accipitrimorphae"], CINZA, CINZA_CLARO)
    ]
    ax.legend(
        marcas,
        [
            "nó da supermatriz (sCF)",
            "alternativa NNI mais frequente",
            "alternativa NNI menos frequente",
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.17),
        ncol=3,
        fontsize=5.8,
        frameon=False,
    )
    fig.subplots_adjust(left=0.25, right=0.77, top=0.93, bottom=0.2)
    destino = SAIDA / "figura2_scf.png"
    fig.savefig(destino, dpi=300)
    plt.close(fig)
    return destino


if __name__ == "__main__":
    SAIDA.mkdir(parents=True, exist_ok=True)
    for f in (figura1(), figura2()):
        print(f.relative_to(RAIZ))
