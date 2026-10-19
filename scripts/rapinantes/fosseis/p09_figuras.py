"""Figuras de resumo: p-AU por hipotese de posicao de fossil e IC por classe de carater."""

from __future__ import annotations

import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import clados as cl
import matplotlib.pyplot as plt

FIG = cl.BASE / "figuras"


def ler(caminho):
    with open(caminho, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def figura_au() -> None:
    linhas = [
        r
        for r in ler(cl.BASE / "au" / "au_resultados.tsv")
        if not r["fossil"].startswith("atuais_")
    ]
    grupos: dict[tuple[str, str], list[tuple[str, float]]] = defaultdict(list)
    for r in linhas:
        if r["hipotese"] == "ML_fosseis_livres":
            continue
        grupos[(r["matriz"], r["fossil"])].append((r["hipotese"], float(r["p-AU"])))
    n = len(grupos)
    fig, axes = plt.subplots(n, 1, figsize=(7.5, 1.05 * n + 0.6), sharex=True)
    for ax, ((m, f), hs) in zip(axes, sorted(grupos.items()), strict=True):
        hs = sorted(hs, key=lambda x: -x[1])
        xs = [p for _, p in hs]
        ax.scatter(
            xs, [0] * len(xs), c=["#b2182b" if p < 0.05 else "#2166ac" for p in xs], s=22, zorder=3
        )
        for i, (h, p) in enumerate(hs[:4]):
            ax.annotate(
                h.replace("_", "+") if h.count("_") == 1 and h[0].isupper() and "_" in h else h,
                (p, 0),
                xytext=(0, 6 + 7 * (i % 2)),
                textcoords="offset points",
                fontsize=5.5,
                ha="center",
            )
        ax.axvline(0.05, color="grey", lw=0.6, ls="--")
        ax.set_yticks([])
        ax.set_ylabel(
            f"{m}\n{f.replace('_', ' ')[:24]}", rotation=0, ha="right", va="center", fontsize=6.5
        )
        ax.set_xscale("log")
        ax.set_xlim(max(1e-7, 0.5 * min(min(p for _, p in v) for v in grupos.values())), 1.5)
    axes[-1].set_xlabel("p-AU (IQ-TREE 3, 10 000 RELL); vermelho = rejeitada a 5%")
    fig.suptitle(
        "Posicao de fosseis-chave: fossil irmao de X, com esqueleto de Stiller et al. 2024",
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(FIG / "au_fosseis.png", dpi=220)
    plt.close(fig)


def figura_ci() -> None:
    caminho = cl.BASE / "homoplasia" / "caracteres_CI_RI.tsv"
    if not caminho.exists():
        return
    rs = [
        r
        for r in ler(caminho)
        if r["versao"] == "completa"
        and r["informativo"] == "TRUE"
        and r["classe"] not in ("NA", "")
    ]
    ids = sorted({r["matriz"] for r in rs})
    fig, axes = plt.subplots(1, len(ids), figsize=(2.3 * len(ids), 3.2), sharey=True)
    if len(ids) == 1:
        axes = [axes]
    for ax, id_ in zip(axes, ids, strict=True):
        dados, rot = [], []
        for c in ("bico", "pe_garra", "tarsometatarso", "outro"):
            for arv, col in (("morf", "CI_morf"), ("mol", "CI_mol")):
                v = [
                    float(r[col])
                    for r in rs
                    if r["matriz"] == id_ and r["classe"] == c and r[col] != "NA"
                ]
                dados.append(v)
                rot.append(f"{c[:4]}\n{arv}\n(n={len(v)})")
        ax.boxplot(dados, tick_labels=rot, widths=0.6, showfliers=False)
        ax.set_title(id_, fontsize=8)
        ax.tick_params(labelsize=5)
    axes[0].set_ylabel("IC por carater (Fitch)")
    fig.suptitle(
        "Homoplasia por classe de carater: arvore morfologica x molecular-restrita", fontsize=8
    )
    fig.tight_layout()
    fig.savefig(FIG / "homoplasia_classes.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    figura_au()
    figura_ci()
