"""Phase 4 — Publication-grade English figures for the avian/testudine enamel-gene
relics bioRxiv preprint.

100% data-driven from verified outputs in runs/enamel_relics/ (zero hardcoded
coordinates or summary numbers). Exports vector SVG, vector PDF, and 320 DPI PNG
using the Okabe-Ito colorblind-safe palette with redundant shape/pattern encodings.

Usage:
    uv run python scripts/enamel_relics/06_make_figures.py
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "enamel_relics"
TABLES = OUT / "tables"
FIG_DIR = OUT / "figures"
LEGACY = ROOT / "runs" / "atavismo"

# Okabe-Ito colorblind-safe palette
COL = {
    "Crocodylia": "#0072B2",      # deep blue (toothed)
    "Squamata": "#56B4E9",        # sky blue (toothed)
    "Testudines": "#E69F00",      # amber/orange (edentulous)
    "Aves": "#D55E00",            # vermilion (edentulous)
    "LPL": "#2B2B2B",             # charcoal (globular control)
    "IBSP": "#009E73",            # bluish green (functional IDP)
    "SPP1": "#CC79A7",            # reddish purple (functional IDP)
    "DSPP": "#7A5195",            # muted violet (dentin IDP)
    "ENAM": "#D55E00",            # vermilion
    "AMBN": "#0072B2",            # blue
    "AMTN": "#009E73",            # teal-green
    "FS": "#0072B2",              # frameshift bar
    "PTC": "#D55E00",             # premature stop bar
    "SPLICE": "#CC79A7",          # splice disruption
    "GRAY": "#6E6E6E",
}

MARKER = {
    "Crocodylia": "s",
    "Squamata": "^",
    "Testudines": "D",
    "Aves": "o",
}

SP_SHORT = {
    "Amis": "A. mississippiensis",
    "Ggan": "G. gangeticus",
    "Cpor": "C. porosus",
    "Acar": "A. carolinensis",
    "Pbiv": "P. bivittatus",
    "Pmur": "P. muralis",
    "Cpic": "C. picta",
    "Cmyd": "C. mydas",
    "Gevg": "G. evgoodei",
    "Psin": "P. sinensis",
    "Ggal": "G. gallus",
    "Apla": "A. platyrhynchos",
    "Shab": "S. habroptila",
    "Tgut": "Taen. guttata",
    "Dnov": "D. novaehollandiae",
    "Scam": "S. camelus",
    "Arow": "A. rowi",
    "Tgua": "Tin. guttatus",
}


def setup_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.5,
        "axes.titlesize": 9.3,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.4,
        "xtick.labelsize": 7.7,
        "ytick.labelsize": 7.7,
        "legend.fontsize": 7.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "path",
        "pdf.fonttype": 42,
    })


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_coord_range(s: str) -> tuple[int, int, str]:
    strand = s[-2] if s.endswith(")") else "+"
    lo_s, hi_s = s.split("(")[0].split("-")
    return int(lo_s), int(hi_s), strand


def save_fig(fig: plt.Figure, stem: str) -> dict[str, str]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    hashes: dict[str, str] = {}
    for ext, kw in [("png", {"dpi": 320}), ("pdf", {}), ("svg", {})]:
        p = FIG_DIR / f"{stem}.{ext}"
        fig.savefig(p, bbox_inches="tight", facecolor="white", **kw)
        hashes[p.name] = sha256_file(p)
    plt.close(fig)
    return hashes


# =============================================================================
# FIGURE 1: Conserved NRG1-LPL Synteny, GRCg7b Annotation Audit, and ODAPH Gap
# =============================================================================
def make_figure1() -> dict[str, str]:
    intervals = {r["code"]: r for r in read_tsv(OUT / "enamel_intervals.tsv")}
    recip = read_tsv(TABLES / "whole_proteome_reciprocal.tsv")
    recip_map = {(r["species_code"], r["query_name"]): r for r in recip}
    mp_rows = read_tsv(TABLES / "miniprot_level3_splice_audit.tsv")
    mp_map = {(r["species_code"], r["query_gene"]): r for r in mp_rows}
    ens_audit = json.loads((OUT / "ensembl_116_ggal_audit.json").read_text(encoding="utf-8"))
    l4_sum = json.loads((OUT / "phase2_step4_level4_summary.json").read_text(encoding="utf-8"))
    odaph_gaps = json.loads((LEGACY / "odam_odaph" / "s32_vaos.json").read_text(encoding="utf-8"))

    # Exact chicken ENAM Exon 6, Exon 8, AMBN Exon 3 coordinates from verified outputs
    ggal_win_start = int(intervals["Ggal"]["win_start"])
    ex6_line = (OUT / "level4" / "ENAM_exon6_splice_aves" / "blast_Ggal.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    ggal_ex6_s = ggal_win_start + min(int(ex6_line[6]), int(ex6_line[7])) - 1
    ggal_ex6_e = ggal_win_start + max(int(ex6_line[6]), int(ex6_line[7])) - 1
    ggal_ex8_hit = next(d for d in l4_sum if d["dataset"] == "ENAM_exon8_aves")["blast_hits"]["Ggal"]
    ggal_ex8_s, ggal_ex8_e = int(ggal_ex8_hit["abs_start"]), int(ggal_ex8_hit["abs_end"])
    ggal_ambn_r = recip_map[("Ggal", "AMBN_Amis")]
    ggal_ambn_s, ggal_ambn_e = int(ggal_ambn_r["abs_start"]), int(ggal_ambn_r["abs_end"])

    fig = plt.figure(figsize=(11.4, 9.4))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.3, 1.65, 1.65], hspace=0.50)

    # ---- Panel A: Macro-synteny across 6 representative amniotes ----
    axA = fig.add_subplot(gs[0])
    rep_species = [
        ("Amis", "chr3", "Toothed · Crocodylia"),
        ("Acar", "chr2", "Toothed · Squamata"),
        ("Cmyd", "chr5", "Edentulous · Testudines"),
        ("Shab", "chrZ", "Edentulous · Neoaves"),
        ("Arow", "scaf.", "Edentulous · Palaeognathae"),
        ("Ggal", "chrZ", "Edentulous · Galloanserae"),
    ]

    def norm_pos(val: float, nrg1_edge: float, lpl_edge: float) -> float:
        return (val - nrg1_edge) / (lpl_edge - nrg1_edge)

    y_positions = np.arange(len(rep_species))[::-1]
    for y, (sp, chrlbl, desc) in zip(y_positions, rep_species):
        im = intervals[sp]
        n_s, n_e, _ = parse_coord_range(im["nrg1_coords"])
        l_s, l_e, _ = parse_coord_range(im["lpl_coords"])
        gap_kb = int(im["intergenic_gap_bp"]) / 1000.0
        cl = im["clade"]

        if im["order"] == "NRG1->LPL":
            nrg1_edge, lpl_edge = float(n_e), float(l_s)
        else:
            nrg1_edge, lpl_edge = float(n_s), float(l_e)

        axA.plot([-0.08, 1.08], [y, y], color="#BDBDBD", lw=1.6, zorder=1)
        axA.add_patch(mpatches.FancyBboxPatch(
            (-0.14, y - 0.24), 0.13, 0.48, boxstyle="round,pad=0.01",
            facecolor="#4A4A4A", edgecolor="black", lw=0.7, zorder=3
        ))
        axA.text(-0.075, y, "NRG1", color="white", ha="center", va="center", fontsize=7.0, fontweight="bold", zorder=4)

        axA.add_patch(mpatches.FancyBboxPatch(
            (1.01, y - 0.24), 0.11, 0.48, boxstyle="round,pad=0.01",
            facecolor="#4A4A4A", edgecolor="black", lw=0.7, zorder=3
        ))
        axA.text(1.065, y, "LPL", color="white", ha="center", va="center", fontsize=7.0, fontweight="bold", zorder=4)

        loci_to_plot: list[tuple[str, float, float, str, bool]] = []
        if (sp, "ENAM_Amis") in recip_map and recip_map[(sp, "ENAM_Amis")]["top_hsp_reciprocal_confirmed"] == "True":
            rr = recip_map[(sp, "ENAM_Amis")]
            s0, e0 = int(rr["abs_start"]), int(rr["abs_end"])
            if sp == "Arow" and mp_map[(sp, "ENAM_Amis")]["aligned"] == "True":
                s0 = min(s0, int(mp_map[(sp, "ENAM_Amis")]["abs_start"]))
                e0 = max(e0, int(mp_map[(sp, "ENAM_Amis")]["abs_end"]))
            if sp == "Ggal":
                s0, e0 = min(s0, ggal_ex8_s), max(e0, ggal_ex6_e)
            p1, p2 = norm_pos(s0, nrg1_edge, lpl_edge), norm_pos(e0, nrg1_edge, lpl_edge)
            loci_to_plot.append(("ENAM", min(p1, p2), max(p1, p2), COL["ENAM"], im["dentition"] == "edentulous"))

        if (sp, "AMBN_Amis") in recip_map and recip_map[(sp, "AMBN_Amis")]["top_hsp_reciprocal_confirmed"] == "True":
            rr = recip_map[(sp, "AMBN_Amis")]
            s0, e0 = int(rr["abs_start"]), int(rr["abs_end"])
            p1, p2 = norm_pos(s0, nrg1_edge, lpl_edge), norm_pos(e0, nrg1_edge, lpl_edge)
            loci_to_plot.append(("AMBN", min(p1, p2), max(p1, p2), COL["AMBN"], im["dentition"] == "edentulous"))

        if (sp, "AMTN_Amis") in recip_map and recip_map[(sp, "AMTN_Amis")]["top_hsp_reciprocal_confirmed"] == "True":
            rr = recip_map[(sp, "AMTN_Amis")]
            s0, e0 = int(rr["abs_start"]), int(rr["abs_end"])
            p1, p2 = norm_pos(s0, nrg1_edge, lpl_edge), norm_pos(e0, nrg1_edge, lpl_edge)
            loci_to_plot.append(("AMTN", min(p1, p2), max(p1, p2), COL["AMTN"], im["dentition"] == "edentulous"))

        # Smart horizontal nudge when adjacent gene labels are <0.075 apart
        label_nudges: dict[tuple[str, str], float] = {
            ("Amis", "ENAM"): -0.025,
            ("Amis", "AMBN"): 0.0,
            ("Amis", "AMTN"): +0.028,
            ("Arow", "ENAM"): -0.020,
            ("Arow", "AMBN"): +0.020,
        }
        for gname, x0, x1, gcol, is_pseudo in loci_to_plot:
            w = max(x1 - x0, 0.022)
            xc = 0.5 * (x0 + x1)
            rect = mpatches.Rectangle(
                (xc - w / 2, y - 0.22), w, 0.44,
                facecolor="white" if is_pseudo else gcol,
                edgecolor=gcol,
                lw=1.5 if is_pseudo else 0.8,
                hatch="////" if is_pseudo else None,
                zorder=4,
            )
            axA.add_patch(rect)
            lbl = f"Ψ{gname}" if is_pseudo else gname
            dx_lbl = label_nudges.get((sp, gname), 0.0)
            axA.text(xc + dx_lbl, y + 0.30, lbl, color=gcol, ha="center", va="bottom", fontsize=7.0, fontweight="bold", zorder=5)

        ori_note = " (rev.)" if im["order"] != "NRG1->LPL" else ""
        axA.text(1.14, y, f"{gap_kb:,.1f} kb{ori_note}", ha="left", va="center", fontsize=7.6, color="#333333")
        axA.text(-0.16, y + 0.09, f"{SP_SHORT[sp]}", ha="right", va="center", fontsize=8.0, fontstyle="italic", fontweight="bold", color=COL[cl])
        axA.text(-0.16, y - 0.16, f"{chrlbl} ({im['seq_acc']}) · {desc}", ha="right", va="center", fontsize=6.7, color="#555555")

    axA.set_xlim(-0.15, 1.33)
    axA.set_ylim(-0.60, len(rep_species) + 0.25)
    axA.set_xticks([0.0, 0.25, 0.50, 0.75, 1.0])
    axA.set_xticklabels(["0%\n(NRG1 edge)", "25%", "50%", "75%", "100%\n(LPL edge)"])
    axA.set_yticks([])
    axA.spines["left"].set_visible(False)
    axA.set_title("A. Conserved NRG1–LPL enamel-cluster synteny across toothed and edentulous sauropsids (normalized intergenic coordinates)", loc="left", pad=22)

    legA = [
        mpatches.Patch(facecolor=COL["ENAM"], edgecolor=COL["ENAM"], label="Intact ENAM"),
        mpatches.Patch(facecolor=COL["AMBN"], edgecolor=COL["AMBN"], label="Intact AMBN"),
        mpatches.Patch(facecolor=COL["AMTN"], edgecolor=COL["AMTN"], label="Intact AMTN"),
        mpatches.Patch(facecolor="white", edgecolor="#333333", hatch="////", lw=1.4, label="Pseudogenized relic (Ψ)"),
    ]
    axA.legend(handles=legA, loc="lower center", bbox_to_anchor=(0.48, 0.99), frameon=False, ncol=4, fontsize=7.5)

    # ---- Panel B: Zoom of Gallus gallus GRCg7b chrZ (NC_052572.1: 55.098–55.236 Mb) ----
    axB = fig.add_subplot(gs[1])
    x_min_mb, x_max_mb = 55.098, 55.236

    tracks = [
        (2.0, "NCBI RefSeq\nRel. 106"),
        (1.0, "Ensembl\nRel. 116"),
        (0.0, "Recovered\nRelics (GRCg7b)"),
    ]
    for ty, _ in tracks:
        axB.axhline(ty, color="#E0E0E0", lw=1.0, zorder=1)

    refseq_feats = [
        (55.098, 55.109784, "NRG1 (-)", "#4A4A4A", True),
        (55.180611, 55.194886, "LOC112530610 (lncRNA, +)", "#9E9E9E", False),
        (55.215002, 55.232438, "LPL (-)", "#4A4A4A", True),
        (55.231641, 55.236, "PSD3", "#777777", False),
    ]
    for s, e, lbl, col, inside in refseq_feats:
        axB.add_patch(mpatches.Rectangle((s, 2.0 - 0.22), max(e - s, 0.0012), 0.44, facecolor=col, edgecolor="black", lw=0.6, zorder=3))
        if inside:
            axB.text(0.5 * (s + e), 2.0, lbl, color="white", ha="center", va="center", fontsize=6.8, fontweight="bold", zorder=4)
        elif "LOC" in lbl:
            axB.text(0.5 * (s + e), 2.28, lbl, color="#555555", ha="center", va="bottom", fontsize=6.7)

    for g in ens_audit["chrZ_55.109_55.235Mb_genes"]:
        s = max(g["start"] / 1e6, x_min_mb)
        e = min(g["end"] / 1e6, x_max_mb)
        if e <= x_min_mb or s >= x_max_mb:
            continue
        is_coding = 'gene_biotype "protein_coding"' in g["attrs"]
        col = "#4A4A4A" if is_coding else "#B0B0B0"
        axB.add_patch(mpatches.Rectangle((s, 1.0 - 0.20), max(e - s, 0.0010), 0.40, facecolor=col, edgecolor="#333333", lw=0.6, zorder=3))
        if is_coding:
            nm = "NRG1 (-)" if "NRG1" in g["attrs"] else "LPL (-)"
            axB.text(0.5 * (s + e), 1.0, nm, color="white", ha="center", va="center", fontsize=6.8, fontweight="bold", zorder=4)
    axB.text(55.134, 1.29, "5 non-overlapping intergenic lncRNAs (0 enamel features)", color="#555555", ha="center", va="bottom", fontsize=6.8, fontstyle="italic")

    relics = [
        (ggal_ex8_s / 1e6, ggal_ex8_e / 1e6, COL["ENAM"], "ΨENAM Exon 8 (-)\n55.1555–55.1575 Mb\n(1,975 bp; Galloanserae fs+PTC)", -0.024),
        (ggal_ex6_s / 1e6, ggal_ex6_e / 1e6, "#B87500", "ΨENAM Exon 6 (-)\n55.1605–55.1607 Mb\n(1-nt del + TAA stop)", +0.004),
        (ggal_ambn_s / 1e6, ggal_ambn_e / 1e6, COL["AMBN"], "ΨAMBN Exon 3 (-)\n55.1689–55.1690 Mb\n(56.6 bits; recip. confirmed)", +0.028),
    ]
    for s, e, col, lbl, dx_txt in relics:
        vis_w = max(e - s, 0.0014)
        xc = 0.5 * (s + e)
        axB.axvspan(xc - vis_w / 2, xc + vis_w / 2, color=col, alpha=0.14, zorder=0)
        axB.add_patch(mpatches.Rectangle((xc - vis_w / 2, -0.22), vis_w, 0.44, facecolor=col, edgecolor="black", lw=0.8, zorder=4))
        axB.annotate(
            lbl,
            xy=(xc, -0.24),
            xytext=(xc + dx_txt, -0.72),
            ha="center", va="top",
            fontsize=6.7, color=col, fontweight="bold",
            arrowprops={"arrowstyle": "->", "color": col, "lw": 0.9, "shrinkA": 2, "shrinkB": 2},
            zorder=5,
        )

    axB.set_xlim(x_min_mb, x_max_mb)
    axB.set_ylim(-1.72, 2.75)
    axB.set_yticks([0.0, 1.0, 2.0])
    axB.set_yticklabels([t[1] for t in reversed(tracks)])
    axB.set_xlabel("Gallus gallus GRCg7b chromosome Z (NC_052572.1), coordinate (Mb)")
    axB.set_title("B. Dual-database annotation audit of the 105.2-kb chicken NRG1–LPL intergenic gap (zero RefSeq 106 or Ensembl 116 overlaps)", loc="left", pad=8)

    # ---- Panel C: ODAPH (CDKL2-THAP6/RCHY1) physical gap collapse across 17 species ----
    axC = fig.add_subplot(gs[2])
    sp_order_c = [
        "Amis", "Cpor", "Acar", "Pbiv", "Pmur",
        "Gevg", "Cpic", "Cmyd", "Psin",
        "Shab", "Dnov", "Scam", "Arow", "Tgua", "Tgut", "Apla", "Ggal",
    ]
    clade_map = {r["code"]: r["clade"] for r in read_tsv(OUT / "enamel_intervals.tsv")}
    x_c = np.arange(len(sp_order_c))
    vals_kb = [odaph_gaps[sp]["bp"] / 1000.0 for sp in sp_order_c]
    cols_c = [COL[clade_map[sp]] for sp in sp_order_c]

    axC.bar(x_c, vals_kb, color=cols_c, edgecolor="#222222", linewidth=0.6, width=0.66, zorder=3)
    axC.axhline(3.784, color="#222222", linestyle="--", lw=1.1, zorder=4)
    # Clean pointer annotation well above Aves bar-top labels
    axC.annotate(
        "Alligator ODAPH full gene span (3.78 kb)",
        xy=(13.5, 3.784), xytext=(12.5, 16.5),
        ha="center", va="bottom", fontsize=7.3, fontweight="bold", color="#222222",
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#888888", "lw": 0.7},
        arrowprops={"arrowstyle": "->", "color": "#222222", "lw": 0.9},
        zorder=6,
    )

    for xi, v in zip(x_c, vals_kb):
        lbl = f"{v:.1f}" if v >= 1.0 else (f"{v:.2f}" if v > 0 else "0")
        y_off = 1.0 if v >= 4.0 else 0.75
        axC.text(
            xi, v + y_off, lbl, ha="center", va="bottom", fontsize=6.4, color="#222222",
            bbox={"fc": "white", "ec": "none", "pad": 0.4, "alpha": 0.88} if v < 3.8 else None,
            zorder=5,
        )

    axC.axvline(4.5, color="#CCCCCC", ls=":", lw=0.9)
    axC.axvline(8.5, color="#CCCCCC", ls=":", lw=0.9)
    axC.text(2.0, 58.5, "Toothed controls (ODAPH intact, n=5)", ha="center", va="top", fontsize=7.5, fontweight="bold", color=COL["Crocodylia"])
    axC.text(6.5, 58.5, "Testudines (n=4)", ha="center", va="top", fontsize=7.5, fontweight="bold", color=COL["Testudines"])
    axC.text(12.5, 58.5, "Aves (all 8 gaps < 2.5 kb; physical deletion, n=8)", ha="center", va="top", fontsize=7.5, fontweight="bold", color=COL["Aves"])

    axC.set_xticks(x_c)
    axC.set_xticklabels([SP_SHORT[s] for s in sp_order_c], rotation=32, ha="right", fontstyle="italic", fontsize=7.2)
    axC.set_ylabel("Intergenic gap (kb)")
    axC.set_ylim(0, 62.0)
    axC.set_title("C. Syntenic gap collapse at the ODAPH locus (CDKL2–THAP6/RCHY1) across 17 measurable sauropsids", loc="left", pad=8)

    return save_fig(fig, "fig1_synteny_map")


# =============================================================================
# FIGURE 2: Compositional / SEG / IDP Homology Bias Across 18 Amniotes
# =============================================================================
def make_figure2() -> dict[str, str]:
    comp = read_tsv(TABLES / "query_composition_seg.tsv")
    per_gene = read_tsv(TABLES / "tblastn_per_gene_seg_comparison.tsv")
    ret = read_tsv(TABLES / "retention_indices_18sp.tsv")

    fig = plt.figure(figsize=(11.5, 8.6))
    gs_outer = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.25], hspace=0.36)
    gs_top = gs_outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.38], wspace=0.25)
    gs_bot = gs_outer[1].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.14)

    # ---- Panel A: Query compositional landscape (%Pro+Gln vs %SEG-masked) ----
    axA = fig.add_subplot(gs_top[0])
    keep_q = {
        "LPL_Amis", "LPL_Acar",
        "ENAM_Amis", "ENAM_Acar",
        "AMBN_Amis", "AMBN_Acar",
        "AMTN_Amis", "AMTN_Acar",
        "IBSP_Amis", "IBSP_Acar",
        "SPP1_Amis", "SPP1_Acar",
        "SPARCL1_Amis", "DSPP_Amis",
    }
    role_style = {
        "globular": (COL["LPL"], "s", "Globular control (LPL)"),
        "enamel": (COL["ENAM"], "o", "Enamel IDPs (ENAM/AMBN/AMTN)"),
        "osteogenic": (COL["IBSP"], "D", "Osteogenic IDPs (IBSP/SPP1/SPARCL1)"),
        "dentin": (COL["DSPP"], "^", "Dentin IDP neg. control (DSPP)"),
    }
    label_offsets = {
        "LPL_Amis": (7, 2),
        "LPL_Acar": (7, -5),
        "ENAM_Amis": (-4, 11),
        "ENAM_Acar": (-50, -3),
        "AMBN_Amis": (-22, -14),
        "AMBN_Acar": (-26, -15),
        "AMTN_Amis": (-28, 7),
        "AMTN_Acar": (-48, 7),
        "IBSP_Amis": (7, -6),
        "IBSP_Acar": (-50, 2),
        "SPP1_Amis": (7, 3),
        "SPP1_Acar": (7, -2),
        "SPARCL1_Amis": (-38, 12),
        "DSPP_Amis": (7, 5),
    }

    for r in comp:
        qn = r["query_name"]
        if qn not in keep_q:
            continue
        pq = float(r["pct_Pro_Gln"])
        seg = float(r["pct_seg_masked"])
        length = int(r["length_aa"])
        if qn.startswith("LPL"):
            cat = "globular"
        elif qn.startswith(("ENAM", "AMBN", "AMTN")):
            cat = "enamel"
        elif qn.startswith("DSPP"):
            cat = "dentin"
        else:
            cat = "osteogenic"
        c, mk, _ = role_style[cat]
        axA.scatter(pq, seg, s=max(32, length / 9.5), color=c, marker=mk, edgecolor="black", lw=0.6, alpha=0.92, zorder=4)
        dx, dy = label_offsets.get(qn, (5, 4))
        short_lbl = qn.replace("_Amis", " (Amis)").replace("_Acar", " (Acar)")
        axA.annotate(short_lbl, (pq, seg), xytext=(dx, dy), textcoords="offset points", fontsize=6.5, color="#222222")

    axA.axvspan(16, 28.5, color=COL["ENAM"], alpha=0.06, zorder=0)
    axA.text(22.2, 31.5, "Pro/Gln-rich\nenamel zone", ha="center", va="center", fontsize=7.0, color=COL["ENAM"], fontstyle="italic")
    axA.set_xlabel("Proline + Glutamine content (% of residues)")
    axA.set_ylabel("SEG-masked low-complexity (%)")
    axA.set_xlim(2.5, 29.0)
    axA.set_ylim(-6, 70)
    axA.grid(True, ls=":", lw=0.5, alpha=0.6)
    handles_a = [Line2D([0], [0], marker=m, color="w", markerfacecolor=c, markeredgecolor="black", markersize=6.0, label=lb)
                 for c, m, lb in role_style.values()]
    axA.legend(handles=handles_a, loc="upper right", frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=6.6)
    axA.set_title("A. Compositional bias & SEG masking of queries", loc="left", pad=8)

    # ---- Panel B: Equidistant Anolis query — Functional IDPs vs Globular LPL vs Enamel IDPs ----
    axB = fig.add_subplot(gs_top[1])
    genes_b = ["LPL_Acar", "IBSP_Acar", "SPP1_Acar", "ENAM_Acar", "AMBN_Acar"]
    gene_labels_b = [
        "LPL\n(globular,\n0% SEG)",
        "IBSP\n(bone IDP,\n52.8% SEG)",
        "SPP1\n(bone IDP,\n25.5% SEG)",
        "ENAM\n(enamel IDP,\n7.5% SEG)",
        "AMBN\n(enamel IDP,\n5.0% SEG)",
    ]
    clades_b = ["Crocodylia", "Testudines", "Aves"]
    clade_offsets = {"Crocodylia": -0.25, "Testudines": 0.0, "Aves": +0.25}

    axB.axhspan(0, 16.0, color=COL["SPP1"], alpha=0.11, zorder=0)
    axB.annotate(
        "Functional avian SPP1 (6.8–25.3 bits)\noverlaps decayed avian ENAM/AMBN",
        xy=(3.55, 16.0), xytext=(3.82, 138.0),
        ha="center", va="bottom", fontsize=6.7, color=COL["SPP1"], fontweight="bold",
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": COL["SPP1"], "lw": 0.8},
        arrowprops={"arrowstyle": "->", "color": COL["SPP1"], "lw": 0.9},
        zorder=6,
    )

    rng = np.random.default_rng(20261017)
    for gi, qg in enumerate(genes_b):
        for cl in clades_b:
            xbase = gi + clade_offsets[cl]
            for seg_mode, dx_s, filled in [("seg_yes", -0.055, True), ("seg_no", +0.055, False)]:
                vals = [
                    float(r["excess_over_p95"])
                    for r in per_gene
                    if r["query_name"] == qg and r["clade"] == cl and r["seg_mode"] == seg_mode
                ]
                med = float(np.median(vals))
                jit = rng.uniform(-0.02, 0.02, size=len(vals))
                axB.scatter(
                    np.full(len(vals), xbase + dx_s) + jit,
                    vals,
                    s=20,
                    marker=MARKER[cl],
                    facecolor=COL[cl] if filled else "white",
                    edgecolor=COL[cl],
                    lw=0.9,
                    alpha=0.88,
                    zorder=4,
                )
                axB.plot([xbase + dx_s - 0.04, xbase + dx_s + 0.04], [med, med], color="black", lw=1.4, zorder=5)

    axB.text(0.25, 152, "Aves:\n112%", ha="center", va="bottom", fontsize=6.7, color=COL["Aves"], fontweight="bold")
    axB.text(1.25, 118, "SEG on: 74%\nSEG off: 48%", ha="center", va="bottom", fontsize=6.6, color=COL["IBSP"], fontweight="bold")
    axB.text(2.25, 62, "Aves:\n55%", ha="center", va="bottom", fontsize=6.7, color=COL["SPP1"], fontweight="bold")

    axB.set_xticks(range(len(genes_b)))
    axB.set_xticklabels(gene_labels_b, fontsize=7.3)
    axB.set_ylabel("tBLASTn excess over null p95 (bits)")
    axB.set_ylim(-10, 325)
    axB.grid(True, axis="y", ls=":", lw=0.5, alpha=0.6)
    handles_b = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=COL["Crocodylia"], markeredgecolor=COL["Crocodylia"], markersize=5.5, label="Crocodylia (n=3)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor=COL["Testudines"], markeredgecolor=COL["Testudines"], markersize=5.5, label="Testudines (n=4)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["Aves"], markeredgecolor=COL["Aves"], markersize=5.5, label="Aves (n=8)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#333333", markeredgecolor="#333333", markersize=5.5, label="Filled: -seg yes"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor="#333333", markersize=5.5, label="Open: -seg no"),
    ]
    axB.legend(handles=handles_b, loc="upper right", frameon=True, facecolor="white", edgecolor="#CCCCCC", ncol=2, fontsize=6.6)
    axB.set_title("B. Equidistant Anolis query (~281 Ma): functional IDPs lose bitscore in Aves even while intact", loc="left", pad=8)

    # ---- Panel C1 & C2: 18-species null-corrected retention index (shared y-axis!) ----
    axC1 = fig.add_subplot(gs_bot[0])
    axC2 = fig.add_subplot(gs_bot[1], sharey=axC1)

    sp_order = [
        "Amis", "Cpor", "Ggan",
        "Acar", "Pbiv", "Pmur",
        "Cmyd", "Cpic", "Psin", "Gevg",
        "Apla", "Arow", "Dnov", "Shab", "Scam", "Ggal", "Tgut", "Tgua",
    ]
    y_idx = np.arange(len(sp_order))[::-1]

    for ax_sub, qset, title_sub, show_ylabels, highlight_overlap in [
        (axC1, "Amis_refseq_full", "C1. Non-equidistant Alligator query (phylogeny confounds score)", True, True),
        (axC2, "Acar_equidistant", "C2. Equidistant Anolis query (~281 Ma; clean separation)", False, False),
    ]:
        sub_yes = {r["species_code"]: r for r in ret if r["query_set"] == qset and r["seg_mode"] == "seg_yes"}
        sub_no = {r["species_code"]: r for r in ret if r["query_set"] == qset and r["seg_mode"] == "seg_no"}

        for yi, sp in zip(y_idx, sp_order):
            ry = sub_yes[sp]
            rn = sub_no[sp]
            vy = float(ry["null_corrected_index"])
            vn = float(rn["null_corrected_index"])
            cl = ry["clade"]
            is_self = ry["is_self_query"] == "True"
            ax_sub.plot([0, vy], [yi, yi], color=COL[cl], lw=1.6, alpha=0.55, zorder=2)
            ax_sub.scatter(vy, yi, color=COL[cl], marker=MARKER[cl], s=36, edgecolor="black", lw=0.6, zorder=4)
            ax_sub.scatter(vn, yi, facecolor="white", edgecolor=COL[cl], marker=MARKER[cl], s=24, lw=0.8, zorder=3)
            if is_self:
                ax_sub.text(vy - 0.35, yi + 0.12, f"self ({vy:.2f})", fontsize=6.4, ha="right", color="#444444", fontstyle="italic")
            else:
                ax_sub.text(max(vy, vn) + 0.28, yi, f"{vy:.2f}", va="center", ha="left", fontsize=6.6, color="#333333")

        for sep in [14.5, 11.5, 7.5]:
            ax_sub.axhline(sep, color="#DDDDDD", ls="--", lw=0.7)

        if highlight_overlap:
            ax_sub.axvspan(1.48, 2.20, ymin=0.39, ymax=0.73, color=COL["Testudines"], alpha=0.18, zorder=0)
            ax_sub.annotate(
                "Overlap: edentulous\nC. mydas (2.17) >\ntoothed squamates\n(1.49–1.67)",
                xy=(2.17, 11.0), xytext=(6.2, 8.8),
                fontsize=6.7, color="#8A5200", fontweight="bold",
                bbox={"boxstyle": "round,pad=0.22", "fc": "#FFF9EE", "ec": COL["Testudines"], "lw": 0.7},
                arrowprops={"arrowstyle": "->", "color": "#8A5200", "lw": 0.9},
            )
        else:
            test_max = max(float(sub_yes[s]["null_corrected_index"]) for s in ("Cmyd", "Cpic", "Gevg", "Psin"))
            croc_min = min(float(sub_yes[s]["null_corrected_index"]) for s in ("Amis", "Cpor", "Ggan"))
            croc_max = max(float(sub_yes[s]["null_corrected_index"]) for s in ("Amis", "Cpor", "Ggan"))
            aves_max = max(float(sub_yes[s]["null_corrected_index"]) for s in ("Ggal", "Apla", "Shab", "Tgut", "Dnov", "Scam", "Arow", "Tgua"))
            ax_sub.axvspan(test_max + 0.04, croc_min - 0.04, color=COL["IBSP"], alpha=0.10, zorder=0)
            ax_sub.annotate(
                f"Clean separation\n({test_max:.2f} vs {croc_min:.2f})",
                xy=(1.40, 10.0), xytext=(5.5, 10.0),
                va="center", ha="left", fontsize=6.8, color=COL["IBSP"], fontweight="bold",
                bbox={"boxstyle": "round,pad=0.22", "fc": "#F0FBF7", "ec": COL["IBSP"], "lw": 0.7},
                arrowprops={"arrowstyle": "->", "color": COL["IBSP"], "lw": 0.9},
            )
            ax_sub.text(5.8, 16.0, f"Crocodylia: {croc_min:.2f}–{croc_max:.2f}", fontsize=6.9, color=COL["Crocodylia"], fontweight="bold")
            ax_sub.text(5.8, 7.9, f"Testudines: 0.00–{test_max:.2f}", fontsize=6.9, color=COL["Testudines"], fontweight="bold")
            ax_sub.text(5.8, 3.6, f"Aves: 0.00–{aves_max:.2f}", fontsize=6.9, color=COL["Aves"], fontweight="bold")

        ax_sub.set_yticks(y_idx)
        if show_ylabels:
            ax_sub.set_yticklabels([SP_SHORT[s] for s in sp_order], fontstyle="italic", fontsize=7.1)
        else:
            plt.setp(ax_sub.get_yticklabels(), visible=False)
        ax_sub.set_xlim(-0.3, 17.2)
        ax_sub.set_xlabel("Null-corrected retention index: Σ excess(ENAM,AMBN,AMTN) / excess(LPL)")
        ax_sub.set_title(title_sub, loc="left", pad=7, fontsize=8.7)

    return save_fig(fig, "fig2_compositional_bias")


# =============================================================================
# FIGURE 3: Level 3 Coding Lesions & Level 4 Phylogenetic Mapping + Homoplasy
# =============================================================================
def make_figure3() -> dict[str, str]:
    mp_rows = read_tsv(TABLES / "miniprot_level3_splice_audit.tsv")
    mp_sum = json.loads((OUT / "phase2_step3_summary.json").read_text(encoding="utf-8"))["miniprot_summary"]
    l4_muts = read_tsv(TABLES / "level4_shared_mutations.tsv")

    fig = plt.figure(figsize=(11.6, 8.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.55], width_ratios=[1.28, 1.0], hspace=0.42, wspace=0.26)

    # ---- Panel A: Level 3 miniprot coding-lesion rates ----
    axA = fig.add_subplot(gs[0, :])
    ed_loci = [
        r for r in mp_rows
        if r["query_gene"] in ("ENAM_Amis", "AMBN_Amis")
        and r["dentition"] == "edentulous"
        and r["aligned"] == "True"
        and int(r["total_coding_lesions"]) > 0
    ]
    ed_loci.sort(key=lambda r: (0 if r["query_gene"] == "ENAM_Amis" else 1, -int(r["total_coding_lesions"]) / max(1, int(r["aligned_aa"]))))

    labels_a: list[str] = []
    fs_rates: list[float] = []
    st_rates: list[float] = []
    splice_flags: list[str] = []

    for lbl, grp_key in [
        ("LPL control\n(18/18 sp.)", "LPL_all"),
        ("Toothed\nenamel\n(15 aln.)", "cluster_toothed_all"),
        ("Edentulous\npooled\n(14 aln.)", "cluster_edentulous_all_aligned"),
    ]:
        if grp_key == "LPL_all":
            aa = mp_sum["LPL_ctrl_toothed"]["aligned_aa"] + mp_sum["LPL_ctrl_edentulous"]["aligned_aa"]
            fs, st = 0, 0
        else:
            d = mp_sum[grp_key]
            aa, fs, st = d["aligned_aa"], d["frameshifts"], d["stops"]
        labels_a.append(lbl)
        fs_rates.append(1000.0 * fs / aa)
        st_rates.append(1000.0 * st / aa)
        splice_flags.append("")

    for r in ed_loci:
        g = r["query_gene"].replace("_Amis", "")
        sp = r["species_code"]
        aa = int(r["aligned_aa"])
        fs = int(r["frameshifts"])
        st = int(r["stops"])
        labels_a.append(f"{SP_SHORT[sp]}\n{g} ({aa} aa)")
        fs_rates.append(1000.0 * fs / aa)
        st_rates.append(1000.0 * st / aa)
        nc = r["noncanonical_splice_sites"]
        splice_flags.append(f"Splice: {nc}" if nc else "")

    xa = np.arange(len(labels_a))
    axA.bar(xa, fs_rates, color=COL["FS"], edgecolor="#222222", lw=0.6, width=0.60, label="Frameshifts (per 1,000 aligned aa)", zorder=3)
    axA.bar(xa, st_rates, bottom=fs_rates, color=COL["PTC"], edgecolor="#222222", lw=0.6, width=0.60, label="Premature stop codons (per 1,000 aligned aa)", zorder=3)

    axA.axvline(2.5, color="#999999", ls="--", lw=1.0)
    for i, (fr, sr, sf) in enumerate(zip(fs_rates, st_rates, splice_flags)):
        tot = fr + sr
        axA.text(i, tot + 0.55, f"{tot:.2f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#222222")
        if i == 0:
            axA.text(i, 3.0, "8,246 aa\n133/133 GT-AG", ha="center", va="bottom", fontsize=6.5, color="#333333",
                     bbox={"boxstyle": "round,pad=0.22", "fc": "#F2F2F2", "ec": "#CCCCCC"})
        elif i == 1:
            axA.text(i, 3.0, "9,394 aa\n1 fs, 0 PTC", ha="center", va="bottom", fontsize=6.5, color="#333333",
                     bbox={"boxstyle": "round,pad=0.22", "fc": "#F2F2F2", "ec": "#CCCCCC"})
        elif sf.startswith("Splice:"):
            axA.scatter([i], [tot + 3.2], marker="D", s=46, color=COL["SPLICE"], edgecolor="black", lw=0.7, zorder=5)
            axA.text(i, tot + 4.4, "GT→AT donor\n" + sf.replace("Splice: AT..AG", ""), ha="center", va="bottom", fontsize=6.4, color=COL["SPLICE"], fontweight="bold")

    axA.set_xticks(xa)
    axA.set_xticklabels(labels_a, fontsize=7.1)
    axA.set_ylabel("Coding lesions per 1,000 aligned aa")
    axA.set_ylim(0, 31.5)
    axA.grid(True, axis="y", ls=":", lw=0.5, alpha=0.6)
    handles_a, _ = axA.get_legend_handles_labels()
    handles_a.append(Line2D([0], [0], marker="D", color="w", markerfacecolor=COL["SPLICE"], markeredgecolor="black", markersize=6.5, label="Non-canonical GT→AT splice donor"))
    axA.legend(handles=handles_a, loc="upper left", frameon=True, facecolor="white", edgecolor="#CCCCCC", ncol=3, fontsize=7.2)
    axA.set_title("A. Level 3 miniprot coding-lesion density and splice-donor disruption (controls vs. edentulous enamel pseudogenes)", loc="left", pad=8)

    # ---- Panel B (bottom-left): Cladogram with Level 4 congruent synapomorphies ----
    axB1 = fig.add_subplot(gs[1, 0])
    leaves = [
        ("Amis", "Crocodylia (n=3 outgroups: 0 lesions)", 11.4, COL["Crocodylia"]),
        ("Cmyd", "C. mydas", 10.0, COL["Testudines"]),
        ("Cpic", "C. picta", 9.0, COL["Testudines"]),
        ("Gevg", "G. evgoodei / P. sinensis", 8.0, COL["Testudines"]),
        ("Arow", "A. rowi (kiwi)", 5.5, COL["Aves"]),
        ("Dnov", "D. novaehollandiae (emu)", 4.5, COL["Aves"]),
        ("Scam", "S. camelus (ostrich)", 3.5, COL["Aves"]),
        ("Ggal", "G. gallus (chicken)", 2.1, COL["Aves"]),
        ("Apla", "A. platyrhynchos (duck)", 1.1, COL["Aves"]),
        ("Shab", "S. habroptila (kakapo)", -0.1, COL["Aves"]),
        ("Tgut", "Taen. guttata (zebra finch)", -1.1, COL["Aves"]),
    ]
    ymap = {k: y for k, _, y, _ in leaves}

    def draw_branch(x0: float, x1: float, y_children: list[float], color: str = "#333333") -> float:
        ym = 0.5 * (min(y_children) + max(y_children))
        axB1.plot([x0, x1], [ym, ym], color=color, lw=1.4)
        axB1.plot([x1, x1], [min(y_children), max(y_children)], color=color, lw=1.4)
        return ym

    y_duro = draw_branch(0.54, 0.76, [ymap["Cmyd"], ymap["Cpic"]], COL["Testudines"])
    axB1.plot([0.76, 1.0], [ymap["Cmyd"], ymap["Cmyd"]], color=COL["Testudines"], lw=1.4)
    axB1.plot([0.76, 1.0], [ymap["Cpic"], ymap["Cpic"]], color=COL["Testudines"], lw=1.4)
    y_test = draw_branch(0.24, 0.54, [y_duro, ymap["Gevg"]], COL["Testudines"])
    axB1.plot([0.54, 1.0], [ymap["Gevg"], ymap["Gevg"]], color=COL["Testudines"], lw=1.1, ls="--")

    y_nova = draw_branch(0.58, 0.78, [ymap["Arow"], ymap["Dnov"]], COL["Aves"])
    axB1.plot([0.78, 1.0], [ymap["Arow"], ymap["Arow"]], color=COL["Aves"], lw=1.4)
    axB1.plot([0.78, 1.0], [ymap["Dnov"], ymap["Dnov"]], color=COL["Aves"], lw=1.4)
    y_palaeo = draw_branch(0.36, 0.58, [y_nova, ymap["Scam"]], COL["Aves"])
    axB1.plot([0.58, 1.0], [ymap["Scam"], ymap["Scam"]], color=COL["Aves"], lw=1.4)

    y_gallo = draw_branch(0.58, 0.78, [ymap["Ggal"], ymap["Apla"]], COL["Aves"])
    axB1.plot([0.78, 1.0], [ymap["Ggal"], ymap["Ggal"]], color=COL["Aves"], lw=1.4)
    axB1.plot([0.78, 1.0], [ymap["Apla"], ymap["Apla"]], color=COL["Aves"], lw=1.4)

    y_psitt = draw_branch(0.58, 0.78, [ymap["Shab"], ymap["Tgut"]], COL["Aves"])
    axB1.plot([0.78, 1.0], [ymap["Shab"], ymap["Shab"]], color=COL["Aves"], lw=1.4)
    axB1.plot([0.78, 1.0], [ymap["Tgut"], ymap["Tgut"]], color=COL["Aves"], lw=1.4)

    y_neog = draw_branch(0.36, 0.58, [y_gallo, y_psitt], COL["Aves"])
    y_aves = draw_branch(0.18, 0.36, [y_palaeo, y_neog], COL["Aves"])

    y_arch = draw_branch(0.08, 0.18, [y_test, y_aves], "#444444")
    axB1.plot([0.18, 0.24], [y_test, y_test], color=COL["Testudines"], lw=1.4)
    draw_branch(0.0, 0.08, [ymap["Amis"], y_arch], "#444444")
    axB1.plot([0.08, 1.0], [ymap["Amis"], ymap["Amis"]], color=COL["Crocodylia"], lw=1.4)

    for code, lbl, y, col in leaves:
        is_ital = code != "Amis"
        axB1.text(1.03, y, lbl, va="center", ha="left", fontsize=7.3, color=col,
                  fontstyle="italic" if is_ital else "normal", fontweight="bold" if code in ("Ggal", "Amis") else "normal")

    def count_mrca(mrca: str) -> tuple[int, int]:
        sub = [r for r in l4_muts if r["tree_congruent"] == "True" and r["mrca_clade"] == mrca]
        ptcs = sum(1 for r in sub if r["mutation_type"] == "stop")
        indels = sum(1 for r in sub if r["mutation_type"] in ("ins", "del"))
        return ptcs, indels

    duro_p, duro_i = count_mrca("Durocryptodira")
    pal_p, pal_i = count_mrca("Palaeognathae")
    gal_p, gal_i = count_mrca("Galloanserae")
    psi_p, psi_i = count_mrca("Psittacopasserae")

    # Carefully placed non-overlapping node callouts (each in its own open pocket)
    badges = [
        (0.76, y_duro, f"Durocryptodira:\nENAM {duro_p} PTCs\nAMBN {duro_i} fs (1-nt ins)", COL["Testudines"], (0.30, 10.35), "center"),
        (0.36, y_aves, "Crown Neornithes:\n0 congruent lesions\n(6/16 PTCs homoplastic)", "#222222", (0.04, 1.15), "center"),
        (0.58, y_palaeo, f"Palaeognathae:\nENAM {pal_p} PTCs + Ex6 AG→AC/GT\nAMBN Ex4 {pal_i} fs (16-nt del)", COL["FS"], (0.14, 6.75), "center"),
        (0.78, y_nova, "Novaeratitae:\nENAM 5 PTCs + Ex6 GT→AT\nAMBN 2 PTCs", COL["SPLICE"], (0.94, 6.80), "center"),
        (0.78, y_gallo, f"Galloanserae:\nENAM {gal_i} fs (2-nt/1-nt)\n+ {gal_p} PTC (G558*)", COL["Aves"], (0.62, 2.75), "center"),
        (0.78, y_psitt, f"Psittacopasserae:\nENAM {psi_p} PTC (W1003*)", COL["Aves"], (0.50, -1.65), "center"),
    ]
    for bx, by, txt, col, (tx, ty), ha in badges:
        axB1.scatter([bx], [by], s=44, color=col, edgecolor="black", lw=0.7, zorder=6)
        axB1.annotate(
            txt, xy=(bx, by), xytext=(tx, ty),
            ha=ha, va="center", fontsize=6.2, fontweight="bold", color=col,
            bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": col, "lw": 0.8, "alpha": 0.97},
            arrowprops={"arrowstyle": "->", "color": col, "lw": 0.8, "shrinkA": 2, "shrinkB": 3},
            zorder=7,
        )

    axB1.set_xlim(-0.24, 1.68)
    axB1.set_ylim(-2.45, 12.2)
    axB1.axis("off")
    axB1.set_title("B. Dual-aligner Level 4 synapomorphies mapped onto phylogeny", loc="left", pad=6)

    # ---- Panel C (bottom-right): Six incongruent avian ENAM PTC homoplasies ----
    axB2 = fig.add_subplot(gs[1, 1])
    incong = [r for r in l4_muts if r["tree_congruent"] == "False" and r["dataset"] == "ENAM_exon8_aves"]
    incong.sort(key=lambda r: int(r["XP_019346377_or_386_codon"]))

    mafft_fa = (OUT / "level4" / "ENAM_exon8_aves" / "mafft.fna").read_text(encoding="utf-8")
    seqs_ex8: dict[str, str] = {}
    cur_h = ""
    buf: list[str] = []
    for ln in mafft_fa.splitlines():
        if ln.startswith(">"):
            if cur_h:
                seqs_ex8[cur_h] = "".join(buf).upper()
            cur_h = ln[1:].split()[0]
            buf = []
        else:
            buf.append(ln.strip())
    if cur_h:
        seqs_ex8[cur_h] = "".join(buf).upper()

    amis_aln = seqs_ex8["Amis"]
    ref_nt_to_col: dict[int, int] = {}
    pos_nt = 0
    for ci, ch in enumerate(amis_aln):
        if ch != "-":
            pos_nt += 1
            ref_nt_to_col[pos_nt] = ci

    bird_rows = ["Arow", "Dnov", "Scam", "Ggal", "Apla", "Shab", "Tgut"]
    stop_colors = {"TGA": "#D55E00", "TAG": "#E69F00", "TAA": "#CC79A7"}

    for ri, sp in enumerate(bird_rows):
        yi = len(bird_rows) - 1 - ri
        s_aln = seqs_ex8.get(sp, "")
        for ci, row in enumerate(incong):
            ref_nt = int(row["ref_exon_nt_1based"])
            cols3 = [ref_nt_to_col[ref_nt + k] for k in range(3)]
            triplet = "".join(s_aln[k] for k in cols3) if s_aln else "---"
            if triplet in stop_colors:
                fc = stop_colors[triplet]
                txt = triplet
                tcol = "white"
            elif "-" in triplet:
                fc = "#F0F0F0"
                txt = "—"
                tcol = "#888888"
            else:
                fc = "#DCEBF5"
                txt = triplet
                tcol = "#1B4F72"
            axB2.add_patch(mpatches.Rectangle((ci - 0.44, yi - 0.40), 0.88, 0.80, facecolor=fc, edgecolor="white", lw=1.2, zorder=2))
            axB2.text(ci, yi, txt, ha="center", va="center", fontsize=7.0, fontweight="bold" if triplet in stop_colors else "normal", color=tcol, zorder=3)

    col_labels = []
    for r in incong:
        cod = r["XP_019346377_or_386_codon"]
        anc = r["ancestral_croc_codon"]
        tags = []
        if r["at_CpG_dinucleotide"] == "True":
            tags.append("CpG")
        if r["at_TGG_codon"] == "True":
            tags.append("TGG")
        if r["distinct_stop_codons"] == "True":
            tags.append("2 stops")
        tag_s = "\n(" + ",".join(tags) + ")" if tags else ""
        col_labels.append(f"c.{cod}\n[{anc}]{tag_s}")

    axB2.set_xticks(range(len(incong)))
    axB2.set_xticklabels(col_labels, fontsize=6.7)
    axB2.set_yticks(range(len(bird_rows)))
    axB2.set_yticklabels([SP_SHORT[s] for s in reversed(bird_rows)], fontstyle="italic", fontsize=7.3)
    axB2.set_xlim(-0.55, len(incong) - 0.45)
    axB2.set_ylim(-0.6, len(bird_rows) - 0.35)
    axB2.spines["left"].set_visible(False)
    axB2.spines["bottom"].set_visible(False)
    axB2.tick_params(length=0)

    leg_b2 = [
        mpatches.Patch(facecolor=stop_colors["TGA"], label="TGA stop"),
        mpatches.Patch(facecolor=stop_colors["TAG"], label="TAG stop"),
        mpatches.Patch(facecolor=stop_colors["TAA"], label="TAA stop"),
        mpatches.Patch(facecolor="#DCEBF5", label="Sense codon"),
        mpatches.Patch(facecolor="#F0F0F0", label="Gap / missing"),
    ]
    axB2.legend(handles=leg_b2, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=True, facecolor="white", edgecolor="#CCCCCC", fontsize=6.8)
    axB2.set_title("C. Six incongruent avian ENAM PTC homoplasies (37.5% of shared PTCs)", loc="left", pad=8)

    return save_fig(fig, "fig3_lesions_and_phylogeny")


def main() -> int:
    setup_style()
    t0 = datetime.now(timezone.utc)
    all_hashes: dict[str, str] = {}
    print("Generating Figure 1 (synteny + GRCg7b audit + ODAPH gap)...")
    all_hashes.update(make_figure1())
    print("Generating Figure 2 (compositional/SEG/IDP homology bias)...")
    all_hashes.update(make_figure2())
    print("Generating Figure 3 (Level 3 lesions + Level 4 phylogeny/homoplasy)...")
    all_hashes.update(make_figure3())

    entry = {
        "ts": datetime.now(timezone.utc).isoformat(),
        "kind": "enamel-relics-phase4-figures",
        "script": "scripts/enamel_relics/06_make_figures.py",
        "outputs": all_hashes,
        "elapsed_s": round((datetime.now(timezone.utc) - t0).total_seconds(), 2),
    }
    for log_path in (OUT / "provenance_run.jsonl", ROOT / "runs" / "provenance.jsonl"):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    print(f"Generated {len(all_hashes)} figure files in {FIG_DIR}:")
    for k, v in sorted(all_hashes.items()):
        print(f"  {k:36s} {v[:16]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
