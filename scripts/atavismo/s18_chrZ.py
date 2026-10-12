"""Ancoras sintenicas ARCOSSAURIAS (nao mamiferas): o cluster de esmalte do jacare e
flanqueado por LPL/PSD3. Em Gallus, LPL esta no cromossomo Z."""

from __future__ import annotations

import json
import time
from pathlib import Path

import requests

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
OUT = Path(r"C:\Users\Manec\bio-stack\runs\atavismo")
SEQ = OUT / "seqs"
CM = {"tool": "bio-stack", "email": "limatheus458@gmail.com"}
S = requests.Session()


def loc(sym, tax=9031, asm="GRCg7b"):
    r = S.get(
        f"{API}/gene/symbol/{sym}/taxon/{tax}", params={"returned_content": "COMPLETE"}, timeout=90
    ).json()
    for rep in r.get("reports", []):
        g = rep["gene"]
        if g["symbol"].upper() != sym.upper():
            continue
        for a in g.get("annotations", []):
            if asm in a.get("assembly_name", ""):
                gl = a["genomic_locations"][0]
                gr = gl["genomic_range"]
                return (
                    gl["genomic_accession_version"],
                    gl.get("sequence_name"),
                    int(gr["begin"]),
                    int(gr["end"]),
                    gr["orientation"],
                )
    return None


print("=== ancoras em Gallus (GRCg7b) ===")
for sym in [
    "LPL",
    "PSD3",
    "ABCA1",
    "SLC44A1",
    "FKTN",
    "ARHGAP6",
    "AMELY",
    "SPP1",
    "IBSP",
    "MEPE",
    "DMP1",
    "MMP20",
    "MMP27",
    "MMP13",
]:
    L = loc(sym)
    print(
        f"  {sym:9} -> {'NAO ANOTADO' if not L else f'{L[0]} (chr {L[1]}) {L[2]:,}-{L[3]:,} {L[4]}'}"
    )
    time.sleep(0.25)
REG = {}
L = loc("LPL")
REG["GGZ_lpl"] = (
    L[0],
    L[2] - 500_000,
    L[3] + 500_000,
    f"Gallus chr{L[1]} ao redor de LPL (ancora arcossauria do cluster de esmalte)",
)
L2 = loc("ARHGAP6")
if L2:
    REG["GG1_arhgap6"] = (
        L2[0],
        L2[2] - 400_000,
        L2[3] + 400_000,
        f"Gallus chr{L2[1]} ao redor de ARHGAP6 (ancora de AMEL)",
    )
L3 = loc("SPP1") or loc("IBSP") or loc("MEPE")
if L3:
    REG["GG4_sibling"] = (
        L3[0],
        L3[2] - 500_000,
        L3[3] + 500_000,
        f"Gallus chr{L3[1]} cluster SIBLING (ancora de DSPP)",
    )
for lab, (acc, s, e, desc) in REG.items():
    s = max(s, 1)
    r = S.get(
        f"{EU}/efetch.fcgi",
        params={
            **CM,
            "db": "nuccore",
            "id": acc,
            "seq_start": s,
            "seq_stop": e,
            "rettype": "fasta",
            "retmode": "text",
        },
        timeout=900,
    )
    r.raise_for_status()
    body = "".join(r.text.splitlines()[1:])
    (SEQ / f"reg_{lab}.fa").write_text(
        f">{lab} {acc}:{s}-{e} {desc}\n"
        + "\n".join(body[i : i + 80] for i in range(0, len(body), 80))
        + "\n",
        encoding="utf-8",
    )
    print(f"\n{lab:14} {acc}:{s:,}-{e:,} {len(body):,} bp N={body.upper().count('N')}  [{desc}]")
    REG[lab] = (acc, s, e, desc)
    time.sleep(0.5)
(OUT / "s18_regioes_chrZ.json").write_text(
    json.dumps(REG, indent=2, ensure_ascii=False), encoding="utf-8"
)
