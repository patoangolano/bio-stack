"""O mesmo locus ancorado em LPL em Taeniopygia e Melopsittacus."""

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
TAX = {
    "Tgut": (59729, "Taeniopygia guttata"),
    "Mund": (13146, "Melopsittacus undulatus"),
    "Mgal": (9103, "Meleagris gallopavo"),
}
info = {}
for k, (tax, nm) in TAX.items():
    r = S.get(
        f"{API}/gene/symbol/LPL/taxon/{tax}", params={"returned_content": "COMPLETE"}, timeout=90
    ).json()
    g = [x["gene"] for x in r.get("reports", []) if x["gene"]["symbol"] == "LPL"]
    if not g:
        print(k, "sem LPL")
        continue
    a = g[0]["annotations"][0]["genomic_locations"][0]
    gr = a["genomic_range"]
    acc = a["genomic_accession_version"]
    b, e = int(gr["begin"]), int(gr["end"])
    lo, hi = max(b - 450_000, 1), e + 120_000
    info[k] = {
        "sp": nm,
        "acc": acc,
        "chr": a.get("sequence_name"),
        "lpl": (b, e),
        "win": (lo, hi),
        "orient": gr["orientation"],
    }
    fr = S.get(
        f"{EU}/efetch.fcgi",
        params={
            **CM,
            "db": "nuccore",
            "id": acc,
            "seq_start": lo,
            "seq_stop": hi,
            "rettype": "fasta",
            "retmode": "text",
        },
        timeout=900,
    )
    fr.raise_for_status()
    body = "".join(fr.text.splitlines()[1:])
    (SEQ / f"reg_{k}_lpl.fa").write_text(
        f">{k}_lpl {acc}:{lo}-{hi} {nm}\n"
        + "\n".join(body[i : i + 80] for i in range(0, len(body), 80))
        + "\n",
        encoding="utf-8",
    )
    print(
        f"{k} {nm:26} LPL em {acc} (chr {a.get('sequence_name')}) {b:,}-{e:,} {gr['orientation']} | janela {lo:,}-{hi:,} ({len(body):,} bp, N={body.upper().count('N')})"
    )
    time.sleep(0.5)
(OUT / "s20_aves.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
