"""Phase 2.1b — Programmatic audit of Ensembl Release 116 Gallus gallus GRCg7b GTF
at the chrZ NRG1-LPL syntenic gap, and reconciliation of Pelodiscus sinensis
current RefSeq representative assembly (GCF_052426305.1, ASM5242630v1).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "enamel_relics"
REF = OUT / "refdata"
GTF_URL = "https://ftp.ensembl.org/pub/release-116/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz"
TOOTH_SYMS = ("ENAM", "AMBN", "AMTN", "AMELX", "MMP20", "DSPP", "ODAM", "ODAPH", "ACP4")


def main() -> int:
    REF.mkdir(parents=True, exist_ok=True)
    dest = REF / "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz"
    if not dest.exists():
        print(f"Downloading {GTF_URL}...")
        r = requests.get(GTF_URL, timeout=180)
        r.raise_for_status()
        dest.write_bytes(r.content)

    sha = hashlib.sha256(dest.read_bytes()).hexdigest()
    feats = []
    tooth_hits = []
    with gzip.open(dest, "rt", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if p[2] != "gene":
                continue
            attrs = p[8]
            for sym in TOOTH_SYMS:
                if f'gene_name "{sym}"' in attrs:
                    tooth_hits.append({"symbol": sym, "chr": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6], "attrs": attrs})
            if p[0] == "Z" and int(p[4]) >= 55_109_000 and int(p[3]) <= 55_235_000:
                feats.append({"chr": p[0], "start": int(p[3]), "end": int(p[4]), "strand": p[6], "attrs": attrs})

    ov_enam = [x for x in feats if x["start"] <= 55_157_241 and x["end"] >= 55_155_480]
    ov_ambn = [x for x in feats if x["start"] <= 55_168_992 and x["end"] >= 55_168_855]

    ens_out = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "ensembl_release": 116,
        "gtf_url": GTF_URL,
        "gtf_file": dest.name,
        "sha256": sha,
        "tooth_panel_symbol_hits_genome_wide": tooth_hits,
        "chrZ_55.109_55.235Mb_genes": feats,
        "overlaps_psiENAM_55155480_55157241": ov_enam,
        "overlaps_psiAMBN_55168855_55168992": ov_ambn,
    }
    (OUT / "ensembl_116_ggal_audit.json").write_text(json.dumps(ens_out, indent=2), encoding="utf-8")
    print(f"Ensembl 116 GTF SHA-256: {sha}")
    print(f"Genome-wide tooth-panel symbol hits: {len(tooth_hits)}")
    print(f"Overlaps with psi-ENAM: {len(ov_enam)} | psi-AMBN: {len(ov_ambn)}")
    print(f"Features in chrZ:55,109,000-55,235,000 ({len(feats)} genes):")
    for x in feats:
        print(f"  Z:{x['start']}-{x['end']} ({x['strand']}) {x['attrs'][:150]}")

    # Reconcile Psin assembly to current RefSeq representative GCF_052426305.1
    cmd = ["wsl", "-d", "Ubuntu-26.04", "--", "bash", "-lc", "datasets summary genome accession GCF_052426305.1"]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rep = json.loads(res.stdout)["reports"][0]
    ai, ast, ann_info = rep["assembly_info"], rep["assembly_stats"], rep["annotation_info"]

    meta_p = OUT / "audit_metadata.json"
    meta = json.loads(meta_p.read_text(encoding="utf-8"))
    meta["ensembl_116_ggal_audit"] = ens_out
    meta["assemblies"]["Psin"].update({
        "assembly_accession": "GCF_052426305.1",
        "assembly_name": ai.get("assembly_name", ""),
        "assembly_level": ai.get("assembly_level", ""),
        "assembly_status": ai.get("assembly_status", ""),
        "release_date": ai.get("release_date", ""),
        "refseq_category": ai.get("refseq_category", ""),
        "total_sequence_length": ast.get("total_sequence_length", ""),
        "contig_n50": ast.get("contig_n50", ""),
        "scaffold_n50": ast.get("scaffold_n50", ""),
        "annotation_name": ann_info.get("name", ""),
        "annotation_provider": ann_info.get("provider", ""),
        "annotation_release_date": ann_info.get("release_date", ""),
        "annotation_status": ann_info.get("status", ""),
    })
    for row in meta["enamel_intervals"]:
        if row["code"] == "Psin":
            row["assembly"] = "GCF_052426305.1"
    meta_p.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

    # Update TSVs cleanly
    with (OUT / "assemblies_audit.tsv").open("w", encoding="utf-8", newline="") as f:
        rows = list(meta["assemblies"].values())
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with (OUT / "enamel_intervals.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(meta["enamel_intervals"][0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(meta["enamel_intervals"])

    print("Updated Psin assembly to GCF_052426305.1 in audit_metadata.json and TSVs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
