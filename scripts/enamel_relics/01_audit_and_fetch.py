"""Phase 2.1 — Deterministic audit of 18 amniote RefSeq assemblies (+ GRCh38),
official NCBI RefSeq & Ensembl annotation checks, syntenic interval extraction
(NRG1-LPL enamel cluster + SPARCL1/SPP1/IBSP osteogenic IDP controls), and
query protein + composition-matched permutation FASTA generation.

Records SHA-256 checksums, UTC timestamps, API versions, and appends an entry
to runs/provenance.jsonl and runs/enamel_relics/provenance_run.jsonl.

Usage:
    uv run python scripts/enamel_relics/01_audit_and_fetch.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
import sys
import time
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "enamel_relics"
INTERVALS_ENAMEL = OUT / "intervals_enamel"
INTERVALS_IDP = OUT / "intervals_idp"
QUERIES_DIR = OUT / "queries"
REFDATA_DIR = OUT / "refdata"

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENS_API = "https://rest.ensembl.org"

SEED = 20261017
N_SHUFFLES = 20
FLANK_ENAMEL = 60_000
FLANK_SINGLE_ANCHOR = 800_000
FLANK_IDP = 100_000
MAX_SPAN = 3_000_000


def load_env_key(name: str) -> str | None:
    val = os.environ.get(name)
    if val:
        return val.strip()
    env_file = ROOT / ".env"
    if env_file.exists():
        for ln in env_file.read_text(encoding="utf-8").splitlines():
            ln = ln.strip()
            if ln.startswith(f"{name}="):
                return ln.split("=", 1)[1].strip().strip('"').strip("'")
    return None


NCBI_KEY = load_env_key("NCBI_API_KEY")
EMAIL = "limatheus458@gmail.com"


def eu_params(**extra: Any) -> dict[str, Any]:
    p: dict[str, Any] = {"tool": "bio-stack-enamel-relics", "email": EMAIL, **extra}
    if NCBI_KEY:
        p["api_key"] = NCBI_KEY
    return p


SESSION = requests.Session()
SESSION.headers["User-Agent"] = "bio-stack-enamel-relics/2.0 (mailto:limatheus458@gmail.com)"
if NCBI_KEY:
    SESSION.headers["api-key"] = NCBI_KEY


@dataclass(frozen=True)
class SpeciesSpec:
    code: str
    scientific_name: str
    common_name: str
    taxid: int
    assembly_acc: str
    clade: str
    subclade: str
    dentition: str  # "toothed" or "edentulous"
    ensembl_name: str | None = None


SPECIES: list[SpeciesSpec] = [
    SpeciesSpec("Amis", "Alligator mississippiensis", "American alligator", 8496, "GCF_030867095.1", "Crocodylia", "Crocodylia", "toothed"),
    SpeciesSpec("Ggan", "Gavialis gangeticus", "Gharial", 94835, "GCF_001723915.1", "Crocodylia", "Crocodylia", "toothed"),
    SpeciesSpec("Cpor", "Crocodylus porosus", "Saltwater crocodile", 8502, "GCF_001723895.1", "Crocodylia", "Crocodylia", "toothed"),
    SpeciesSpec("Acar", "Anolis carolinensis", "Green anole", 28377, "GCF_035594765.1", "Squamata", "Squamata", "toothed", "anolis_carolinensis"),
    SpeciesSpec("Pbiv", "Python bivittatus", "Burmese python", 176946, "GCF_059951505.1", "Squamata", "Squamata", "toothed"),
    SpeciesSpec("Pmur", "Podarcis muralis", "Common wall lizard", 64176, "GCF_964188315.1", "Squamata", "Squamata", "toothed", "podarcis_muralis"),
    SpeciesSpec("Cpic", "Chrysemys picta", "Painted turtle", 8478, "GCF_011386835.1", "Testudines", "Testudines", "edentulous", "chrysemys_picta_bellii"),
    SpeciesSpec("Cmyd", "Chelonia mydas", "Green sea turtle", 8469, "GCF_015237465.2", "Testudines", "Testudines", "edentulous"),
    SpeciesSpec("Gevg", "Gopherus evgoodei", "Goode's thornscrub tortoise", 1825980, "GCF_007399415.2", "Testudines", "Testudines", "edentulous", "gopherus_evgoodei"),
    SpeciesSpec("Psin", "Pelodiscus sinensis", "Chinese softshell turtle", 13735, "GCF_049634645.1", "Testudines", "Testudines", "edentulous", "pelodiscus_sinensis"),
    SpeciesSpec("Ggal", "Gallus gallus", "Chicken", 9031, "GCF_016699485.2", "Aves", "Galloanserae", "edentulous", "gallus_gallus"),
    SpeciesSpec("Apla", "Anas platyrhynchos", "Mallard duck", 8839, "GCF_047663525.1", "Aves", "Galloanserae", "edentulous", "anas_platyrhynchos"),
    SpeciesSpec("Shab", "Strigops habroptila", "Kakapo", 2489341, "GCF_004027225.2", "Aves", "Psittacopasserae", "edentulous", "strigops_habroptila"),
    SpeciesSpec("Tgut", "Taeniopygia guttata", "Zebra finch", 59729, "GCF_048771995.1", "Aves", "Psittacopasserae", "edentulous", "taeniopygia_guttata"),
    SpeciesSpec("Dnov", "Dromaius novaehollandiae", "Emu", 8790, "GCF_036370855.1", "Aves", "Palaeognathae", "edentulous", "dromaius_novaehollandiae"),
    SpeciesSpec("Scam", "Struthio camelus", "Common ostrich", 8801, "GCF_040807025.1", "Aves", "Palaeognathae", "edentulous", "struthio_camelus_australis"),
    SpeciesSpec("Arow", "Apteryx rowi", "Okarito brown kiwi", 308060, "GCF_003343035.1", "Aves", "Palaeognathae", "edentulous", "apteryx_rowi"),
    SpeciesSpec("Tgua", "Tinamus guttatus", "White-throated tinamou", 94827, "GCF_059628265.1", "Aves", "Palaeognathae", "edentulous"),
]

HUMAN_REF = SpeciesSpec("Hsap", "Homo sapiens", "Human (outgroup ref)", 9606, "GCF_000001405.40", "Mammalia", "Mammalia", "toothed", "homo_sapiens")

PANEL_TOOTH = ["ENAM", "AMBN", "AMTN", "AMELX", "MMP20", "DSPP", "ODAM", "ODAPH", "ACP4"]
PANEL_SYNTENY_CTRL = ["NRG1", "LPL", "PSD3", "IBSP", "SPP1", "SPARCL1", "DMP1", "MEPE", "CDKL2", "THAP6", "RCHY1"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def append_provenance(event: dict[str, Any]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    line = json.dumps(event, ensure_ascii=False) + "\n"
    with (OUT / "provenance_run.jsonl").open("a", encoding="utf-8") as f:
        f.write(line)
    with (ROOT / "runs" / "provenance.jsonl").open("a", encoding="utf-8") as f:
        f.write(line)


def get_json(url: str, params: dict[str, Any] | None = None, timeout: int = 60, retries: int = 3) -> dict[str, Any] | None:
    for attempt in range(retries):
        try:
            r = SESSION.get(url, params=params, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(1.5 * (attempt + 1))
                continue
            return None
        except requests.RequestException:
            time.sleep(1.5 * (attempt + 1))
    return None


def fetch_assembly_reports(accs: list[str]) -> dict[str, dict[str, Any]]:
    """Fetch official assembly & annotation metadata via NCBI Datasets CLI (WSL) or REST v2."""
    reports: dict[str, dict[str, Any]] = {}
    url = f"{DS_API}/genome/accession/{','.join(accs)}/dataset_report"
    data = get_json(url, timeout=90)
    if data and "reports" in data:
        for rep in data["reports"]:
            acc = rep.get("accession")
            if acc:
                reports[acc] = rep
        return reports
    # Fallback to WSL datasets CLI
    cmd = ["wsl", "-d", "Ubuntu-26.04", "--", "bash", "-lc",
           "datasets summary genome accession " + " ".join(accs)]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode == 0 and res.stdout.strip():
        data = json.loads(res.stdout)
        for rep in data.get("reports", []):
            acc = rep.get("accession")
            if acc:
                reports[acc] = rep
    return reports


def query_gene_batch_by_taxon(symbols: list[str], taxid: int, target_asm: str) -> dict[str, dict[str, Any] | None]:
    """Query NCBI Datasets v2 gene endpoint for a batch of symbols in one taxon,
    preferring genomic locations matching target_asm."""
    url = f"{DS_API}/gene/symbol/{','.join(symbols)}/taxon/{taxid}"
    data = get_json(url, params={"returned_content": "COMPLETE"}, timeout=90)
    out: dict[str, dict[str, Any] | None] = {s: None for s in symbols}
    if not data:
        return out
    sym_upper = {s.upper(): s for s in symbols}
    for rep in data.get("reports", []):
        g = rep.get("gene", {})
        gsym = g.get("symbol", "").upper()
        if gsym not in sym_upper:
            continue
        key = sym_upper[gsym]
        anns = g.get("annotations", [])
        chosen = None
        # First pass: exact assembly accession match
        for a in anns:
            if a.get("assembly_accession") == target_asm:
                locs = a.get("genomic_locations", [])
                if locs:
                    chosen = (a, locs[0])
                    break
        # Fallback: first annotation with genomic_locations
        if chosen is None:
            for a in anns:
                locs = a.get("genomic_locations", [])
                if locs:
                    chosen = (a, locs[0])
                    break
        if chosen is not None:
            a, gl = chosen
            gr = gl.get("genomic_range", {})
            orient = gr.get("orientation", "")
            strand = "+" if orient == "plus" else "-" if orient == "minus" else orient
            out[key] = {
                "gene_id": g.get("gene_id"),
                "symbol": g.get("symbol"),
                "description": g.get("description", ""),
                "type": g.get("type", ""),
                "assembly_accession": a.get("assembly_accession", ""),
                "assembly_name": a.get("assembly_name", ""),
                "annotation_release": a.get("release_name", ""),
                "acc": gl.get("genomic_accession_version", ""),
                "seq": gl.get("sequence_name", ""),
                "ini": int(gr.get("begin", 0)),
                "fim": int(gr.get("end", 0)),
                "strand": strand,
                "proteins": [p.get("accession_version") for p in g.get("proteins", []) if p.get("accession_version")],
            }
    return out


def fetch_fasta_slice(acc: str, ini: int, fim: int, label: str, dest: Path) -> int:
    """Download genomic slice via NCBI efetch if not already present."""
    if dest.exists() and dest.stat().st_size > 100:
        seq = "".join(ln.strip() for ln in dest.read_text(encoding="utf-8").splitlines() if not ln.startswith(">"))
        return len(seq)
    ini = max(1, ini)
    for attempt in range(4):
        try:
            r = SESSION.get(
                f"{EU}/efetch.fcgi",
                params=eu_params(db="nuccore", id=acc, rettype="fasta", retmode="text",
                                 seq_start=ini, seq_stop=fim),
                timeout=300,
            )
            if r.status_code == 200 and r.text.strip().startswith(">"):
                body = r.text.strip().split("\n", 1)[1]
                dest.write_text(f">{label} {acc}:{ini}-{fim}\n{body}\n", encoding="utf-8", newline="\n")
                return len(body.replace("\n", "").replace("\r", ""))
            time.sleep(2.0 * (attempt + 1))
        except requests.RequestException:
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"Failed efetch for {label} ({acc}:{ini}-{fim})")


def fetch_protein_seq(acc: str) -> tuple[str, str]:
    """Return (header_line, amino_acid_sequence) from NCBI protein efetch."""
    for attempt in range(4):
        try:
            r = SESSION.get(
                f"{EU}/efetch.fcgi",
                params=eu_params(db="protein", id=acc, rettype="fasta", retmode="text"),
                timeout=120,
            )
            if r.status_code == 200 and r.text.strip().startswith(">"):
                lines = r.text.strip().splitlines()
                return lines[0][1:], "".join(ln.strip() for ln in lines[1:])
            time.sleep(1.5 * (attempt + 1))
        except requests.RequestException:
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Failed protein efetch for {acc}")


def find_refseq_protein_for_gene(gene_id: str | int) -> str | None:
    """Resolve RefSeq XP_/NP_ protein accession from NCBI Gene via Datasets v2."""
    url = f"{DS_API}/gene/id/{gene_id}/product_report"
    data = get_json(url, timeout=60)
    if not data:
        return None
    best: tuple[int, str] | None = None
    for rep in data.get("reports", []):
        prod = rep.get("product", {})
        for tr in prod.get("transcripts", []):
            prot = tr.get("protein", {})
            acc = prot.get("accession_version", "")
            plen = int(prot.get("length", 0) or 0)
            if acc.startswith(("NP_", "XP_")):
                # Prefer NP_ over XP_, then longest isoform
                score = (1_000_000 if acc.startswith("NP_") else 0) + plen
                if best is None or score > best[0]:
                    best = (score, acc)
    return best[1] if best else None


def check_ensembl_chicken() -> dict[str, Any]:
    """Query Ensembl REST API for current release and genes in GRCg7b chrZ:55109784-55214993."""
    result: dict[str, Any] = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint": ENS_API,
        "status": "unknown",
        "ensembl_releases": None,
        "chicken_assembly": None,
        "genes_in_NRG1_LPL_gap": [],
        "symbol_lookup": {},
    }
    try:
        r = requests.get(f"{ENS_API}/info/data/?content-type=application/json", timeout=25)
        if r.status_code == 200:
            result["ensembl_releases"] = r.json().get("releases")
        r_asm = requests.get(f"{ENS_API}/info/assembly/gallus_gallus?content-type=application/json", timeout=25)
        if r_asm.status_code == 200:
            d = r_asm.json()
            result["chicken_assembly"] = {
                "assembly_name": d.get("assembly_name"),
                "assembly_accession": d.get("assembly_accession"),
                "genebuild_last_geneset_update": d.get("genebuild_last_geneset_update"),
            }
        r_ov = requests.get(
            f"{ENS_API}/overlap/region/gallus_gallus/Z:55109784-55214993?feature=gene;content-type=application/json",
            timeout=30,
        )
        if r_ov.status_code == 200:
            genes = r_ov.json()
            result["genes_in_NRG1_LPL_gap"] = [
                {
                    "id": g.get("id"),
                    "external_name": g.get("external_name"),
                    "biotype": g.get("biotype"),
                    "start": g.get("start"),
                    "end": g.get("end"),
                    "strand": g.get("strand"),
                }
                for g in genes
            ]
        for sym in ["ENAM", "AMBN", "AMTN", "AMELX", "NRG1", "LPL"]:
            rs = requests.get(
                f"{ENS_API}/xrefs/symbol/gallus_gallus/{sym}?content-type=application/json",
                timeout=20,
            )
            result["symbol_lookup"][sym] = rs.json() if rs.status_code == 200 else f"HTTP {rs.status_code}"
        result["status"] = "ok"
    except requests.RequestException as exc:
        result["status"] = f"network_error: {exc}"
    return result


def download_refdata_wsl() -> dict[str, str]:
    """Download Alligator proteome + GFF3 (GCF_030867095.1) and Chicken GFF3 (GCF_016699485.2) via WSL datasets CLI."""
    REFDATA_DIR.mkdir(parents=True, exist_ok=True)
    hashes: dict[str, str] = {}

    amis_prot = REFDATA_DIR / "Amis_protein.faa"
    amis_gff = REFDATA_DIR / "Amis_genomic.gff"
    ggal_gff = REFDATA_DIR / "Ggal_genomic.gff"

    wsl_ref = "/mnt/c/Users/user/bio-stack/runs/enamel_relics/refdata"

    if not (amis_prot.exists() and amis_gff.exists()):
        print("[refdata] Downloading Alligator GCF_030867095.1 (protein + gff3)...")
        cmd = [
            "wsl", "-d", "Ubuntu-26.04", "--", "bash", "-lc",
            f"mkdir -p {wsl_ref} && datasets download genome accession GCF_030867095.1 "
            f"--include protein,gff3 --filename {wsl_ref}/amis_ref.zip",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode != 0:
            raise RuntimeError(f"datasets Amis failed: {res.stderr}")
        with zipfile.ZipFile(REFDATA_DIR / "amis_ref.zip") as zf:
            for member in zf.namelist():
                if member.endswith("protein.faa"):
                    amis_prot.write_bytes(zf.read(member))
                elif member.endswith("genomic.gff"):
                    amis_gff.write_bytes(zf.read(member))
        (REFDATA_DIR / "amis_ref.zip").unlink(missing_ok=True)

    if not ggal_gff.exists():
        print("[refdata] Downloading Gallus GCF_016699485.2 (gff3)...")
        cmd = [
            "wsl", "-d", "Ubuntu-26.04", "--", "bash", "-lc",
            f"datasets download genome accession GCF_016699485.2 "
            f"--include gff3 --filename {wsl_ref}/ggal_ref.zip",
        ]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode != 0:
            raise RuntimeError(f"datasets Ggal failed: {res.stderr}")
        with zipfile.ZipFile(REFDATA_DIR / "ggal_ref.zip") as zf:
            for member in zf.namelist():
                if member.endswith("genomic.gff"):
                    ggal_gff.write_bytes(zf.read(member))
        (REFDATA_DIR / "ggal_ref.zip").unlink(missing_ok=True)

    for p in (amis_prot, amis_gff, ggal_gff):
        hashes[p.name] = sha256_file(p)
    return hashes


def audit_gff_overlaps(ggal_gff: Path, amis_gff: Path) -> dict[str, Any]:
    """Programmatically verify RefSeq GFF3 features in GRCg7b chrZ NRG1-LPL gap
    and exact Alligator ENAM/AMBN/AMTN/ODAPH exon coordinates."""
    ggal_gap_features: list[dict[str, Any]] = []
    ggal_release = ""
    with ggal_gff.open(encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("#!annotation-source") or ln.startswith("#!genome-build"):
                ggal_release += ln.strip() + " | "
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != "NC_052572.1":
                continue
            ftype, s, e, strand, attrs = parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]
            if ftype in ("gene", "pseudogene", "lnc_RNA", "mRNA", "region") and ftype != "region":
                # Gap between NRG1 (ends 55109784) and LPL (starts 55214993)
                if e >= 55_109_000 and s <= 55_235_000 and ftype in ("gene", "pseudogene"):
                    ggal_gap_features.append({
                        "type": ftype,
                        "start": s,
                        "end": e,
                        "strand": strand,
                        "attrs": attrs[:180],
                    })

    # Check overlap specifically with psi-ENAM (55155480-55157241) and psi-AMBN (55168855-55168992)
    enam_overlaps = [x for x in ggal_gap_features if x["start"] <= 55_157_241 and x["end"] >= 55_155_480]
    ambn_overlaps = [x for x in ggal_gap_features if x["start"] <= 55_168_992 and x["end"] >= 55_168_855]

    # Parse Alligator ENAM (XP_019346377.1), AMBN (XP_019346386.1), AMTN (XP_014451634.1), ODAPH (LOC106740521 / XP_019349100.2)
    amis_genes: dict[str, dict[str, Any]] = {}
    amis_cds: dict[str, list[tuple[str, int, int, str]]] = {
        "XP_019346377.1": [],
        "XP_019346386.1": [],
        "XP_014451634.1": [],
        "XP_019349100.2": [],
        "XP_019346320.1": [],
    }
    amis_release = ""
    with amis_gff.open(encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("#!annotation-source") or ln.startswith("#!genome-build"):
                amis_release += ln.strip() + " | "
            if ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, ftype, s, e, strand, attrs = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6], parts[8]
            if ftype == "gene":
                for target in ("Name=ENAM;", "Name=AMBN;", "Name=AMTN;", "Name=LPL;", "Name=NRG1;", "Name=PSD3;", "Name=LOC106740521;"):
                    if target in attrs:
                        gname = target.split("=")[1].rstrip(";")
                        amis_genes[gname] = {"seqid": seqid, "start": s, "end": e, "strand": strand}
            elif ftype == "CDS":
                for pid in amis_cds:
                    if f"protein_id={pid}" in attrs or f"Name={pid}" in attrs:
                        amis_cds[pid].append((seqid, s, e, strand))

    return {
        "ggal_gff_header": ggal_release.strip(" | "),
        "ggal_chrZ_55.109_55.235Mb_gene_features": ggal_gap_features,
        "ggal_overlaps_psiENAM": enam_overlaps,
        "ggal_overlaps_psiAMBN": ambn_overlaps,
        "amis_gff_header": amis_release.strip(" | "),
        "amis_genes": amis_genes,
        "amis_cds_exons": {k: sorted(v, key=lambda x: x[1]) for k, v in amis_cds.items()},
    }


def main() -> int:
    t0 = datetime.now(timezone.utc)
    for d in (OUT, INTERVALS_ENAMEL, INTERVALS_IDP, QUERIES_DIR, REFDATA_DIR):
        d.mkdir(parents=True, exist_ok=True)

    print("=== [1/6] Fetching official NCBI Datasets assembly reports (18 species + GRCh38) ===")
    all_specs = SPECIES + [HUMAN_REF]
    asm_reports = fetch_assembly_reports([s.assembly_acc for s in all_specs])
    print(f"  Retrieved {len(asm_reports)}/{len(all_specs)} assembly reports")

    print("=== [2/6] Auditing NCBI RefSeq annotations across all 18 species + GRCh38 ===")
    all_symbols = PANEL_TOOTH + PANEL_SYNTENY_CTRL
    annotations: dict[str, dict[str, Any]] = {}
    for sp in all_specs:
        # Query in two batches to keep URL clean
        res1 = query_gene_batch_by_taxon(PANEL_TOOTH, sp.taxid, sp.assembly_acc)
        time.sleep(0.15)
        res2 = query_gene_batch_by_taxon(PANEL_SYNTENY_CTRL, sp.taxid, sp.assembly_acc)
        time.sleep(0.15)
        annotations[sp.code] = {**res1, **res2}
        n_tooth = sum(1 for g in PANEL_TOOTH if res1[g] is not None)
        print(f"  {sp.code:4s} ({sp.scientific_name:28s}): tooth_panel={n_tooth}/9 | "
              f"NRG1={'Y' if res2['NRG1'] else '-'} LPL={'Y' if res2['LPL'] else '-'} "
              f"IBSP={'Y' if res2['IBSP'] else '-'} SPP1={'Y' if res2['SPP1'] else '-'} "
              f"SPARCL1={'Y' if res2['SPARCL1'] else '-'}")

    print("=== [3/6] Checking Ensembl REST API for GRCg7b chrZ NRG1-LPL annotation ===")
    ensembl_audit = check_ensembl_chicken()
    print(f"  Ensembl status: {ensembl_audit['status']} | releases={ensembl_audit['ensembl_releases']} | "
          f"genes in gap={len(ensembl_audit['genes_in_NRG1_LPL_gap'])}")

    print("=== [4/6] Downloading Alligator proteome/GFF3 & Chicken GFF3 + verifying GFF features ===")
    ref_hashes = download_refdata_wsl()
    gff_audit = audit_gff_overlaps(REFDATA_DIR / "Ggal_genomic.gff", REFDATA_DIR / "Amis_genomic.gff")
    print(f"  Ggal GFF overlaps with psi-ENAM: {len(gff_audit['ggal_overlaps_psiENAM'])}, "
          f"psi-AMBN: {len(gff_audit['ggal_overlaps_psiAMBN'])}")
    print(f"  Amis ENAM (XP_019346377.1) CDS exons: {len(gff_audit['amis_cds_exons']['XP_019346377.1'])}, "
          f"AMBN CDS exons: {len(gff_audit['amis_cds_exons']['XP_019346386.1'])}, "
          f"ODAPH CDS exons: {len(gff_audit['amis_cds_exons']['XP_019349100.2'])}")
    print(f"  Amis ODAPH gene coords: {gff_audit['amis_genes'].get('LOC106740521')}")

    print("=== [5/6] Downloading syntenic genomic intervals (18 ENAMEL NRG1-LPL + 18 SIBLING IDP) ===")
    enamel_intervals_meta: list[dict[str, Any]] = []
    idp_intervals_meta: list[dict[str, Any]] = []

    for sp in SPECIES:
        ann = annotations[sp.code]
        # --- Enamel cluster interval (NRG1..LPL) ---
        nrg1, lpl, psd3 = ann["NRG1"], ann["LPL"], ann["PSD3"]
        if nrg1 and lpl and nrg1["acc"] == lpl["acc"]:
            lo = min(nrg1["ini"], nrg1["fim"], lpl["ini"], lpl["fim"]) - FLANK_ENAMEL
            hi = max(nrg1["ini"], nrg1["fim"], lpl["ini"], lpl["fim"]) + FLANK_ENAMEL
            # True non-negative intergenic gap regardless of orientation
            if max(nrg1["ini"], nrg1["fim"]) < min(lpl["ini"], lpl["fim"]):
                gap_bp = min(lpl["ini"], lpl["fim"]) - max(nrg1["ini"], nrg1["fim"])
                order = "NRG1->LPL"
            else:
                gap_bp = min(nrg1["ini"], nrg1["fim"]) - max(lpl["ini"], lpl["fim"])
                order = "LPL->NRG1"
            acc, seq_name, anchor_mode = lpl["acc"], lpl["seq"], "pair_NRG1_LPL"
        elif nrg1 or lpl:
            p = lpl or nrg1
            assert p is not None
            mid = (p["ini"] + p["fim"]) // 2
            lo, hi = max(1, mid - FLANK_SINGLE_ANCHOR), mid + FLANK_SINGLE_ANCHOR
            gap_bp, order = None, f"single_{'LPL' if lpl else 'NRG1'}"
            acc, seq_name, anchor_mode = p["acc"], p["seq"], order
        else:
            raise RuntimeError(f"No NRG1 or LPL anchor for {sp.code}")

        dest_e = INTERVALS_ENAMEL / f"interval_{sp.code}.fna"
        bp_e = fetch_fasta_slice(acc, max(1, lo), hi, f"{sp.code}_NRG1_LPL", dest_e)
        enamel_intervals_meta.append({
            "code": sp.code,
            "species": sp.scientific_name,
            "clade": sp.clade,
            "dentition": sp.dentition,
            "assembly": sp.assembly_acc,
            "seq_acc": acc,
            "seq_name": seq_name,
            "anchor_mode": anchor_mode,
            "order": order,
            "nrg1_coords": f"{nrg1['ini']}-{nrg1['fim']}({nrg1['strand']})" if nrg1 else "NA",
            "lpl_coords": f"{lpl['ini']}-{lpl['fim']}({lpl['strand']})" if lpl else "NA",
            "psd3_coords": f"{psd3['ini']}-{psd3['fim']}({psd3['strand']})" if psd3 else "NA",
            "intergenic_gap_bp": gap_bp if gap_bp is not None else "NA",
            "win_start": max(1, lo),
            "win_end": hi,
            "win_bp": bp_e,
            "fasta": str(dest_e.relative_to(ROOT)),
            "sha256": sha256_file(dest_e),
        })
        print(f"  [ENAMEL] {sp.code:4s} {acc} (chr/scaf={seq_name}) win={max(1,lo)}-{hi} ({bp_e:,} bp, gap={gap_bp})")

        # --- Osteogenic SIBLING/IDP cluster interval (SPARCL1..SPP1..MEPE..IBSP..DMP1) ---
        idp_genes = [ann[g] for g in ("SPARCL1", "SPP1", "MEPE", "IBSP", "DMP1") if ann[g] is not None]
        # Group by chromosome/scaffold accession (usually all 3-5 sit within ~300-600 kb on 1 chromosome!)
        by_acc: dict[str, list[dict[str, Any]]] = {}
        for rec in idp_genes:
            by_acc.setdefault(rec["acc"], []).append(rec)
        # Fetch one combined FASTA containing all needed windows for this species' IDP controls
        idp_dest = INTERVALS_IDP / f"idp_{sp.code}.fna"
        idp_slices_desc: list[str] = []
        if not (idp_dest.exists() and idp_dest.stat().st_size > 100):
            chunks: list[str] = []
            for idx_a, (iacc, recs) in enumerate(by_acc.items()):
                # Cluster nearby genes within 1.5 Mb; otherwise slice per gene
                recs_sorted = sorted(recs, key=lambda r: r["ini"])
                clusters: list[list[dict[str, Any]]] = []
                for r in recs_sorted:
                    if not clusters or r["ini"] - clusters[-1][-1]["fim"] > 1_500_000:
                        clusters.append([r])
                    else:
                        clusters[-1].append(r)
                for idx_c, cl in enumerate(clusters):
                    ilo = max(1, min(r["ini"] for r in cl) - FLANK_IDP)
                    ihi = max(r["fim"] for r in cl) + FLANK_IDP
                    tmp_p = INTERVALS_IDP / f"_tmp_{sp.code}_{idx_a}_{idx_c}.fna"
                    syms = "_".join(r["symbol"] for r in cl)
                    fetch_fasta_slice(iacc, ilo, ihi, f"{sp.code}_IDP_{syms}", tmp_p)
                    chunks.append(tmp_p.read_text(encoding="utf-8"))
                    tmp_p.unlink(missing_ok=True)
                    idp_slices_desc.append(f"{iacc}:{ilo}-{ihi}({syms})")
            idp_dest.write_text("".join(chunks), encoding="utf-8", newline="\n")
        else:
            for ln in idp_dest.read_text(encoding="utf-8").splitlines():
                if ln.startswith(">"):
                    idp_slices_desc.append(ln[1:].strip())
        idp_bp = sum(len(ln.strip()) for ln in idp_dest.read_text(encoding="utf-8").splitlines() if not ln.startswith(">"))
        idp_intervals_meta.append({
            "code": sp.code,
            "species": sp.scientific_name,
            "clade": sp.clade,
            "dentition": sp.dentition,
            "genes_present": ",".join(sorted({r["symbol"] for r in idp_genes})),
            "slices": ";".join(idp_slices_desc),
            "total_bp": idp_bp,
            "fasta": str(idp_dest.relative_to(ROOT)),
            "sha256": sha256_file(idp_dest),
        })
        print(f"  [IDP   ] {sp.code:4s} genes={','.join(sorted({r['symbol'] for r in idp_genes}))} ({idp_bp:,} bp)")

    print("=== [6/6] Fetching query proteins & generating composition-matched shuffles (seed=20261017) ===")
    # Resolve Anolis SPP1 and SPARCL1 RefSeq accessions dynamically from verified annotations
    acar_spp1_acc = find_refseq_protein_for_gene(annotations["Acar"]["SPP1"]["gene_id"]) if annotations["Acar"]["SPP1"] else None
    acar_sparcl1_acc = find_refseq_protein_for_gene(annotations["Acar"]["SPARCL1"]["gene_id"]) if annotations["Acar"]["SPARCL1"] else None
    amis_ibsp_acc = find_refseq_protein_for_gene(annotations["Amis"]["IBSP"]["gene_id"]) if annotations["Amis"]["IBSP"] else None
    amis_spp1_acc = find_refseq_protein_for_gene(annotations["Amis"]["SPP1"]["gene_id"]) if annotations["Amis"]["SPP1"] else None
    amis_sparcl1_acc = find_refseq_protein_for_gene(annotations["Amis"]["SPARCL1"]["gene_id"]) if annotations["Amis"]["SPARCL1"] else None

    query_catalog: dict[str, tuple[str, str, str]] = {
        # name: (accession, source_species, role)
        "ENAM_Amis": ("XP_019346377.1", "Amis", "enamel_target_full_refseq"),
        "ENAM_Amis_legacy": ("QOX06214.1", "Amis", "enamel_target_legacy_partial"),
        "AMBN_Amis": ("XP_019346386.1", "Amis", "enamel_target"),
        "AMTN_Amis": ("XP_014451634.1", "Amis", "enamel_target_refseq"),
        "AMTN_Amis_gb": ("KYO35092.1", "Amis", "enamel_target_genbank"),
        "LPL_Amis": ("XP_019346320.1", "Amis", "globular_positive_control"),
        "DSPP_Amis": ("XP_059578269.1", "Amis", "cross_chr_negative_control"),
        "ENAM_Acar": ("XP_062825318.1", "Acar", "enamel_target_equidistant"),
        "AMBN_Acar": ("NP_001315156.1", "Acar", "enamel_target_equidistant"),
        "AMTN_Acar": ("NP_001315157.1", "Acar", "enamel_target_equidistant"),
        "LPL_Acar": ("XP_062827593.1", "Acar", "globular_positive_control_equidistant"),
        "IBSP_Acar": ("XP_003221343.2", "Acar", "idp_osteogenic_control_equidistant"),
    }
    if acar_spp1_acc:
        query_catalog["SPP1_Acar"] = (acar_spp1_acc, "Acar", "idp_osteogenic_control_equidistant")
    if acar_sparcl1_acc:
        query_catalog["SPARCL1_Acar"] = (acar_sparcl1_acc, "Acar", "idp_scpp_progenitor_control_equidistant")
    if amis_ibsp_acc:
        query_catalog["IBSP_Amis"] = (amis_ibsp_acc, "Amis", "idp_osteogenic_control")
    if amis_spp1_acc:
        query_catalog["SPP1_Amis"] = (amis_spp1_acc, "Amis", "idp_osteogenic_control")
    if amis_sparcl1_acc:
        query_catalog["SPARCL1_Amis"] = (amis_sparcl1_acc, "Amis", "idp_scpp_progenitor_control")

    fetched_queries: dict[str, dict[str, Any]] = {}
    enamel_real_lines: list[str] = []
    enamel_with_shuf_lines: list[str] = []
    idp_with_shuf_lines: list[str] = []

    rng = np.random.default_rng(SEED)
    for qname, (acc, src, role) in query_catalog.items():
        hdr, seq = fetch_protein_seq(acc)
        fetched_queries[qname] = {
            "query_name": qname,
            "accession": acc,
            "source_species": src,
            "role": role,
            "length_aa": len(seq),
            "header": hdr,
        }
        print(f"  [QUERY ] {qname:18s} {acc:16s} {len(seq):4d} aa  ({role})")
        real_block = f">{qname} {acc}\n{seq}\n"
        shuf_blocks = []
        for i in range(N_SHUFFLES):
            perm = "".join(rng.permutation(list(seq)))
            shuf_blocks.append(f">{qname}~shuf{i:02d}\n{perm}\n")

        if any(qname.startswith(p) for p in ("ENAM_", "AMBN_", "AMTN_", "LPL_", "DSPP_")):
            enamel_real_lines.append(real_block)
            enamel_with_shuf_lines.append(real_block)
            enamel_with_shuf_lines.extend(shuf_blocks)
        else:
            idp_with_shuf_lines.append(real_block)
            idp_with_shuf_lines.extend(shuf_blocks)

    (QUERIES_DIR / "enamel_real.faa").write_text("".join(enamel_real_lines), encoding="utf-8", newline="\n")
    (QUERIES_DIR / "enamel_with_shuffles.faa").write_text("".join(enamel_with_shuf_lines), encoding="utf-8", newline="\n")
    (QUERIES_DIR / "idp_with_shuffles.faa").write_text("".join(idp_with_shuf_lines), encoding="utf-8", newline="\n")

    # Individual miniprot query FASTAs (Alligator full-length RefSeq + legacy + Anolis)
    for qname in ("ENAM_Amis", "ENAM_Amis_legacy", "AMBN_Amis", "AMTN_Amis", "AMTN_Amis_gb", "LPL_Amis"):
        acc = fetched_queries[qname]["accession"]
        _, seq = fetch_protein_seq(acc)
        (QUERIES_DIR / f"mp_{qname}.faa").write_text(f">{qname} {acc}\n{seq}\n", encoding="utf-8", newline="\n")

    # Write structured audit outputs
    audit_payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "n_shuffles": N_SHUFFLES,
        "assemblies": {},
        "annotations": annotations,
        "ensembl_chicken_audit": ensembl_audit,
        "gff_audit": gff_audit,
        "enamel_intervals": enamel_intervals_meta,
        "idp_intervals": idp_intervals_meta,
        "queries": fetched_queries,
        "refdata_hashes": ref_hashes,
    }
    for sp in all_specs:
        rep = asm_reports.get(sp.assembly_acc, {})
        ai = rep.get("assembly_info", {})
        ast = rep.get("assembly_stats", {})
        ann_info = rep.get("annotation_info", {})
        audit_payload["assemblies"][sp.code] = {
            "code": sp.code,
            "scientific_name": sp.scientific_name,
            "common_name": sp.common_name,
            "taxid": sp.taxid,
            "clade": sp.clade,
            "subclade": sp.subclade,
            "dentition": sp.dentition,
            "assembly_accession": sp.assembly_acc,
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
        }

    audit_json = OUT / "audit_metadata.json"
    audit_json.write_text(json.dumps(audit_payload, indent=2, ensure_ascii=False), encoding="utf-8")

    # Export clean TSVs immediately for assemblies & intervals
    with (OUT / "assemblies_audit.tsv").open("w", encoding="utf-8", newline="") as f:
        rows = list(audit_payload["assemblies"].values())
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    with (OUT / "enamel_intervals.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(enamel_intervals_meta[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(enamel_intervals_meta)

    with (OUT / "idp_intervals.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(idp_intervals_meta[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(idp_intervals_meta)

    append_provenance({
        "ts": datetime.now(timezone.utc).isoformat(),
        "kind": "enamel-relics-phase2-step1",
        "script": "scripts/enamel_relics/01_audit_and_fetch.py",
        "seed": SEED,
        "n_shuffles": N_SHUFFLES,
        "n_species": len(SPECIES),
        "outputs": {
            "audit_metadata.json": sha256_file(audit_json),
            "enamel_with_shuffles.faa": sha256_file(QUERIES_DIR / "enamel_with_shuffles.faa"),
            "idp_with_shuffles.faa": sha256_file(QUERIES_DIR / "idp_with_shuffles.faa"),
        },
        "elapsed_s": round((datetime.now(timezone.utc) - t0).total_seconds(), 2),
    })
    print(f"\nDone in {(datetime.now(timezone.utc) - t0).total_seconds():.1f}s -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
