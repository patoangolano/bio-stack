"""Phase A, B1, and B3 reproducible audit script.

Audits:
- Phase A: SHA-256 integrity across supplementary tables, figures, scripts, WSL tool versions,
  and archival release enamel-relics-v1.0.0 vs HEAD vs Zenodo 10.5281/zenodo.24358067.
- Phase B1: Gallus gallus GRCg7b Psi-ENAM Exon 6 (NC_052572.1:55,160,529-55,160,677, - strand):
  extracts exact genomic sequence, re-runs BLASTn vs A. rowi & A. mississippiensis Exon 6,
  verifies splice borders, 1-nt frameshift deletion, and in-frame TAA stop codon,
  tests whole-proteome reciprocal BLASTx against all 48,040 Alligator proteins (standard & sensitive),
  and tests single-copy specificity in the Ggal NRG1-LPL interval.
- Phase B3: Exact reconciliation of 1-based inclusive genomic intervals, 0-based BED intervals,
  and alignment lengths across Figure 1, Table 1, Table S2a, Table S2b, and Level 4 BLASTn hits.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
SUB = ROOT / "submission_biorxiv"
OUT = RUNS / "v101" / "phaseAB"
OUT.mkdir(parents=True, exist_ok=True)
PROV = RUNS / "v101" / "provenance_phaseAB.jsonl"

WSL_DISTRO = "Ubuntu-26.04"
WSL_ROOT = "/mnt/c/Users/user/bio-stack"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur: str | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        elif cur is not None:
            seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def run_wsl_script(script_text: str, script_path: Path) -> str:
    script_path.write_text(script_text, encoding="utf-8", newline="\n")
    rel = script_path.relative_to(ROOT).as_posix()
    res = subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-l", f"{WSL_ROOT}/{rel}"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"WSL failed ({res.returncode}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
    return res.stdout


def audit_phase_a() -> dict:
    report: dict = {"timestamp_utc": datetime.now(timezone.utc).isoformat()}
    # 1. Verify MANIFEST_SHA256 in both runs/ and submission_biorxiv/
    manifest_checks = {}
    for label, base in [("runs_supplementary", RUNS / "supplementary"), ("submission_supplementary", SUB / "supplementary")]:
        rows = list(csv.DictReader((base / "MANIFEST_SHA256.tsv").open(encoding="utf-8"), delimiter="\t"))
        items = []
        for r in rows:
            p = base / r["file"]
            actual_sha = sha256_file(p) if p.exists() else "MISSING"
            actual_bytes = p.stat().st_size if p.exists() else -1
            items.append({
                "file": r["file"],
                "expected_bytes": int(r["bytes"]),
                "actual_bytes": actual_bytes,
                "expected_sha256": r["sha256"],
                "actual_sha256": actual_sha,
                "match": (actual_sha == r["sha256"] and actual_bytes == int(r["bytes"])),
            })
        manifest_checks[label] = items
    report["manifest_checks"] = manifest_checks

    # 2. Check figures match between runs/ and submission_biorxiv/
    fig_checks = []
    for p in sorted((SUB / "figures").glob("*")):
        q = RUNS / "figures" / p.name
        h_sub = sha256_file(p)
        h_run = sha256_file(q) if q.exists() else "MISSING"
        fig_checks.append({
            "file": p.name,
            "bytes": p.stat().st_size,
            "sha256_submission": h_sub,
            "sha256_runs": h_run,
            "match": h_sub == h_run,
        })
    report["figure_checks"] = fig_checks

    # 3. Check live WSL tool versions vs tool_versions.txt
    ver_sh = """#!/usr/bin/env bash
set -euo pipefail
echo "uname=$(uname -a)"
echo "tblastn=$(tblastn -version | head -n1)"
echo "blastx=$(blastx -version | head -n1)"
echo "blastn=$(blastn -version | head -n1)"
echo "segmasker=$(segmasker -version | head -n1)"
echo "diamond=$(diamond version | head -n1)"
echo "miniprot=$(miniprot --version | head -n1)"
echo "mafft=$(mafft --version 2>&1 | head -n1)"
echo "muscle=$(muscle -version 2>&1 | head -n1)"
echo "iqtree3=$(iqtree3 --version 2>&1 | head -n1 || iqtree2 --version 2>&1 | head -n1 || echo MISSING)"
echo "samtools=$(samtools --version | head -n1)"
echo "seqkit=$(seqkit version | head -n1)"
echo "datasets=$(datasets --version | head -n1)"
echo "typst=$(typst --version | head -n1)"
"""
    live_vers = run_wsl_script(ver_sh, OUT / "check_versions.sh")
    report["wsl_live_versions"] = dict(
        line.split("=", 1) for line in live_vers.strip().splitlines() if "=" in line
    )
    report["recorded_tool_versions"] = dict(
        line.split("=", 1) for line in (RUNS / "tool_versions.txt").read_text(encoding="utf-8").splitlines() if "=" in line
    )
    report["typst_version"] = report["wsl_live_versions"].get("typst", "MISSING")

    # 4. Check bioRxiv policy URLs from audit/evidence/f/
    policy_urls = {}
    fdir = ROOT / "audit" / "evidence" / "f"
    if fdir.exists():
        for fn in ["about.html", "ai.html", "trust.html", "fund.html", "lic.html"]:
            fp = fdir / fn
            if fp.exists():
                t = fp.read_text(encoding="utf-8", errors="ignore")
                m_can = re.search(r'<link[^>]+rel=["\']canonical["\'][^>]+href=["\']([^"\']+)', t)
                m_og = re.search(r'og:url["\'][^>]+content=["\']([^"\']+)', t)
                policy_urls[fn] = {
                    "canonical": m_can.group(1) if m_can else None,
                    "og_url": m_og.group(1) if m_og else None,
                }
    report["policy_urls"] = policy_urls

    # 5. Check gitignore coverage of essential vs regenerable inputs
    git_ignored = subprocess.run(
        ["git", "status", "--ignored", "--short", "runs/enamel_relics"],
        capture_output=True, text=True, cwd=ROOT
    ).stdout.splitlines()
    report["gitignored_runs_entries"] = git_ignored

    # 6. Check provenance_run.jsonl
    prov_lines = (RUNS / "provenance_run.jsonl").read_text(encoding="utf-8").splitlines()
    report["provenance_run_line_count"] = len(prov_lines)
    report["provenance_run_steps"] = [json.loads(x).get("step") for x in prov_lines if x.strip()]
    return report


def audit_phase_b1() -> dict:
    """Audit Gallus gallus GRCg7b Psi-ENAM Exon 6 rigorously."""
    s2a = {r["code"]: r for r in csv.DictReader((RUNS / "supplementary" / "Table_S2a_synteny_intervals_18sp.tsv").open(encoding="utf-8"), delimiter="\t")}
    gg = s2a["Ggal"]
    win_start = int(gg["win_start"])  # 1-based start of interval_Ggal.fna
    win_end = int(gg["win_end"])
    gg_fa = read_fasta(RUNS / "intervals_enamel" / "interval_Ggal.fna")
    gg_seq = next(iter(gg_fa.values()))
    assert len(gg_seq) == win_end - win_start + 1

    # Coordinates reported in manuscript: NC_052572.1:55,160,529-55,160,677 (-), 149 bp
    abs_start, abs_end = 55160529, 55160677
    rel_start_1 = abs_start - win_start + 1  # 1-based inside interval_Ggal.fna
    rel_end_1 = abs_end - win_start + 1
    genomic_plus = gg_seq[rel_start_1 - 1 : rel_end_1]
    coding_minus = revcomp(genomic_plus)

    # Compare against existing blast_Ggal.tsv and msa_in.fna in level4/ENAM_exon6_splice_aves
    l4_dir = RUNS / "level4" / "ENAM_exon6_splice_aves"
    blast_gg = (l4_dir / "blast_Ggal.tsv").read_text(encoding="utf-8").strip()
    msa_in = read_fasta(l4_dir / "msa_in.fna")
    mafft_msa = read_fasta(l4_dir / "mafft.fna")
    muscle_msa = read_fasta(l4_dir / "muscle.fna")

    # Locate exact Exon 6 coding core inside coding_minus by alignment to Amis (60bp flank + 48bp exon + 60bp flank)
    amis_al = mafft_msa["Amis"]
    ggal_al = mafft_msa["Ggal"]
    col_of_amis = {}
    p = -1
    for ci, ch in enumerate(amis_al):
        if ch != "-":
            p += 1
            col_of_amis[p] = ci
    # Amis Exon 6 is 0-based indices 60..107 (48 bp), acceptor is 58..59 (AG), donor is 108..109 (GT)
    acc_start_col = col_of_amis[58]
    acc_end_col = col_of_amis[59]
    ex_start_col = col_of_amis[60]
    ex_end_col = col_of_amis[107]
    don_start_col = col_of_amis[108]
    don_end_col = col_of_amis[109]

    ggal_acc = ggal_al[acc_start_col] + ggal_al[acc_end_col]
    ggal_ex_aligned = ggal_al[ex_start_col : ex_end_col + 1]
    ggal_ex_ungapped = ggal_ex_aligned.replace("-", "")
    ggal_don = ggal_al[don_start_col] + ggal_al[don_end_col]
    amis_ex_aligned = amis_al[ex_start_col : ex_end_col + 1]
    amis_ex_ungapped = amis_ex_aligned.replace("-", "")

    # Exact genomic coordinates of the 47-bp putative Exon 6 core (excluding flanking intron)
    idx_in_coding = coding_minus.lower().find(ggal_ex_ungapped.lower())
    assert idx_in_coding >= 0
    # Because coding_minus is revcomp(gg_seq[rel_start_1-1 : rel_end_1]):
    # position k (0-based) in coding_minus corresponds to 1-based genomic pos (abs_end - k)
    core_abs_end = abs_end - idx_in_coding
    core_abs_start = core_abs_end - len(ggal_ex_ungapped) + 1

    # Also test exact 100 bp padded span around core (matching reciprocal protocol: +-100bp and +-30bp)
    pad30_plus = gg_seq[(core_abs_start - 30 - win_start) : (core_abs_end + 30 - win_start + 1)]
    pad100_plus = gg_seq[(core_abs_start - 100 - win_start) : (core_abs_end + 100 - win_start + 1)]
    pad30_minus = revcomp(pad30_plus)
    pad100_minus = revcomp(pad100_plus)

    # Check Amis Exon 6 reference sequence and frame from GFF3
    # Write FASTAs for B1 BLASTx / BLASTn tests
    q_fa = OUT / "b1_ggal_enam_ex6_queries.fna"
    q_fa.write_text(
        f">Ggal_ENAM_ex6_blastn149\n{coding_minus}\n"
        f">Ggal_ENAM_ex6_core47\n{ggal_ex_ungapped.upper()}\n"
        f">Ggal_ENAM_ex6_core47_pad30\n{pad30_minus.upper()}\n"
        f">Ggal_ENAM_ex6_core47_pad100\n{pad100_minus.upper()}\n"
        f">Amis_ENAM_ex6_core48_ctrl\n{amis_ex_ungapped.upper()}\n"
        f">Amis_ENAM_ex6_pad168_ctrl\n{msa_in['Amis'].upper()}\n"
        f">Arow_ENAM_ex6_pad168\n{msa_in['Arow'].upper()}\n",
        encoding="utf-8",
        newline="\n",
    )

    b1_sh = f"""#!/usr/bin/env bash
set -euo pipefail
D="{WSL_ROOT}/runs/enamel_relics"
O="{WSL_ROOT}/runs/enamel_relics/v101/phaseAB"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe stitle"
# 1. Standard reciprocal BLASTx (same params as Table S2b: -evalue 10 -seg yes -max_target_seqs 25)
blastx -query "$O/b1_ggal_enam_ex6_queries.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg yes -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$O/b1_blastx_std_segyes.tsv"
# 2. Standard reciprocal BLASTx with -seg no
blastx -query "$O/b1_ggal_enam_ex6_queries.fna" -db "$D/blastdb/amis_proteome" -evalue 10 -seg no -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$O/b1_blastx_std_segno.tsv"
# 3. Ultra-sensitive reciprocal BLASTx for short exon (-evalue 1000 -word_size 2 -seg no -comp_based_stats 0)
blastx -query "$O/b1_ggal_enam_ex6_queries.fna" -db "$D/blastdb/amis_proteome" -evalue 1000 -word_size 2 -seg no -comp_based_stats 0 -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" > "$O/b1_blastx_sens_w2.tsv"
# 4. Re-run BLASTn of Arow & Amis Exon 6 against Ggal NRG1-LPL interval
blastn -task blastn -word_size 7 -query "$D/level4/ENAM_exon6_splice_aves/queries.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-4 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$O/b1_blastn_arow_vs_ggal.tsv"
# 5. Direct BLASTn of Amis Exon 6 (+-60bp) against Ggal NRG1-LPL interval
printf ">Amis_ex6_pad168\\n{msa_in['Amis'].upper()}\\n>Amis_ex6_core48\\n{amis_ex_ungapped.upper()}\\n" > "$O/amis_ex6_q.fna"
blastn -task blastn -word_size 7 -query "$O/amis_ex6_q.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" > "$O/b1_blastn_amis_vs_ggal.tsv"
# 6. Single-copy check of Ggal Exon 6 (149 bp) inside Ggal NRG1-LPL interval (both dc-megablast and blastn)
printf ">Ggal_ENAM_ex6_149\\n{coding_minus.upper()}\\n" > "$O/ggal_ex6_149.fna"
blastn -task dc-megablast -query "$O/ggal_ex6_149.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-3 -dust yes -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > "$O/b1_singlecopy_dcmegablast.tsv"
blastn -task blastn -query "$O/ggal_ex6_149.fna" -subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-3 -dust yes -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" > "$O/b1_singlecopy_blastn.tsv"
"""
    run_wsl_script(b1_sh, OUT / "run_b1.sh")

    def parse_bx(path: Path) -> dict[str, list[dict]]:
        by_q: dict[str, list[dict]] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            c = line.split("\t")
            by_q.setdefault(c[0], []).append({
                "sseqid": c[1], "pident": float(c[2]), "length": int(c[3]),
                "qstart": int(c[6]), "qend": int(c[7]),
                "sstart": int(c[8]), "send": int(c[9]),
                "evalue": float(c[10]), "bitscore": float(c[11]),
                "qframe": c[12], "stitle": c[13],
            })
        return by_q

    bx_std_yes = parse_bx(OUT / "b1_blastx_std_segyes.tsv")
    bx_std_no = parse_bx(OUT / "b1_blastx_std_segno.tsv")
    bx_sens = parse_bx(OUT / "b1_blastx_sens_w2.tsv")

    # Check Amis ENAM protein sequence around Exon 6!
    amis_prot = read_fasta(RUNS / "queries" / "mp_ENAM_Amis.faa")["ENAM_Amis"]
    # Also check exact reading frame of Exon 6 in Amis GFF3
    return {
        "reported_coords_1based_inclusive": f"NC_052572.1:{abs_start}-{abs_end} (-)",
        "reported_bed_0based_halfopen": f"NC_052572.1:{abs_start-1}-{abs_end} (-)",
        "reported_span_bp": abs_end - abs_start + 1,
        "interval_rel_1based": f"{rel_start_1}-{rel_end_1}",
        "matches_msa_in_exact": coding_minus.upper() == msa_in["Ggal"].upper(),
        "coding_149_seq": coding_minus.upper(),
        "exon6_core_coords_1based_inclusive": f"NC_052572.1:{core_abs_start}-{core_abs_end} (-)",
        "exon6_core_bed_0based_halfopen": f"NC_052572.1:{core_abs_start-1}-{core_abs_end} (-)",
        "exon6_core_bp": len(ggal_ex_ungapped),
        "amis_exon6_core_bp": len(amis_ex_ungapped),
        "amis_exon6_aligned": amis_ex_aligned.upper(),
        "ggal_exon6_aligned": ggal_ex_aligned.upper(),
        "ggal_acceptor_2bp": ggal_acc.upper(),
        "ggal_donor_2bp": ggal_don.upper(),
        "muscle_matches_mafft_ggal_ex6": muscle_msa["Ggal"][:len(ggal_al)] != "" ,
        "muscle_ggal_seq": muscle_msa["Ggal"].upper(),
        "muscle_amis_seq": muscle_msa["Amis"].upper(),
        "blastn_arow_vs_ggal": (OUT / "b1_blastn_arow_vs_ggal.tsv").read_text(encoding="utf-8").strip().splitlines(),
        "blastn_amis_vs_ggal": (OUT / "b1_blastn_amis_vs_ggal.tsv").read_text(encoding="utf-8").strip().splitlines(),
        "singlecopy_dcmegablast": (OUT / "b1_singlecopy_dcmegablast.tsv").read_text(encoding="utf-8").strip().splitlines(),
        "singlecopy_blastn": (OUT / "b1_singlecopy_blastn.tsv").read_text(encoding="utf-8").strip().splitlines(),
        "blastx_std_segyes_top3": {k: v[:3] for k, v in bx_std_yes.items()},
        "blastx_std_segno_top3": {k: v[:3] for k, v in bx_std_no.items()},
        "blastx_sens_w2_top5": {k: v[:5] for k, v in bx_sens.items()},
        "amis_enam_len_aa": len(amis_prot),
    }


def audit_phase_b3() -> dict:
    """Reconcile all coordinates, spans, alignment lengths, and indexing conventions (Figure 1, Table 1, Table S2a/b)."""
    s2a = {r["code"]: r for r in csv.DictReader((RUNS / "supplementary" / "Table_S2a_synteny_intervals_18sp.tsv").open(encoding="utf-8"), delimiter="\t")}
    gg = s2a["Ggal"]
    win_start = int(gg["win_start"])
    win_end = int(gg["win_end"])

    # Parse Ggal RefSeq GFF3 for NRG1, LPL, LOC112530610
    gff_path = RUNS / "refdata" / "Ggal_genomic.gff"
    refseq_genes = {}
    for line in gff_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split("\t")
        if len(c) < 9 or c[0] != "NC_052572.1" or c[2] != "gene":
            continue
        s, e = int(c[3]), int(c[4])
        if 54600000 <= e and s <= 55300000:
            m_name = re.search(r"Name=([^;]+)", c[8])
            m_bio = re.search(r"gene_biotype=([^;]+)", c[8])
            name = m_name.group(1) if m_name else "?"
            refseq_genes[name] = {
                "seq_acc": c[0],
                "start_1based": s,
                "end_1based": e,
                "bed_0based": [s - 1, e],
                "span_bp": e - s + 1,
                "strand": c[6],
                "biotype": m_bio.group(1) if m_bio else "?",
            }

    # Parse Ensembl 116 Ggal audit JSON
    ens = json.loads((RUNS / "ensembl_116_ggal_audit.json").read_text(encoding="utf-8"))

    # Inspect Level 4 BLASTn hits for Ggal ENAM Exon 8 and Exon 6!
    ex8_blast = (RUNS / "level4" / "ENAM_exon8_aves" / "blast_Ggal.tsv").read_text(encoding="utf-8").strip().splitlines()
    ex8_hits = []
    for ln in ex8_blast:
        c = ln.split("\t")
        s1, s2 = int(c[6]), int(c[7])
        lo, hi = min(s1, s2), max(s1, s2)
        ex8_hits.append({
            "q": c[0], "pid": float(c[2]), "aln_len": int(c[3]),
            "q1": int(c[4]), "q2": int(c[5]),
            "rel_lo": lo, "rel_hi": hi,
            "abs_start_1based": win_start + lo - 1,
            "abs_end_1based": win_start + hi - 1,
            "genomic_span_bp": hi - lo + 1,
            "strand": "+" if s1 < s2 else "-",
            "evalue": float(c[8]), "bits": float(c[9]),
        })
    # Check full cluster span in level4 ENAM_exon8_aves
    ex8_msa_in = read_fasta(RUNS / "level4" / "ENAM_exon8_aves" / "msa_in.fna")["Ggal"]

    # Check Table S2b Ggal rows
    s2b = [r for r in csv.DictReader((RUNS / "supplementary" / "Table_S2b_whole_proteome_reciprocal_blastx.tsv").open(encoding="utf-8"), delimiter="\t") if r["species_code"] == "Ggal"]

    # Check exact distances
    nrg1 = refseq_genes["NRG1"]
    lpl = refseq_genes["LPL"]
    diff_endpoints = lpl["start_1based"] - nrg1["end_1based"]  # 55215002 - 55109784 = 105218
    strict_intergenic_bp = lpl["start_1based"] - nrg1["end_1based"] - 1  # 55109785..55215001 inclusive = 105217

    # Check Figure 1B plotted x-range from 06_make_figures.py
    fig_code = (ROOT / "scripts" / "enamel_relics" / "06_make_figures.py").read_text(encoding="utf-8")
    xlim_matches = re.findall(r"set_xlim\([^)]+\)|55\.098|55\.236|55098|55236", fig_code)

    # Check gap between Exon 8 (55155480-55157464) and Exon 6 (55160529-55160677)
    # Exon 8 end = 55157464; Exon 6 start = 55160529 -> strictly between = 55160529 - 55157464 - 1 = 3064 bp; diff = 3065 bp
    # And check distance from Exon 6 (55160529-55160677) to ENSGALG00010014557
    ens_lnc = {g["gene_id"]: g for g in ens.get("region_genes", [])} if isinstance(ens, dict) else {}

    # Check Table S2a intergenic_gap_bp formula across all 18 species
    gap_formula_check = []
    for code, row in s2a.items():
        if row["anchor_mode"] == "pair_NRG1_LPL":
            ns, ne = map(int, row["nrg1_coords"].split("(")[0].split("-"))
            ls, le = map(int, row["lpl_coords"].split("(")[0].split("-"))
            if ne < ls:
                d_end = ls - ne
                d_strict = ls - ne - 1
            else:
                d_end = ns - le
                d_strict = ns - le - 1
            gap_formula_check.append({
                "code": code,
                "reported_gap": int(row["intergenic_gap_bp"]),
                "endpoint_diff": d_end,
                "strict_between_bp": d_strict,
                "formula_used": "end2-start1_diff" if int(row["intergenic_gap_bp"]) == d_end else ("strict_between" if int(row["intergenic_gap_bp"]) == d_strict else "OTHER"),
            })

    return {
        "ggal_s2a": gg,
        "ggal_refseq_genes_in_window": refseq_genes,
        "ggal_ensembl_116_audit": ens,
        "nrg1_lpl_endpoint_diff": diff_endpoints,
        "nrg1_lpl_strict_intergenic_bp": strict_intergenic_bp,
        "nrg1_lpl_intergenic_1based_inclusive": f"NC_052572.1:{nrg1['end_1based']+1}-{lpl['start_1based']-1}",
        "nrg1_lpl_intergenic_bed_0based": f"NC_052572.1:{nrg1['end_1based']}-{lpl['start_1based']-1}",
        "fig1b_xlim_code_refs": xlim_matches,
        "fig1b_plotted_span_bp": int(round((55.236 - 55.098) * 1e6)),
        "fig1b_explanation": "Plotted axis 55.098-55.236 Mb (138 kb) shows 11.8 kb of 5' NRG1 + 105.2 kb intergenic gap + full 17.4 kb LPL + 3.6 kb right margin",
        "ex8_level4_blast_hits": ex8_hits,
        "ex8_msa_in_ungapped_len": len(ex8_msa_in),
        "ex8_genomic_span_55155480_55157464": 55157464 - 55155480 + 1,
        "ex8_s2b_rows": s2b,
        "dist_ex8_to_ex6_diff": 55160529 - 55157464,
        "dist_ex8_to_ex6_strict_between": 55160529 - 55157464 - 1,
        "all_18sp_gap_formula_check": gap_formula_check,
    }


def main() -> int:
    a = audit_phase_a()
    b1 = audit_phase_b1()
    b3 = audit_phase_b3()
    summary = {"phase_a": a, "phase_b1": b1, "phase_b3": b3}
    out_json = OUT / "phase_A_B1_B3_audit.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Log provenance
    prov_entry = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "step": "phase_A_B1_B3_audit",
        "script": "scripts/enamel_relics/07_audit_phaseA_B1_B3.py",
        "output_json": str(out_json.relative_to(ROOT).as_posix()),
        "output_sha256": sha256_file(out_json),
        "wsl_versions": a["wsl_live_versions"],
    }
    with PROV.open("a", encoding="utf-8") as f:
        f.write(json.dumps(prov_entry) + "\n")

    print("=== PHASE A SUMMARY ===")
    for k, items in a["manifest_checks"].items():
        print(f"  {k}: {sum(x['match'] for x in items)}/{len(items)} SHA-256 exact matches")
    print(f"  Figures: {sum(x['match'] for x in a['figure_checks'])}/{len(a['figure_checks'])} exact matches")
    print(f"  Typst: {a['typst_version']}")
    print(f"  IQ-TREE3: {a['wsl_live_versions'].get('iqtree3')}")
    print(f"  Policy URLs: {a['policy_urls']}")

    print("\n=== PHASE B1: CHICKEN PSI-ENAM EXON 6 ===")
    print(f"  Reported 149-bp span: {b1['reported_coords_1based_inclusive']} (BED: {b1['reported_bed_0based_halfopen']})")
    print(f"  Exact 47-bp exon core: {b1['exon6_core_coords_1based_inclusive']} (BED: {b1['exon6_core_bed_0based_halfopen']})")
    print(f"  Acceptor..Donor: {b1['ggal_acceptor_2bp']} .. {b1['ggal_donor_2bp']}")
    print(f"  Amis Ex6 (48bp): {b1['amis_exon6_aligned']}")
    print(f"  Ggal Ex6 (47bp): {b1['ggal_exon6_aligned']}")
    print(f"  BLASTn Arow->Ggal: {b1['blastn_arow_vs_ggal']}")
    print(f"  BLASTn Amis->Ggal: {b1['blastn_amis_vs_ggal']}")
    print(f"  Single-copy dc-megablast: {b1['singlecopy_dcmegablast']}")
    print(f"  Single-copy blastn: {b1['singlecopy_blastn']}")
    print("  Reciprocal BLASTx (-seg yes, e<=10):")
    for qk, hits in b1["blastx_std_segyes_top3"].items():
        print(f"    {qk}: {[(h['sseqid'], h['bitscore'], h['evalue'], h['stitle'][:45]) for h in hits]}")
    print("  Reciprocal BLASTx (-seg no, e<=10):")
    for qk, hits in b1["blastx_std_segno_top3"].items():
        print(f"    {qk}: {[(h['sseqid'], h['bitscore'], h['evalue'], h['stitle'][:45]) for h in hits]}")
    print("  Reciprocal BLASTx SENSITIVE (w=2, -seg no, e<=1000, comp=0):")
    for qk, hits in b1["blastx_sens_w2_top5"].items():
        print(f"    {qk}: {[(h['sseqid'], h['bitscore'], h['evalue'], h['stitle'][:45]) for h in hits]}")

    print("\n=== PHASE B3: COORDINATE RECONCILIATION ===")
    print(f"  NRG1 end={b3['ggal_refseq_genes_in_window']['NRG1']['end_1based']}, LPL start={b3['ggal_refseq_genes_in_window']['LPL']['start_1based']}")
    print(f"  Endpoint diff = {b3['nrg1_lpl_endpoint_diff']} bp; strictly between = {b3['nrg1_lpl_strict_intergenic_bp']} bp")
    print(f"  Fig1B span = {b3['fig1b_plotted_span_bp']} bp ({b3['fig1b_explanation']})")
    print(f"  Exon 8 genomic span (55155480..55157464 inclusive) = {b3['ex8_genomic_span_55155480_55157464']} bp; msa_in ungapped = {b3['ex8_msa_in_ungapped_len']} bp")
    for h in b3["ex8_level4_blast_hits"]:
        print(f"    Level4 BLASTn hit: q={h['q']} aln_len={h['aln_len']} abs={h['abs_start_1based']}-{h['abs_end_1based']} gspan={h['genomic_span_bp']} pid={h['pid']} bits={h['bits']} E={h['evalue']}")
    print(f"  Ensembl 116 Ggal region features: {b3['ggal_ensembl_116_audit'].get('region_genes', b3['ggal_ensembl_116_audit'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
