"""Phase 2.3 — Whole-proteome reciprocal BLASTx verification (GCF_030867095.1),
compositional/IDP bias quantification (-seg yes vs -seg no, 20 shuffles, SEG masking,
Shannon entropy across 18 species), and Level 3 miniprot lesion + canonical GT-AG
splice-site audit (full-length RefSeq ENAM XP_019346377.1 vs legacy QOX06214.1,
AMBN, AMTN, and positive control LPL).

Usage:
    uv run python scripts/enamel_relics/03_reciprocal_and_compositional.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs" / "enamel_relics"
RECIP_DIR = OUT / "reciprocal"
TABLES_DIR = OUT / "tables"

N_SHUFFLES = 20
STOPS = {"taa", "tag", "tga"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def append_provenance(event: dict[str, Any]) -> None:
    line = json.dumps(event, ensure_ascii=False) + "\n"
    with (OUT / "provenance_run.jsonl").open("a", encoding="utf-8") as f:
        f.write(line)
    with (ROOT / "runs" / "provenance.jsonl").open("a", encoding="utf-8") as f:
        f.write(line)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    cur = ""
    for ln in path.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            cur = ln[1:].split()[0]
            seqs[cur] = []
        elif cur:
            seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def parse_seg_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Parse segmasker -outfmt interval output (0-based inclusive start - end)."""
    out: dict[str, list[tuple[int, int]]] = {}
    cur = ""
    if not path.exists():
        return out
    for ln in path.read_text(encoding="utf-8").splitlines():
        s = ln.strip()
        if not s:
            continue
        if s.startswith(">"):
            cur = s[1:].split()[0]
            out.setdefault(cur, [])
        elif cur and " - " in s:
            a, b = s.split(" - ", 1)
            out[cur].append((int(a.strip()), int(b.strip())))
    return out


def seq_entropy(seq: str) -> float:
    n = len(seq)
    if n == 0:
        return 0.0
    counts = Counter(seq.upper())
    return -sum((c / n) * math.log2(c / n) for c in counts.values() if c > 0)


def compute_query_composition(meta: dict[str, Any]) -> list[dict[str, Any]]:
    seqs = read_fasta(OUT / "queries" / "all_real_queries.faa")
    segs = parse_seg_intervals(OUT / "segmask" / "seg_intervals.txt")
    rows: list[dict[str, Any]] = []
    for qname, qinfo in meta["queries"].items():
        seq = seqs.get(qname, "").upper()
        n = len(seq)
        if n == 0:
            continue
        masked = set()
        for a, b in segs.get(qname, []):
            for pos in range(a, min(b + 1, n)):
                masked.add(pos)
        c = Counter(seq)
        pct_p = 100.0 * c["P"] / n
        pct_q = 100.0 * c["Q"] / n
        pct_s = 100.0 * c["S"] / n
        pct_ed = 100.0 * (c["E"] + c["D"]) / n
        pct_disorder_scpp = 100.0 * sum(c[aa] for aa in "PQSDEGN") / n
        pct_seg = 100.0 * len(masked) / n
        rows.append({
            "query_name": qname,
            "accession": qinfo["accession"],
            "source_species": qinfo["source_species"],
            "role": qinfo["role"],
            "length_aa": n,
            "pct_Pro": round(pct_p, 2),
            "pct_Gln": round(pct_q, 2),
            "pct_Pro_Gln": round(pct_p + pct_q, 2),
            "pct_Ser": round(pct_s, 2),
            "pct_Glu_Asp": round(pct_ed, 2),
            "pct_SCPP_disorder_PQSDEGN": round(pct_disorder_scpp, 2),
            "shannon_entropy_bits": round(seq_entropy(seq), 3),
            "seg_masked_aa": len(masked),
            "pct_seg_masked": round(pct_seg, 2),
        })
    return rows


def analyze_tblastn_and_nulls(meta: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Compute per-species per-query obs, null mean/p95/max, excess for BOTH -seg yes and -seg no,
    plus retention indices across all 18 species."""
    sp_order = [r["code"] for r in meta["enamel_intervals"]]
    asm_map = meta["assemblies"]

    per_gene_rows: list[dict[str, Any]] = []
    index_rows: list[dict[str, Any]] = []

    enamel_qnames = [k for k in meta["queries"] if any(k.startswith(p) for p in ("ENAM_", "AMBN_", "AMTN_", "LPL_", "DSPP_"))]
    idp_qnames = [k for k in meta["queries"] if k not in enamel_qnames]

    # Cache best bitscore per (domain, seg_mode, sp, qid)
    best_cache: dict[tuple[str, str, str], dict[str, float]] = {}
    for domain, folder_yes, folder_no in [
        ("enamel", "tblastn_enamel_seg_yes", "tblastn_enamel_seg_no"),
        ("idp", "tblastn_idp_seg_yes", "tblastn_idp_seg_no"),
    ]:
        for seg_mode, folder in [("seg_yes", folder_yes), ("seg_no", folder_no)]:
            for sp in sp_order:
                p = OUT / folder / f"{sp}.tsv"
                best: dict[str, float] = {}
                if p.exists():
                    with p.open(encoding="utf-8") as f:
                        for c in csv.reader(f, delimiter="\t"):
                            if len(c) >= 12:
                                qid, bs = c[0], float(c[11])
                                if bs > best.get(qid, 0.0):
                                    best[qid] = bs
                best_cache[(domain, seg_mode, sp)] = best

    # Per-gene rows
    for sp in sp_order:
        sp_meta = asm_map[sp]
        for seg_mode in ("seg_yes", "seg_no"):
            for domain, qlist in [("enamel", enamel_qnames), ("idp", idp_qnames)]:
                best = best_cache[(domain, seg_mode, sp)]
                for q in qlist:
                    obs = best.get(q, 0.0)
                    shufs = [best.get(f"{q}~shuf{i:02d}", 0.0) for i in range(N_SHUFFLES)]
                    n_mean = float(np.mean(shufs))
                    n_p95 = float(np.percentile(shufs, 95))
                    n_max = float(np.max(shufs))
                    exc = max(0.0, obs - n_p95)
                    per_gene_rows.append({
                        "species_code": sp,
                        "scientific_name": sp_meta["scientific_name"],
                        "clade": sp_meta["clade"],
                        "dentition": sp_meta["dentition"],
                        "domain": domain,
                        "seg_mode": seg_mode,
                        "query_name": q,
                        "query_source": meta["queries"][q]["source_species"],
                        "query_role": meta["queries"][q]["role"],
                        "obs_bits": round(obs, 2),
                        "null_mean_bits": round(n_mean, 2),
                        "null_p95_bits": round(n_p95, 2),
                        "null_max_bits": round(n_max, 2),
                        "excess_over_p95": round(exc, 2),
                    })

    # Retention index sets:
    index_configs = [
        ("Amis_refseq_full", "Amis", "ENAM_Amis", "AMBN_Amis", "AMTN_Amis", "LPL_Amis"),
        ("Amis_legacy", "Amis", "ENAM_Amis_legacy", "AMBN_Amis", "AMTN_Amis_gb", "LPL_Amis"),
        ("Acar_equidistant", "Acar", "ENAM_Acar", "AMBN_Acar", "AMTN_Acar", "LPL_Acar"),
    ]
    lookup = {(r["species_code"], r["seg_mode"], r["query_name"]): r for r in per_gene_rows}
    for cfg_name, src, q_enam, q_ambn, q_amtn, q_lpl in index_configs:
        for seg_mode in ("seg_yes", "seg_no"):
            for sp in sp_order:
                sp_meta = asm_map[sp]
                re_ = lookup[(sp, seg_mode, q_enam)]
                ra_ = lookup[(sp, seg_mode, q_ambn)]
                rt_ = lookup[(sp, seg_mode, q_amtn)]
                rl_ = lookup[(sp, seg_mode, q_lpl)]
                raw_idx = (re_["obs_bits"] + ra_["obs_bits"] + rt_["obs_bits"]) / rl_["obs_bits"] if rl_["obs_bits"] > 0 else float("nan")
                exc_idx = (re_["excess_over_p95"] + ra_["excess_over_p95"] + rt_["excess_over_p95"]) / rl_["excess_over_p95"] if rl_["excess_over_p95"] > 0 else float("nan")
                index_rows.append({
                    "query_set": cfg_name,
                    "query_source": src,
                    "seg_mode": seg_mode,
                    "species_code": sp,
                    "scientific_name": sp_meta["scientific_name"],
                    "clade": sp_meta["clade"],
                    "dentition": sp_meta["dentition"],
                    "is_self_query": sp == src,
                    "ENAM_obs": re_["obs_bits"],
                    "ENAM_p95": re_["null_p95_bits"],
                    "ENAM_exc": re_["excess_over_p95"],
                    "AMBN_obs": ra_["obs_bits"],
                    "AMBN_p95": ra_["null_p95_bits"],
                    "AMBN_exc": ra_["excess_over_p95"],
                    "AMTN_obs": rt_["obs_bits"],
                    "AMTN_p95": rt_["null_p95_bits"],
                    "AMTN_exc": rt_["excess_over_p95"],
                    "LPL_obs": rl_["obs_bits"],
                    "LPL_p95": rl_["null_p95_bits"],
                    "LPL_exc": rl_["excess_over_p95"],
                    "raw_index": round(raw_idx, 3),
                    "null_corrected_index": round(exc_idx, 3),
                })

    return per_gene_rows, index_rows


def run_whole_proteome_reciprocal(meta: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Extract candidate loci (bitscore >= 28 or top clustered HSPs) from tblastn_enamel_seg_yes
    across all 18 species and run BLASTx against the ENTIRE Alligator proteome (~48k proteins)."""
    RECIP_DIR.mkdir(parents=True, exist_ok=True)
    interval_map = {r["code"]: r for r in meta["enamel_intervals"]}
    targets = ["ENAM_Amis", "ENAM_Amis_legacy", "AMBN_Amis", "AMTN_Amis", "AMTN_Amis_gb"]
    expected_keywords = {
        "ENAM_Amis": ("enamelin", "XP_019346377", "XP_019346375", "XP_019346376"),
        "ENAM_Amis_legacy": ("enamelin", "XP_019346377", "XP_019346375", "XP_019346376"),
        "AMBN_Amis": ("ameloblastin", "XP_019346386", "XP_019346383", "XP_019346384", "XP_019346385"),
        "AMTN_Amis": ("amelotin", "XP_014451634"),
        "AMTN_Amis_gb": ("amelotin", "XP_014451634"),
    }

    cand_meta: dict[str, dict[str, Any]] = {}
    cand_fasta_lines: list[str] = []
    top_hsp_fasta_lines: list[str] = []

    for sp, imeta in interval_map.items():
        tsv = OUT / "tblastn_enamel_seg_yes" / f"{sp}.tsv"
        seq_dict = read_fasta(ROOT / imeta["fasta"])
        seq = next(iter(seq_dict.values()))
        win_start = int(imeta["win_start"])
        seq_acc = imeta["seq_acc"]

        hsps_by_q: dict[str, list[dict[str, Any]]] = {t: [] for t in targets}
        with tsv.open(encoding="utf-8") as f:
            for c in csv.reader(f, delimiter="\t"):
                if len(c) >= 12 and c[0] in hsps_by_q:
                    s1, s2 = int(c[8]), int(c[9])
                    hsps_by_q[c[0]].append({
                        "pid": float(c[2]),
                        "alen": int(c[3]),
                        "qstart": int(c[6]),
                        "qend": int(c[7]),
                        "s_lo": min(s1, s2),
                        "s_hi": max(s1, s2),
                        "strand": "+" if s1 < s2 else "-",
                        "evalue": float(c[10]),
                        "bits": float(c[11]),
                    })

        for q in targets:
            hsps = [h for h in hsps_by_q[q] if h["bits"] >= 28.0 or h["evalue"] <= 1e-2]
            if not hsps:
                continue
            # Cluster HSPs within 15 kb of the top-scoring HSP on the same strand
            hsps.sort(key=lambda h: -h["bits"])
            top = hsps[0]
            cluster = [h for h in hsps if h["strand"] == top["strand"] and abs(h["s_lo"] - top["s_lo"]) <= 15_000]
            c_lo = min(h["s_lo"] for h in cluster)
            c_hi = max(h["s_hi"] for h in cluster)
            pad_lo = max(1, c_lo - 100)
            pad_hi = min(len(seq), c_hi + 100)
            cid = f"{sp}__{q}"
            subseq = seq[pad_lo - 1 : pad_hi]
            cand_fasta_lines.append(f">{cid}\n{subseq}\n")
            # Also extract the single top-scoring HSP (+-30 bp) so multi-kb introns cannot
            # mask a short orthologous exon in whole-proteome reciprocal BLASTx
            t_lo = max(1, top["s_lo"] - 30)
            t_hi = min(len(seq), top["s_hi"] + 30)
            top_hsp_fasta_lines.append(f">{cid}\n{seq[t_lo - 1 : t_hi]}\n")
            cand_meta[cid] = {
                "candidate_id": cid,
                "species_code": sp,
                "query_name": q,
                "seq_acc": seq_acc,
                "seq_name": imeta["seq_name"],
                "strand": top["strand"],
                "n_clustered_hsps": len(cluster),
                "total_hsps_ge28": len(hsps),
                "top_tblastn_pid": top["pid"],
                "top_tblastn_len": top["alen"],
                "top_tblastn_qspan": f"{min(h['qstart'] for h in cluster)}-{max(h['qend'] for h in cluster)}",
                "top_tblastn_evalue": top["evalue"],
                "top_tblastn_bits": top["bits"],
                "abs_start": win_start + c_lo - 1,
                "abs_end": win_start + c_hi - 1,
                "span_bp": c_hi - c_lo + 1,
                "extracted_start": win_start + pad_lo - 1,
                "extracted_end": win_start + pad_hi - 1,
                "top_hsp_start": win_start + top["s_lo"] - 1,
                "top_hsp_end": win_start + top["s_hi"] - 1,
            }

    cand_fa = RECIP_DIR / "candidates.fna"
    cand_fa.write_text("".join(cand_fasta_lines), encoding="utf-8", newline="\n")
    top_fa = RECIP_DIR / "candidates_top_hsp.fna"
    top_fa.write_text("".join(top_hsp_fasta_lines), encoding="utf-8", newline="\n")

    # Run whole-proteome BLASTx in WSL (-seg yes AND -seg no, full-span + top-HSP) + chicken single-copy dc-megablast
    if (
        not (RECIP_DIR / "blastx_proteome_seg_yes.tsv").exists()
        or not (RECIP_DIR / "blastx_top_hsp_seg_yes.tsv").exists()
        or not (RECIP_DIR / "ggal_single_copy.tsv").exists()
    ):
        wsl_out = "/mnt/c/Users/user/bio-stack/runs/enamel_relics"
        script_sh = RECIP_DIR / "run_reciprocal.sh"
        script_sh.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            f"D={wsl_out}\n"
            'FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle"\n'
            'blastx -query "$D/reciprocal/candidates.fna" -db "$D/blastdb/amis_proteome" '
            '-evalue 10 -seg yes -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" '
            '> "$D/reciprocal/blastx_proteome_seg_yes.tsv"\n'
            'blastx -query "$D/reciprocal/candidates.fna" -db "$D/blastdb/amis_proteome" '
            '-evalue 10 -seg no -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" '
            '> "$D/reciprocal/blastx_proteome_seg_no.tsv"\n'
            'blastx -query "$D/reciprocal/candidates_top_hsp.fna" -db "$D/blastdb/amis_proteome" '
            '-evalue 10 -seg yes -max_target_seqs 25 -num_threads 8 -outfmt "$FMT" '
            '> "$D/reciprocal/blastx_top_hsp_seg_yes.tsv"\n'
            '# Single-copy test for Gallus gallus psi-ENAM and psi-AMBN against Ggal NRG1-LPL interval\n'
            'awk \'/^>Ggal__ENAM_Amis$|^>Ggal__AMBN_Amis$/{p=1; print; next} /^>/{p=0} p\' '
            '"$D/reciprocal/candidates.fna" > "$D/reciprocal/ggal_relics.fna"\n'
            'blastn -task dc-megablast -query "$D/reciprocal/ggal_relics.fna" '
            '-subject "$D/intervals_enamel/interval_Ggal.fna" -evalue 1e-3 -dust yes '
            '-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" '
            '> "$D/reciprocal/ggal_single_copy.tsv"\n',
            encoding="utf-8",
            newline="\n",
        )
        res = subprocess.run(
            ["wsl", "-d", "Ubuntu-26.04", "--", "bash", "-l", f"{wsl_out}/reciprocal/run_reciprocal.sh"],
            capture_output=True,
            text=True,
            check=False,
        )
        if res.returncode != 0:
            raise RuntimeError(f"Reciprocal BLASTx failed:\n{res.stderr}\n{res.stdout}")

    # Parse whole-proteome BLASTx results (full span + top HSP)
    def load_blastx(path: Path) -> dict[str, list[dict[str, Any]]]:
        out: dict[str, list[dict[str, Any]]] = {}
        with path.open(encoding="utf-8") as f:
            for c in csv.reader(f, delimiter="\t"):
                if len(c) >= 11:
                    out.setdefault(c[0], []).append({
                        "sseqid": c[1],
                        "pid": float(c[2]),
                        "alen": int(c[3]),
                        "evalue": float(c[8]),
                        "bits": float(c[9]),
                        "stitle": c[10],
                    })
        return out

    hits_by_cid = load_blastx(RECIP_DIR / "blastx_proteome_seg_yes.tsv")
    top_hsp_hits = load_blastx(RECIP_DIR / "blastx_top_hsp_seg_yes.tsv")

    recip_rows: list[dict[str, Any]] = []
    for cid, cm in cand_meta.items():
        q = cm["query_name"]
        kws = expected_keywords[q]
        hits = sorted(hits_by_cid.get(cid, []), key=lambda h: (-h["bits"], h["evalue"]))
        thits = sorted(top_hsp_hits.get(cid, []), key=lambda h: (-h["bits"], h["evalue"]))

        def is_ortholog(h: dict[str, Any]) -> bool:
            low = (h["sseqid"] + " " + h["stitle"]).lower()
            return any(kw.lower() in low for kw in kws)

        top_hit = hits[0] if hits else None
        runner_up = next((h for h in hits if not is_ortholog(h)), None)
        confirmed = bool(top_hit and is_ortholog(top_hit))

        th_top = thits[0] if thits else None
        th_ru = next((h for h in thits if not is_ortholog(h)), None)
        th_confirmed = bool(th_top and is_ortholog(th_top))

        recip_rows.append({
            **cm,
            "proteome_size": 48040,  # overwritten below with exact count
            "reciprocal_confirmed_span": confirmed,
            "top_hsp_reciprocal_confirmed": th_confirmed,
            "top_proteome_acc": top_hit["sseqid"] if top_hit else "NONE",
            "top_proteome_title": top_hit["stitle"][:90] if top_hit else "NONE",
            "top_proteome_pid": top_hit["pid"] if top_hit else "",
            "top_proteome_len": top_hit["alen"] if top_hit else "",
            "top_proteome_evalue": top_hit["evalue"] if top_hit else "",
            "top_proteome_bits": top_hit["bits"] if top_hit else "",
            "runnerup_nonortholog_acc": runner_up["sseqid"] if runner_up else "NONE",
            "runnerup_nonortholog_title": runner_up["stitle"][:80] if runner_up else "NONE",
            "runnerup_nonortholog_evalue": runner_up["evalue"] if runner_up else "",
            "runnerup_nonortholog_bits": runner_up["bits"] if runner_up else "",
            "top_hsp_proteome_acc": th_top["sseqid"] if th_top else "NONE",
            "top_hsp_proteome_evalue": th_top["evalue"] if th_top else "",
            "top_hsp_proteome_bits": th_top["bits"] if th_top else "",
            "top_hsp_runnerup_acc": th_ru["sseqid"] if th_ru else "NONE",
            "top_hsp_runnerup_bits": th_ru["bits"] if th_ru else "",
        })

    # Exact proteome count
    n_prot = sum(1 for ln in (OUT / "refdata" / "Amis_protein.faa").open(encoding="utf-8") if ln.startswith(">"))
    for r in recip_rows:
        r["proteome_size"] = n_prot

    # Parse chicken single-copy test
    ggal_sc: dict[str, int] = Counter()
    with (RECIP_DIR / "ggal_single_copy.tsv").open(encoding="utf-8") as f:
        for c in csv.reader(f, delimiter="\t"):
            if c:
                ggal_sc[c[0]] += 1

    return recip_rows, {"proteome_size": n_prot, "ggal_single_copy_hits": dict(ggal_sc)}


@dataclass
class MpAln:
    qlen: int
    qini: int
    qfim: int
    strand: str
    sini: int
    sfim: int
    score: int
    fs: int
    st: int
    da: int
    do: int
    identity: float
    lesions: list[tuple[str, int]] = field(default_factory=list)
    introns: list[tuple[str, int, str]] = field(default_factory=list)  # (donor, len, acceptor)

    @property
    def coverage(self) -> float:
        return (self.qfim - self.qini) / self.qlen if self.qlen else 0.0


def parse_miniprot_paf(path: Path, max_intron_bp: int = 25_000) -> list[MpAln]:
    """Parse miniprot PAF alignments. If an alignment contains an implausibly large
    intron (> max_intron_bp, e.g., Cpic ENAM 58-kb spurious jump to protein SON),
    split/trim to the primary contiguous orthologous block (largest non-jump segment)
    so distal low-complexity chaining artifacts never inflate lesion counts."""
    alns: list[MpAln] = []
    if not path.exists():
        return alns
    for ln in path.read_text(encoding="utf-8").splitlines():
        c = ln.split("\t")
        if len(c) < 12:
            continue
        tags = {t[:2]: t[5:] for t in c[12:] if len(t) > 5}
        cg, cs = tags.get("cg", ""), tags.get("cs", "")
        strand = c[4]
        qlen, qini, qfim = int(c[1]), int(c[2]), int(c[3])
        sini, sfim = int(c[7]), int(c[8])
        # Walk cg & cs, segmenting at any intron > max_intron_bp
        cg_ops = [(int(n_s), op) for n_s, op in re.findall(r"(\d+)([A-Z])", cg)]
        cs_toks = re.findall(r":\d+|\*[a-z]{3}[A-Z*]|\+[A-Z*]+|-[a-z]+|~[a-z]{2}\d+[a-z]{2}", cs)

        # Segment cg_ops by introns > max_intron_bp
        cg_segs: list[tuple[int, int, int, int, list[tuple[int, int]]]] = []
        # track (q_start, q_end, t_lo, t_hi, fs_list)
        qpos = qini
        tpos = sini if strand == "+" else sfim
        sgn = 1 if strand == "+" else -1
        seg_q0, seg_t0 = qpos, tpos
        seg_fs: list[tuple[str, int]] = []
        for n, op in cg_ops:
            if op in "NUV" and n > max_intron_bp:
                cg_segs.append((seg_q0, qpos, min(seg_t0, tpos), max(seg_t0, tpos), list(seg_fs)))
                if op == "N":
                    tpos += sgn * n
                elif op == "U":
                    tpos += sgn * n
                    qpos += 1
                elif op == "V":
                    tpos += sgn * n
                    qpos += 1
                seg_q0, seg_t0 = qpos, tpos
                seg_fs = []
                continue
            if op == "M":
                qpos += n
                tpos += sgn * 3 * n
            elif op == "D":
                tpos += sgn * 3 * n
            elif op == "I":
                qpos += n
            elif op == "F":
                seg_fs.append(("fs", qpos))
                tpos += sgn * n
            elif op == "G":
                seg_fs.append(("fs", qpos))
                qpos += 1
                tpos += sgn * n
            elif op == "N":
                tpos += sgn * n
            elif op in "UV":
                qpos += 1
                tpos += sgn * n
        cg_segs.append((seg_q0, qpos, min(seg_t0, tpos), max(seg_t0, tpos), list(seg_fs)))

        # Pick largest segment by query aa span
        best_seg_idx = max(range(len(cg_segs)), key=lambda idx: cg_segs[idx][1] - cg_segs[idx][0])
        b_q0, b_q1, b_tlo, b_thi, b_fs = cg_segs[best_seg_idx]

        # Segment cs_toks identically on introns > max_intron_bp
        cs_segs: list[tuple[int, int, list[tuple[str, int]], list[tuple[str, int, str]]]] = []
        pos = qini
        eq, sub = 0, 0
        seg_st: list[tuple[str, int]] = []
        seg_introns: list[tuple[str, int, str]] = []
        for tok in cs_toks:
            if tok.startswith("~"):
                m = re.match(r"~([a-z]{2})(\d+)([a-z]{2})", tok)
                if m:
                    ilen = int(m.group(2))
                    if ilen > max_intron_bp:
                        cs_segs.append((eq, sub, list(seg_st), list(seg_introns)))
                        eq, sub, seg_st, seg_introns = 0, 0, [], []
                        continue
                    seg_introns.append((m.group(1).upper(), ilen, m.group(3).upper()))
            elif tok.startswith(":"):
                val = int(tok[1:])
                eq += val
                pos += val
            elif tok.startswith("*"):
                if tok[1:4] in STOPS:
                    seg_st.append(("stop", pos))
                sub += 1
                pos += 1
            elif tok.startswith("+"):
                pos += len(tok) - 1
        cs_segs.append((eq, sub, list(seg_st), list(seg_introns)))

        # Match best segment
        if best_seg_idx < len(cs_segs):
            b_eq, b_sub, b_st, b_introns = cs_segs[best_seg_idx]
        else:
            b_eq, b_sub, b_st, b_introns = cs_segs[-1]

        lesions = sorted(b_fs + b_st, key=lambda x: x[1])
        alns.append(
            MpAln(
                qlen=qlen,
                qini=b_q0,
                qfim=b_q1,
                strand=strand,
                sini=b_tlo,
                sfim=b_thi,
                score=int(tags.get("AS", 0)),
                fs=len(b_fs),
                st=len(b_st),
                da=int(tags.get("da", 0)),
                do=int(tags.get("do", 0)),
                identity=b_eq / (b_eq + b_sub) if (b_eq + b_sub) else 0.0,
                lesions=lesions,
                introns=b_introns,
            )
        )
    return alns


def audit_miniprot_level3(meta: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    interval_map = {r["code"]: r for r in meta["enamel_intervals"]}
    asm_map = meta["assemblies"]
    genes = ["ENAM_Amis", "ENAM_Amis_legacy", "AMBN_Amis", "AMTN_Amis", "AMTN_Amis_gb", "LPL_Amis"]
    rows: list[dict[str, Any]] = []

    for g in genes:
        for sp, imeta in interval_map.items():
            sp_meta = asm_map[sp]
            paf = OUT / "miniprot" / f"mp_{g}_{sp}.paf"
            alns = parse_miniprot_paf(paf)
            if not alns:
                rows.append({
                    "query_gene": g,
                    "species_code": sp,
                    "scientific_name": sp_meta["scientific_name"],
                    "clade": sp_meta["clade"],
                    "dentition": sp_meta["dentition"],
                    "aligned": False,
                    "q_ini": "", "q_fim": "", "q_len": "", "aligned_aa": 0,
                    "coverage": 0.0, "identity": 0.0, "score_AS": 0,
                    "frameshifts": 0, "stops": 0, "total_coding_lesions": 0,
                    "n_introns": 0, "n_gtag": 0, "n_gcag": 0, "n_noncanonical_introns": 0,
                    "da_tag": 0, "do_tag": 0, "noncanonical_splice_sites": "",
                    "seq_acc": imeta["seq_acc"], "strand": "", "abs_start": "", "abs_end": "",
                    "lesion_coords": "",
                })
                continue
            a = max(alns, key=lambda x: x.score)
            win_start = int(imeta["win_start"])
            n_gtag = sum(1 for d, _, ac in a.introns if d == "GT" and ac == "AG")
            n_gcag = sum(1 for d, _, ac in a.introns if d == "GC" and ac == "AG")
            noncanon = [f"{d}..{ac}({ln}bp)" for d, ln, ac in a.introns if (d, ac) not in (("GT", "AG"), ("GC", "AG"))]
            rows.append({
                "query_gene": g,
                "species_code": sp,
                "scientific_name": sp_meta["scientific_name"],
                "clade": sp_meta["clade"],
                "dentition": sp_meta["dentition"],
                "aligned": True,
                "q_ini": a.qini + 1,
                "q_fim": a.qfim,
                "q_len": a.qlen,
                "aligned_aa": a.qfim - a.qini,
                "coverage": round(a.coverage, 3),
                "identity": round(a.identity, 3),
                "score_AS": a.score,
                "frameshifts": a.fs,
                "stops": a.st,
                "total_coding_lesions": a.fs + a.st,
                "n_introns": len(a.introns),
                "n_gtag": n_gtag,
                "n_gcag": n_gcag,
                "n_noncanonical_introns": len(noncanon),
                "da_tag": a.da,
                "do_tag": a.do,
                "noncanonical_splice_sites": ";".join(noncanon),
                "seq_acc": imeta["seq_acc"],
                "strand": a.strand,
                "abs_start": win_start + a.sini,
                "abs_end": win_start + a.sfim - 1,
                "lesion_coords": ";".join(f"{t}@{p+1}" for t, p in a.lesions),
            })

    # Clean summary stats (strictly primary RefSeq cluster: ENAM_Amis + AMBN_Amis + AMTN_Amis, vs LPL_Amis positive control)
    primary_genes = {"ENAM_Amis", "AMBN_Amis", "AMTN_Amis"}
    summary: dict[str, Any] = {}
    for label, gset, dent, only_lesioned in [
        ("cluster_toothed_all", primary_genes, "toothed", False),
        ("cluster_edentulous_with_lesions", primary_genes, "edentulous", True),
        ("cluster_edentulous_all_aligned", primary_genes, "edentulous", False),
        ("LPL_ctrl_toothed", {"LPL_Amis"}, "toothed", False),
        ("LPL_ctrl_edentulous", {"LPL_Amis"}, "edentulous", False),
    ]:
        sub = [
            r for r in rows
            if r["query_gene"] in gset and r["dentition"] == dent and r["aligned"]
            and (not only_lesioned or r["total_coding_lesions"] > 0)
        ]
        tot_aa = sum(int(r["aligned_aa"]) for r in sub)
        tot_les = sum(int(r["total_coding_lesions"]) for r in sub)
        tot_fs = sum(int(r["frameshifts"]) for r in sub)
        tot_st = sum(int(r["stops"]) for r in sub)
        tot_int = sum(int(r["n_introns"]) for r in sub)
        tot_gtag = sum(int(r["n_gtag"]) for r in sub)
        tot_noncan = sum(int(r["n_noncanonical_introns"]) for r in sub)
        summary[label] = {
            "n_alignments": len(sub),
            "aligned_aa": tot_aa,
            "frameshifts": tot_fs,
            "stops": tot_st,
            "total_coding_lesions": tot_les,
            "lesions_per_1000aa": round(1000.0 * tot_les / tot_aa, 2) if tot_aa else 0.0,
            "total_introns": tot_int,
            "canonical_gtag_introns": tot_gtag,
            "noncanonical_introns": tot_noncan,
        }
    return rows, summary


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    t0 = datetime.now(timezone.utc)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    meta = json.loads((OUT / "audit_metadata.json").read_text(encoding="utf-8"))

    print("=== [1/4] Query compositional & SEG masking profile ===")
    comp_rows = compute_query_composition(meta)
    write_tsv(TABLES_DIR / "query_composition_seg.tsv", comp_rows)
    for r in comp_rows:
        print(f"  {r['query_name']:18s} len={r['length_aa']:4d} Pro+Gln={r['pct_Pro_Gln']:5.1f}% "
              f"SCPP_dis={r['pct_SCPP_disorder_PQSDEGN']:5.1f}% H={r['shannon_entropy_bits']:.2f} "
              f"SEG={r['pct_seg_masked']:5.1f}%")

    print("\n=== [2/4] tBLASTn SEG-on vs SEG-off & 20-shuffle null analysis (18 species) ===")
    per_gene_rows, index_rows = analyze_tblastn_and_nulls(meta)
    write_tsv(TABLES_DIR / "tblastn_per_gene_seg_comparison.tsv", per_gene_rows)
    write_tsv(TABLES_DIR / "retention_indices_18sp.tsv", index_rows)

    # Print key IDP vs Globular summary for equidistant Anolis query
    print("\n  --- Equidistant Anolis query: Clade medians/ranges (SEG YES vs SEG NO) ---")
    for seg_mode in ("seg_yes", "seg_no"):
        print(f"  [{seg_mode}]")
        for q in ("LPL_Acar", "IBSP_Acar", "SPP1_Acar", "ENAM_Acar", "AMBN_Acar", "AMTN_Acar"):
            parts = []
            for cl in ("Crocodylia", "Testudines", "Aves"):
                sub = [r for r in per_gene_rows if r["seg_mode"] == seg_mode and r["query_name"] == q and r["clade"] == cl]
                excs = [r["excess_over_p95"] for r in sub]
                p95s = [r["null_p95_bits"] for r in sub]
                parts.append(f"{cl[:4]}={min(excs):.1f}..{max(excs):.1f}(med {np.median(excs):.1f}, n95={np.median(p95s):.1f})")
            print(f"    {q:10s}: " + " | ".join(parts))

    print("\n=== [3/4] Whole-proteome reciprocal BLASTx (GCF_030867095.1) ===")
    recip_rows, recip_summary = run_whole_proteome_reciprocal(meta)
    write_tsv(TABLES_DIR / "whole_proteome_reciprocal.tsv", recip_rows)
    print(f"  Proteome size: {recip_summary['proteome_size']} proteins | Ggal single-copy: {recip_summary['ggal_single_copy_hits']}")
    for r in recip_rows:
        if r["species_code"] in ("Ggal", "Shab", "Dnov", "Apla", "Tgut", "Scam", "Arow", "Tgua", "Cpic", "Cmyd", "Gevg", "Psin") and r["query_name"] in ("ENAM_Amis", "AMBN_Amis", "AMTN_Amis"):
            print(f"  {r['species_code']:4s} {r['query_name']:10s} {r['seq_acc']}:{r['abs_start']}-{r['abs_end']}({r['strand']}) "
                  f"q={r['top_tblastn_qspan']} tBLASTn={r['top_tblastn_bits']}b(e={r['top_tblastn_evalue']:.2g}) "
                  f"RECIP_SPAN={r['reciprocal_confirmed_span']} RECIP_HSP={r['top_hsp_reciprocal_confirmed']} "
                  f"hsp_top={r['top_hsp_proteome_acc']}({r['top_hsp_proteome_bits']}b) "
                  f"hsp_next={r['top_hsp_runnerup_acc']}({r['top_hsp_runnerup_bits']}b)")

    print("\n=== [4/4] Level 3 miniprot coding lesions & GT-AG splice-site audit ===")
    mp_rows, mp_summary = audit_miniprot_level3(meta)
    write_tsv(TABLES_DIR / "miniprot_level3_splice_audit.tsv", mp_rows)
    for k, v in mp_summary.items():
        print(f"  {k:32s}: n={v['n_alignments']:2d} aa={v['aligned_aa']:5d} fs={v['frameshifts']:2d} st={v['stops']:2d} "
              f"rate={v['lesions_per_1000aa']:5.2f}/1000 | introns={v['total_introns']} GT-AG={v['canonical_gtag_introns']} noncanon={v['noncanonical_introns']}")

    print("\n  --- Newly aligned or lesioned miniprot hits (ENAM_Amis full vs legacy, AMBN_Amis, AMTN_Amis) ---")
    for r in mp_rows:
        if r["aligned"] and r["dentition"] == "edentulous" and r["query_gene"] != "LPL_Amis":
            print(f"  {r['query_gene']:16s} {r['species_code']:4s} q={r['q_ini']}-{r['q_fim']}/{r['q_len']} "
                  f"cov={r['coverage']:.2f} id={r['identity']:.2f} fs={r['frameshifts']} st={r['stops']} "
                  f"introns={r['n_introns']}(GT-AG={r['n_gtag']}, noncan={r['n_noncanonical_introns']}:{r['noncanonical_splice_sites']}) "
                  f"loc={r['seq_acc']}:{r['abs_start']}-{r['abs_end']}({r['strand']})")

    summary_json = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reciprocal_summary": recip_summary,
        "miniprot_summary": mp_summary,
    }
    (OUT / "phase2_step3_summary.json").write_text(json.dumps(summary_json, indent=2), encoding="utf-8")

    append_provenance({
        "ts": datetime.now(timezone.utc).isoformat(),
        "kind": "enamel-relics-phase2-step3",
        "script": "scripts/enamel_relics/03_reciprocal_and_compositional.py",
        "outputs": {
            "query_composition_seg.tsv": sha256_file(TABLES_DIR / "query_composition_seg.tsv"),
            "retention_indices_18sp.tsv": sha256_file(TABLES_DIR / "retention_indices_18sp.tsv"),
            "whole_proteome_reciprocal.tsv": sha256_file(TABLES_DIR / "whole_proteome_reciprocal.tsv"),
            "miniprot_level3_splice_audit.tsv": sha256_file(TABLES_DIR / "miniprot_level3_splice_audit.tsv"),
        },
        "elapsed_s": round((datetime.now(timezone.utc) - t0).total_seconds(), 2),
    })
    print(f"\nDone in {(datetime.now(timezone.utc) - t0).total_seconds():.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
