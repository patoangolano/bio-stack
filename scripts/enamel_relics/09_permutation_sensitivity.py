"""Phase C: High-permutation sensitivity benchmark (N = 20, 200, 500, 1000) + block-shuffle null
+ length/composition/SEG/uncertainty/LPL-index statistical re-evaluation across all 18 sauropsids.

Deterministic seeds:
- Global composition-preserving shuffles (N=1000 per query): numpy.random.default_rng(20261119)
- Local 20-aa block-preserving shuffles (N=200 per equidistant Anolis query): numpy.random.default_rng(20261120)
- Bootstrap CIs (B=2000): numpy.random.default_rng(20261121)

Outputs in runs/enamel_relics/v101/stats/:
- perm_sensitivity_per_gene_v101.tsv
- perm_sensitivity_summary_v101.tsv
- retention_indices_1000perm_v101.tsv
- clade_medians_equidistant_1000perm_v101.tsv
- regression_length_composition_v101.json
- phaseC_summary.json
- runs/enamel_relics/v101/provenance_phaseC.jsonl
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
OUT = RUNS / "v101" / "stats"
OUT.mkdir(parents=True, exist_ok=True)
PROV = RUNS / "v101" / "provenance_phaseC.jsonl"

WSL_DISTRO = "Ubuntu-26.04"
WSL_ROOT = "/mnt/c/Users/user/bio-stack"

N_MAX = 1000
N_LEVELS = (20, 200, 500, 1000)
N_BLOCK = 200
SEED_GLOBAL = 20261119
SEED_BLOCK = 20261120
SEED_BOOT = 20261121

ENAMEL_QUERIES = [
    "ENAM_Acar", "AMBN_Acar", "AMTN_Acar", "LPL_Acar",
    "ENAM_Amis", "AMBN_Amis", "AMTN_Amis", "LPL_Amis", "DSPP_Amis",
]
IDP_QUERIES = [
    "IBSP_Acar", "SPP1_Acar",
    "IBSP_Amis", "SPP1_Amis", "SPARCL1_Amis",
]
BLOCK_NULL_QUERIES_ENAMEL = ["ENAM_Acar", "AMBN_Acar", "AMTN_Acar", "LPL_Acar"]
BLOCK_NULL_QUERIES_IDP = ["IBSP_Acar", "SPP1_Acar", "SPARCL1_Amis"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def write_fasta(path: Path, seqs: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for name, seq in seqs:
            f.write(f">{name}\n{seq}\n")


def block_shuffle(seq: str, rng: np.random.Generator, block_size: int = 20) -> str:
    """Shuffle amino acids independently within non-overlapping 20-aa windows,
    preserving local compositional clustering (low-complexity blocks)."""
    chars = list(seq)
    for start in range(0, len(chars), block_size):
        sub = chars[start : start + block_size]
        rng.shuffle(sub)
        chars[start : start + block_size] = sub
    return "".join(chars)


def run_wsl_script(script_text: str, script_path: Path) -> str:
    script_path.write_text(script_text, encoding="utf-8", newline="\n")
    rel = script_path.relative_to(ROOT).as_posix()
    res = subprocess.run(
        ["wsl", "-d", WSL_DISTRO, "--", "bash", "-l", f"{WSL_ROOT}/{rel}"],
        capture_output=True, text=True, encoding="utf-8", check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(f"WSL failed ({res.returncode}):\n{res.stderr[:1200]}")
    return res.stdout


def bootstrap_p95_ci(vals: np.ndarray, rng: np.random.Generator, n_boot: int = 2000) -> tuple[float, float]:
    n = len(vals)
    idx = rng.integers(0, n, size=(n_boot, n))
    p95s = np.percentile(vals[idx], 95, axis=1)
    lo, hi = np.percentile(p95s, [2.5, 97.5])
    return float(lo), float(hi)


def main() -> int:
    t0 = time.time()
    meta = json.loads((RUNS / "audit_metadata.json").read_text(encoding="utf-8"))
    asm_map = meta["assemblies"]
    sp_order = [r["code"] for r in meta["enamel_intervals"]]

    # Read original real queries
    all_real = read_fasta(RUNS / "queries" / "all_real_queries.faa")
    comp_rows = {r["query_name"]: r for r in csv.DictReader((RUNS / "supplementary" / "Table_S3a_query_composition_seg.tsv").open(encoding="utf-8"), delimiter="\t")}
    orig_s3c = {(r["species_code"], r["seg_mode"], r["query_name"]): r for r in csv.DictReader((RUNS / "supplementary" / "Table_S3c_per_gene_seg_null_18sp.tsv").open(encoding="utf-8"), delimiter="\t")}

    # 1. Generate 1,000 global composition-preserving shuffles + 200 local-block shuffles
    qdir = OUT / "queries"
    qdir.mkdir(parents=True, exist_ok=True)
    rng_g = np.random.default_rng(SEED_GLOBAL)
    rng_b = np.random.default_rng(SEED_BLOCK)

    enamel_records: list[tuple[str, str]] = []
    for q in ENAMEL_QUERIES:
        seq = all_real[q]
        enamel_records.append((q, seq))
        arr = list(seq)
        for i in range(N_MAX):
            perm = arr.copy()
            rng_g.shuffle(perm)
            enamel_records.append((f"{q}~g{i:04d}", "".join(perm)))
        if q in BLOCK_NULL_QUERIES_ENAMEL:
            for i in range(N_BLOCK):
                enamel_records.append((f"{q}~b{i:03d}", block_shuffle(seq, rng_b, 20)))

    idp_records: list[tuple[str, str]] = []
    for q in IDP_QUERIES:
        seq = all_real[q]
        idp_records.append((q, seq))
        arr = list(seq)
        for i in range(N_MAX):
            perm = arr.copy()
            rng_g.shuffle(perm)
            idp_records.append((f"{q}~g{i:04d}", "".join(perm)))
        if q in BLOCK_NULL_QUERIES_IDP:
            for i in range(N_BLOCK):
                idp_records.append((f"{q}~b{i:03d}", block_shuffle(seq, rng_b, 20)))

    enamel_qfa = qdir / "enamel_1000shuf.faa"
    idp_qfa = qdir / "idp_1000shuf.faa"
    write_fasta(enamel_qfa, enamel_records)
    write_fasta(idp_qfa, idp_records)
    print(f"[1/4] Generated queries: enamel={len(enamel_records)} seqs, idp={len(idp_records)} seqs ({time.time()-t0:.1f}s)", flush=True)

    # 2. Run tblastn in WSL using GNU parallel across all 18 species (1 job per species -> 18 cores in parallel!)
    # Copy inputs to WSL local /tmp/enamel_c for blazing fast I/O, keep exact -subject per species and -evalue 10
    raw_out = OUT / "tblastn_raw"
    raw_out.mkdir(parents=True, exist_ok=True)
    # Only output qseqid, pident, length, bitscore, evalue -> huge speed & small disk, and collapse to max bitscore per qseqid in awk right in WSL!
    sh_text = f"""#!/usr/bin/env bash
set -euo pipefail
SRC="{WSL_ROOT}/runs/enamel_relics"
O="{WSL_ROOT}/runs/enamel_relics/v101/stats"
TMP="/tmp/enamel_phaseC_$$"
mkdir -p "$TMP/enam" "$TMP/idp" "$TMP/out"
cp "$SRC/intervals_enamel"/interval_*.fna "$TMP/enam/"
cp "$SRC/intervals_idp"/idp_*.fna "$TMP/idp/"
cp "$O/queries/enamel_1000shuf.faa" "$TMP/enamel_q.faa"
cp "$O/queries/idp_1000shuf.faa" "$TMP/idp_q.faa"

run_one() {{
  local dom="$1" seg="$2" sp="$3"
  local qfa subj outf
  if [[ "$dom" == "enamel" ]]; then
    qfa="$TMP/enamel_q.faa"
    subj="$TMP/enam/interval_${{sp}}.fna"
  else
    qfa="$TMP/idp_q.faa"
    subj="$TMP/idp/idp_${{sp}}.fna"
  fi
  outf="$TMP/out/${{dom}}_${{seg}}_${{sp}}.max.tsv"
  tblastn -query "$qfa" -subject "$subj" -evalue 10 -seg "$seg" -max_target_seqs 50 -num_threads 1 \
          -outfmt '6 qseqid pident length evalue bitscore' \
    | awk -F'\\t' '{{
        q=$1; b=$5+0;
        if (!(q in best) || b > best[q]) {{ best[q]=b; pid[q]=$2; len[q]=$3; ev[q]=$4; }}
      }} END {{
        for (q in best) printf "%s\\t%s\\t%s\\t%s\\t%.2f\\n", q, pid[q], len[q], ev[q], best[q];
      }}' > "$outf"
}}
export -f run_one
export TMP

SPS=({' '.join(sp_order)})
echo "[WSL] Starting 72 parallel tblastn jobs (18 sp x 2 domains x 2 SEG modes)..."
parallel -j 16 run_one {{1}} {{2}} {{3}} ::: enamel idp ::: yes no ::: "${{SPS[@]}}"
cp "$TMP/out"/*.max.tsv "$O/tblastn_raw/"
rm -rf "$TMP"
echo "[WSL] All 72 jobs completed."
"""
    t1 = time.time()
    out_log = run_wsl_script(sh_text, OUT / "run_phaseC_tblastn.sh")
    print(out_log.strip(), flush=True)
    print(f"[2/4] Completed all 72 x 1000-perm tblastn runs in {time.time()-t1:.1f}s", flush=True)

    # 3. Parse all max-bitscore files and compute full statistics
    best_bits: dict[tuple[str, str, str], dict[str, tuple[float, float, int, float]]] = {}
    for dom in ("enamel", "idp"):
        for seg in ("yes", "no"):
            for sp in sp_order:
                p = raw_out / f"{dom}_{seg}_{sp}.max.tsv"
                d: dict[str, tuple[float, float, int, float]] = {}
                for ln in p.read_text(encoding="utf-8").splitlines():
                    if not ln.strip():
                        continue
                    q, pid, alen, ev, bs = ln.split("\t")
                    d[q] = (float(bs), float(pid), int(alen), float(ev))
                best_bits[(dom, seg, sp)] = d

    rng_boot = np.random.default_rng(SEED_BOOT)
    per_gene_rows: list[dict[str, Any]] = []
    # Store full 1000-perm arrays for downstream index bootstrap
    perm_arrays: dict[tuple[str, str, str], np.ndarray] = {}  # (sp, seg, q) -> shape (1000,)
    obs_map: dict[tuple[str, str, str], float] = {}

    for sp in sp_order:
        sp_meta = asm_map[sp]
        for seg in ("yes", "no"):
            seg_mode = f"seg_{seg}"
            for dom, qlist, blist in [
                ("enamel", ENAMEL_QUERIES, BLOCK_NULL_QUERIES_ENAMEL),
                ("idp", IDP_QUERIES, BLOCK_NULL_QUERIES_IDP),
            ]:
                d = best_bits[(dom, seg, sp)]
                for q in qlist:
                    obs_t = d.get(q, (0.0, 0.0, 0, 10.0))
                    obs, top_pid, top_len, top_ev = obs_t
                    obs_map[(sp, seg_mode, q)] = obs
                    g_vals = np.array([d.get(f"{q}~g{i:04d}", (0.0, 0.0, 0, 10.0))[0] for i in range(N_MAX)], dtype=float)
                    perm_arrays[(sp, seg_mode, q)] = g_vals

                    orig = orig_s3c.get((sp, seg_mode, q), {})
                    orig_obs = float(orig.get("obs_bits", 0.0))
                    orig_p95_20 = float(orig.get("null_p95_bits", 0.0))
                    orig_exc_20 = float(orig.get("excess_over_p95", 0.0))

                    p95_by_n = {}
                    exc_by_n = {}
                    for n in N_LEVELS:
                        sub = g_vals[:n]
                        p95 = float(np.percentile(sub, 95))
                        p95_by_n[n] = p95
                        exc_by_n[n] = max(0.0, obs - p95)

                    ci_lo_20, ci_hi_20 = bootstrap_p95_ci(g_vals[:20], rng_boot, 1000)
                    ci_lo_1000, ci_hi_1000 = bootstrap_p95_ci(g_vals, rng_boot, 2000)
                    emp_p_1000 = (1.0 + float(np.sum(g_vals >= obs))) / (N_MAX + 1.0) if obs > 0 else 1.0

                    # Block-shuffle null (N=200) if run
                    b_p95 = ""
                    b_exc = ""
                    b_emp_p = ""
                    if q in blist:
                        b_vals = np.array([d.get(f"{q}~b{i:03d}", (0.0, 0.0, 0, 10.0))[0] for i in range(N_BLOCK)], dtype=float)
                        bp95 = float(np.percentile(b_vals, 95))
                        b_p95 = round(bp95, 2)
                        b_exc = round(max(0.0, obs - bp95), 2)
                        b_emp_p = round((1.0 + float(np.sum(b_vals >= obs))) / (N_BLOCK + 1.0) if obs > 0 else 1.0, 4)

                    per_gene_rows.append({
                        "species_code": sp,
                        "scientific_name": sp_meta["scientific_name"],
                        "clade": sp_meta["clade"],
                        "dentition": sp_meta["dentition"],
                        "domain": dom,
                        "seg_mode": seg_mode,
                        "query_name": q,
                        "obs_bits": round(obs, 2),
                        "orig_obs_bits": round(orig_obs, 2),
                        "obs_exact_match": abs(obs - orig_obs) < 0.05,
                        "top_hsp_pident": round(top_pid, 2),
                        "top_hsp_len_aa": top_len,
                        "top_hsp_evalue": f"{top_ev:.2e}",
                        "orig_p95_N20": round(orig_p95_20, 2),
                        "orig_exc_N20": round(orig_exc_20, 2),
                        "new_p95_N20": round(p95_by_n[20], 2),
                        "new_p95_N20_ci_lo": round(ci_lo_20, 2),
                        "new_p95_N20_ci_hi": round(ci_hi_20, 2),
                        "new_p95_N200": round(p95_by_n[200], 2),
                        "new_p95_N500": round(p95_by_n[500], 2),
                        "new_p95_N1000": round(p95_by_n[1000], 2),
                        "new_p95_N1000_ci_lo": round(ci_lo_1000, 2),
                        "new_p95_N1000_ci_hi": round(ci_hi_1000, 2),
                        "null_mean_N1000": round(float(np.mean(g_vals)), 2),
                        "null_max_N1000": round(float(np.max(g_vals)), 2),
                        "exc_N20": round(exc_by_n[20], 2),
                        "exc_N200": round(exc_by_n[200], 2),
                        "exc_N500": round(exc_by_n[500], 2),
                        "exc_N1000": round(exc_by_n[1000], 2),
                        "exc_N1000_ci_lo": round(max(0.0, obs - ci_hi_1000), 2),
                        "exc_N1000_ci_hi": round(max(0.0, obs - ci_lo_1000), 2),
                        "delta_p95_N1000_vs_orig20": round(p95_by_n[1000] - orig_p95_20, 2),
                        "delta_exc_N1000_vs_orig20": round(exc_by_n[1000] - orig_exc_20, 2),
                        "empirical_p_N1000": round(emp_p_1000, 4),
                        "block20_p95_N200": b_p95,
                        "block20_exc_N200": b_exc,
                        "block20_emp_p_N200": b_emp_p,
                    })

    # Write per-gene table
    pg_path = OUT / "perm_sensitivity_per_gene_v101.tsv"
    with pg_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(per_gene_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(per_gene_rows)

    # 4. Retention indices with N=1000 + bootstrap 95% CI propagation!
    idx_configs = [
        ("Amis_refseq_full", "Amis", "ENAM_Amis", "AMBN_Amis", "AMTN_Amis", "LPL_Amis"),
        ("Acar_equidistant", "Acar", "ENAM_Acar", "AMBN_Acar", "AMTN_Acar", "LPL_Acar"),
    ]
    pg_lookup = {(r["species_code"], r["seg_mode"], r["query_name"]): r for r in per_gene_rows}
    idx_rows = []
    for cfg_name, src, qe, qa, qt, ql in idx_configs:
        for seg_mode in ("seg_yes", "seg_no"):
            for sp in sp_order:
                sp_meta = asm_map[sp]
                re_ = pg_lookup[(sp, seg_mode, qe)]
                ra_ = pg_lookup[(sp, seg_mode, qa)]
                rt_ = pg_lookup[(sp, seg_mode, qt)]
                rl_ = pg_lookup[(sp, seg_mode, ql)]
                # Point estimate at N=1000
                num = re_["exc_N1000"] + ra_["exc_N1000"] + rt_["exc_N1000"]
                den = rl_["exc_N1000"]
                idx_1000 = num / den if den > 0 else float("nan")
                # Orig N=20 index
                num_o = re_["orig_exc_N20"] + ra_["orig_exc_N20"] + rt_["orig_exc_N20"]
                den_o = rl_["orig_exc_N20"]
                idx_orig20 = num_o / den_o if den_o > 0 else float("nan")
                # Bootstrap 95% CI propagating p95 uncertainty across all 4 queries!
                B = 2000
                boot_vals = []
                pe_arr = perm_arrays[(sp, seg_mode, qe)]
                pa_arr = perm_arrays[(sp, seg_mode, qa)]
                pt_arr = perm_arrays[(sp, seg_mode, qt)]
                pl_arr = perm_arrays[(sp, seg_mode, ql)]
                oe, oa, ot, ol = obs_map[(sp, seg_mode, qe)], obs_map[(sp, seg_mode, qa)], obs_map[(sp, seg_mode, qt)], obs_map[(sp, seg_mode, ql)]
                ie = np.maximum(0.0, oe - np.percentile(pe_arr[rng_boot.integers(0, N_MAX, size=(B, N_MAX))], 95, axis=1))
                ia = np.maximum(0.0, oa - np.percentile(pa_arr[rng_boot.integers(0, N_MAX, size=(B, N_MAX))], 95, axis=1))
                it = np.maximum(0.0, ot - np.percentile(pt_arr[rng_boot.integers(0, N_MAX, size=(B, N_MAX))], 95, axis=1))
                il = np.maximum(0.0, ol - np.percentile(pl_arr[rng_boot.integers(0, N_MAX, size=(B, N_MAX))], 95, axis=1))
                ratios = np.where(il > 0, (ie + ia + it) / il, np.nan)
                ci_lo, ci_hi = np.nanpercentile(ratios, [2.5, 97.5])
                idx_rows.append({
                    "query_set": cfg_name,
                    "seg_mode": seg_mode,
                    "species_code": sp,
                    "scientific_name": sp_meta["scientific_name"],
                    "clade": sp_meta["clade"],
                    "dentition": sp_meta["dentition"],
                    "ENAM_exc_1000": re_["exc_N1000"],
                    "AMBN_exc_1000": ra_["exc_N1000"],
                    "AMTN_exc_1000": rt_["exc_N1000"],
                    "LPL_exc_1000": rl_["exc_N1000"],
                    "orig_index_N20": round(idx_orig20, 3),
                    "index_N1000": round(idx_1000, 3),
                    "index_N1000_ci_lo": round(float(ci_lo), 3),
                    "index_N1000_ci_hi": round(float(ci_hi), 3),
                    "delta_vs_orig20": round(idx_1000 - idx_orig20, 3),
                })

    idx_path = OUT / "retention_indices_1000perm_v101.tsv"
    with idx_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(idx_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(idx_rows)

    # 5. Clade medians table for equidistant Anolis queries (Table 2 replacement / comparison!) + SPARCL1_Amis
    clade_rows = []
    for q in ["LPL_Acar", "IBSP_Acar", "SPP1_Acar", "ENAM_Acar", "AMBN_Acar", "AMTN_Acar", "SPARCL1_Amis", "IBSP_Amis", "SPP1_Amis"]:
        for seg_mode in ("seg_yes", "seg_no"):
            cstats = {}
            for cl in ("Crocodylia", "Testudines", "Aves", "Squamata"):
                sub = [r for r in per_gene_rows if r["query_name"] == q and r["seg_mode"] == seg_mode and r["clade"] == cl]
                e_o20 = [r["orig_exc_N20"] for r in sub]
                e_1k = [r["exc_N1000"] for r in sub]
                p95_o20 = [r["orig_p95_N20"] for r in sub]
                p95_1k = [r["new_p95_N1000"] for r in sub]
                pids = [r["top_hsp_pident"] for r in sub if r["obs_bits"] > 0]
                b_exc = [r["block20_exc_N200"] for r in sub if r["block20_exc_N200"] != ""]
                cstats[cl] = {
                    "orig20_med": round(float(np.median(e_o20)), 2),
                    "orig20_min": round(float(np.min(e_o20)), 2),
                    "orig20_max": round(float(np.max(e_o20)), 2),
                    "n1000_med": round(float(np.median(e_1k)), 2),
                    "n1000_min": round(float(np.min(e_1k)), 2),
                    "n1000_max": round(float(np.max(e_1k)), 2),
                    "p95_orig20_med": round(float(np.median(p95_o20)), 2),
                    "p95_1000_med": round(float(np.median(p95_1k)), 2),
                    "top_pident_med": round(float(np.median(pids)), 1) if pids else 0.0,
                    "block20_med": round(float(np.median(b_exc)), 2) if b_exc else None,
                    "n_emp_p_le_005": sum(1 for r in sub if r["empirical_p_N1000"] <= 0.05),
                    "n_emp_p_le_001": sum(1 for r in sub if r["empirical_p_N1000"] <= 0.01),
                    "n_total": len(sub),
                }
            aves_croc_orig = round(100.0 * cstats["Aves"]["orig20_med"] / cstats["Crocodylia"]["orig20_med"], 1) if cstats["Crocodylia"]["orig20_med"] > 0 else float("nan")
            aves_croc_1k = round(100.0 * cstats["Aves"]["n1000_med"] / cstats["Crocodylia"]["n1000_med"], 1) if cstats["Crocodylia"]["n1000_med"] > 0 else float("nan")
            clade_rows.append({
                "query_name": q,
                "seg_mode": seg_mode,
                "pct_seg": comp_rows[q]["pct_seg_masked"],
                "len_aa": comp_rows[q]["length_aa"],
                "croc_orig20": f"{cstats['Crocodylia']['orig20_med']} ({cstats['Crocodylia']['orig20_min']}-{cstats['Crocodylia']['orig20_max']})",
                "croc_1000": f"{cstats['Crocodylia']['n1000_med']} ({cstats['Crocodylia']['n1000_min']}-{cstats['Crocodylia']['n1000_max']})",
                "turt_1000": f"{cstats['Testudines']['n1000_med']} ({cstats['Testudines']['n1000_min']}-{cstats['Testudines']['n1000_max']})",
                "aves_orig20": f"{cstats['Aves']['orig20_med']} ({cstats['Aves']['orig20_min']}-{cstats['Aves']['orig20_max']})",
                "aves_1000": f"{cstats['Aves']['n1000_med']} ({cstats['Aves']['n1000_min']}-{cstats['Aves']['n1000_max']})",
                "aves_croc_pct_orig20": aves_croc_orig,
                "aves_croc_pct_1000": aves_croc_1k,
                "croc_pident_med": cstats["Crocodylia"]["top_pident_med"],
                "turt_pident_med": cstats["Testudines"]["top_pident_med"],
                "aves_pident_med": cstats["Aves"]["top_pident_med"],
                "aves_block20_med": cstats["Aves"]["block20_med"],
                "aves_p_le_001": f"{cstats['Aves']['n_emp_p_le_001']}/{cstats['Aves']['n_total']}",
                "turt_p_le_001": f"{cstats['Testudines']['n_emp_p_le_001']}/{cstats['Testudines']['n_total']}",
            })

    cm_path = OUT / "clade_medians_equidistant_1000perm_v101.tsv"
    with cm_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(clade_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(clade_rows)

    # 6. Regression/correlation of null_p95_N1000 and observed bitscores against length, SEG%, Pro+Gln%, entropy
    reg_summary: dict[str, Any] = {}
    for seg_mode in ("seg_yes", "seg_no"):
        # Aggregate per query (median across 18 species for null_p95_N1000)
        q_stats = []
        for q in ENAMEL_QUERIES + IDP_QUERIES:
            sub = [r for r in per_gene_rows if r["query_name"] == q and r["seg_mode"] == seg_mode]
            c = comp_rows[q]
            q_stats.append({
                "q": q,
                "len_aa": float(c["length_aa"]),
                "pct_seg": float(c["pct_seg_masked"]),
                "pct_pq": float(c["pct_Pro_Gln"]),
                "pct_scpp": float(c["pct_SCPP_disorder_PQSDEGN"]),
                "entropy": float(c["shannon_entropy_bits"]),
                "med_p95_1000": float(np.median([r["new_p95_N1000"] for r in sub])),
                "med_p95_ci_width_20": float(np.median([r["new_p95_N20_ci_hi"] - r["new_p95_N20_ci_lo"] for r in sub])),
                "med_p95_ci_width_1000": float(np.median([r["new_p95_N1000_ci_hi"] - r["new_p95_N1000_ci_lo"] for r in sub])),
            })
        y = np.array([x["med_p95_1000"] for x in q_stats])
        corrs = {}
        for feat in ("len_aa", "pct_seg", "pct_pq", "pct_scpp", "entropy"):
            xv = np.array([x[feat] for x in q_stats])
            r_p, p_p = stats.pearsonr(xv, y)
            r_s, p_s = stats.spearmanr(xv, y)
            corrs[feat] = {"pearson_r": round(float(r_p), 3), "pearson_p": round(float(p_p), 4), "spearman_rho": round(float(r_s), 3), "spearman_p": round(float(p_s), 4)}
        reg_summary[seg_mode] = {
            "n_queries": len(q_stats),
            "null_p95_correlations": corrs,
            "median_ci95_width_N20_bits": round(float(np.median([x["med_p95_ci_width_20"] for x in q_stats])), 2),
            "median_ci95_width_N1000_bits": round(float(np.median([x["med_p95_ci_width_1000"] for x in q_stats])), 2),
            "per_query": q_stats,
        }

    # Global stability summary across all 18 sp x 14 queries x 2 SEG modes = 504 cells
    deltas_p95 = [abs(r["delta_p95_N1000_vs_orig20"]) for r in per_gene_rows]
    deltas_exc = [abs(r["delta_exc_N1000_vs_orig20"]) for r in per_gene_rows]
    ci_w20 = [r["new_p95_N20_ci_hi"] - r["new_p95_N20_ci_lo"] for r in per_gene_rows]
    ci_w1000 = [r["new_p95_N1000_ci_hi"] - r["new_p95_N1000_ci_lo"] for r in per_gene_rows]
    # Check status flips: orig_exc_N20 > 0 vs exc_N1000 > 0 (and empirical_p <= 0.05)
    flips = [
        {
            "sp": r["species_code"], "seg": r["seg_mode"], "q": r["query_name"],
            "obs": r["obs_bits"], "orig_p95_20": r["orig_p95_N20"], "p95_1000": r["new_p95_N1000"],
            "orig_exc_20": r["orig_exc_N20"], "exc_1000": r["exc_N1000"], "emp_p": r["empirical_p_N1000"],
        }
        for r in per_gene_rows
        if (r["orig_exc_N20"] > 0) != (r["exc_N1000"] > 0)
    ]

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_sec": round(time.time() - t0, 1),
        "seeds": {"global": SEED_GLOBAL, "block20": SEED_BLOCK, "bootstrap": SEED_BOOT},
        "n_cells": len(per_gene_rows),
        "obs_exact_matches": f"{sum(r['obs_exact_match'] for r in per_gene_rows)}/{len(per_gene_rows)}",
        "p95_abs_diff_N1000_vs_orig20": {
            "median": round(float(np.median(deltas_p95)), 2),
            "p95": round(float(np.percentile(deltas_p95, 95)), 2),
            "max": round(float(np.max(deltas_p95)), 2),
        },
        "exc_abs_diff_N1000_vs_orig20": {
            "median": round(float(np.median(deltas_exc)), 2),
            "p95": round(float(np.percentile(deltas_exc, 95)), 2),
            "max": round(float(np.max(deltas_exc)), 2),
        },
        "ci95_width_N20_median_bits": round(float(np.median(ci_w20)), 2),
        "ci95_width_N1000_median_bits": round(float(np.median(ci_w1000)), 2),
        "zero_vs_nonzero_flips_count": len(flips),
        "zero_vs_nonzero_flips": flips,
        "clade_medians": clade_rows,
        "regression_summary": reg_summary,
    }
    sum_path = OUT / "phaseC_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    with PROV.open("a", encoding="utf-8") as f:
        f.write(json.dumps({
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "step": "phase_C_1000perm_sensitivity",
            "script": "scripts/enamel_relics/09_permutation_sensitivity.py",
            "seeds": summary["seeds"],
            "outputs": {
                p.name: sha256_file(p) for p in [pg_path, idx_path, cm_path, sum_path]
            },
        }) + "\n")

    print(f"\n[3/4] === PHASE C SUMMARY ({summary['elapsed_sec']}s) ===")
    print(f"  Observed bitscore exact reproduction: {summary['obs_exact_matches']}")
    print(f"  |p95(N=1000) - orig p95(N=20)|: median={summary['p95_abs_diff_N1000_vs_orig20']['median']} b, p95={summary['p95_abs_diff_N1000_vs_orig20']['p95']} b, max={summary['p95_abs_diff_N1000_vs_orig20']['max']} b")
    print(f"  Bootstrap 95% CI width on p95: N=20 median={summary['ci95_width_N20_median_bits']} b -> N=1000 median={summary['ci95_width_N1000_median_bits']} b")
    print(f"  Zero/nonzero excess flips (N=1000 vs orig N=20): {len(flips)}")
    for fl in flips[:20]:
        print(f"    FLIP: {fl}")
    print("\n  Clade Medians (Equidistant Anolis + Amis IDPs, N=1000 vs orig N=20):")
    for r in clade_rows:
        print(f"    {r['query_name']:13s} {r['seg_mode']:7s} | Croc={r['croc_1000']:22s} | Turt={r['turt_1000']:18s} | Aves={r['aves_1000']:18s} | Aves/Croc: {r['aves_croc_pct_orig20']}% -> {r['aves_croc_pct_1000']}% | pid C/T/A={r['croc_pident_med']}/{r['turt_pident_med']}/{r['aves_pident_med']}% | blk20={r['aves_block20_med']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
