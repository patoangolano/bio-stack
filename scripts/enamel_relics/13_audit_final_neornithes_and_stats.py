"""Phase A–C pre-submission audit script: Crown-Neornithes detectability limits,
cross-clade codon coverage breakdown, conditional ancestral-stop survival/detection
probabilities under empirical coverage masks, and exploratory p-value multiplicity audit.

Run: uv run python scripts/enamel_relics/13_audit_final_neornithes_and_stats.py
"""

from __future__ import annotations

import csv
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
SUPP = ROOT / "submission_biorxiv" / "supplementary"
OUT_DIR = RUNS / "v102"
OUT_JSON = OUT_DIR / "neornithes_detectability_audit.json"

BIRDS_LEGACY = ["Scam", "Arow", "Dnov", "Ggal", "Apla", "Shab", "Tgut"]
BIRDS_6 = ["Strcam", "Aptrow", "Dronov", "Galgal", "Anapla", "Strhab", "Taegut"]
PALAEO_IDX = [0, 1, 2]
NEOG_IDX = [3, 4, 5, 6]


def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip().upper()
    return seqs


def extract_ref_aligned_codons(aln: dict[str, str], ref_key: str = "Amis") -> tuple[int, dict[str, list[str]]]:
    ref_seq = aln[ref_key]
    ref_cols = [i for i, ch in enumerate(ref_seq) if ch != "-"]
    n_cod = len(ref_cols) // 3
    codons_by_sp: dict[str, list[str]] = {}
    for sp, s in aln.items():
        proj = "".join(s[i] for i in ref_cols[: n_cod * 3])
        codons_by_sp[sp] = [proj[3 * c : 3 * c + 3] for c in range(n_cod)]
    return n_cod, codons_by_sp


def is_valid_codon(c: str) -> bool:
    return len(c) == 3 and all(b in "ACGT" for b in c)


def k80_prob_matrix(t: float, kappa: float) -> np.ndarray:
    if t <= 0.0:
        return np.eye(4, dtype=np.float64)
    e1 = np.exp(-4.0 * t / (kappa + 2.0))
    e2 = np.exp(-2.0 * t * (kappa + 1.0) / (kappa + 2.0))
    p_same = 0.25 + 0.25 * e1 + 0.5 * e2
    p_ts = 0.25 + 0.25 * e1 - 0.5 * e2
    p_tv = 0.25 - 0.25 * e1
    # Order: A=0, C=1, G=2, T=3; transitions: A<->G (0<->2), C<->T (1<->3)
    m = np.full((4, 4), p_tv, dtype=np.float64)
    np.fill_diagonal(m, p_same)
    m[0, 2] = m[2, 0] = p_ts
    m[1, 3] = m[3, 1] = p_ts
    return m


def k80_cdf(t: float, kappa: float) -> np.ndarray:
    m = k80_prob_matrix(t, kappa)
    cdf = np.cumsum(m, axis=1)
    cdf[:, -1] = 1.0
    return cdf


def main() -> int:
    aln_dir = RUNS / "level4"
    exons = ["ENAM_exon8_aves", "ENAM_exon6_splice_aves", "AMBN_exon3_aves", "AMBN_exon4_aves", "AMBN_exon9_aves"]
    exon_cov_summary: dict[str, Any] = {}

    cov_dual_enam8 = None
    cov_mafft_enam8 = None

    for ex in exons:
        m_path = aln_dir / ex / "mafft.fna"
        u_path = aln_dir / ex / "muscle.fna"
        if not (m_path.exists() and u_path.exists()):
            continue
        aln_m = parse_fasta(m_path)
        aln_u = parse_fasta(u_path)
        n_cod, cod_m = extract_ref_aligned_codons(aln_m, "Amis")
        _, cod_u = extract_ref_aligned_codons(aln_u, "Amis")

        cov_m = np.zeros((7, n_cod), dtype=bool)
        cov_d = np.zeros((7, n_cod), dtype=bool)
        for bi, sp in enumerate(BIRDS_LEGACY):
            if sp in cod_m and sp in cod_u:
                for c in range(n_cod):
                    vm = is_valid_codon(cod_m[sp][c])
                    vu = is_valid_codon(cod_u[sp][c])
                    cov_m[bi, c] = vm
                    cov_d[bi, c] = vm and vu and (cod_m[sp][c] == cod_u[sp][c])

        if ex == "ENAM_exon8_aves":
            cov_dual_enam8 = cov_d
            cov_mafft_enam8 = cov_m

        pal_cnt_d = cov_d[PALAEO_IDX, :].sum(axis=0)
        neog_cnt_d = cov_d[NEOG_IDX, :].sum(axis=0)
        tot_cnt_d = cov_d.sum(axis=0)

        pal_cnt_m = cov_m[PALAEO_IDX, :].sum(axis=0)
        neog_cnt_m = cov_m[NEOG_IDX, :].sum(axis=0)
        tot_cnt_m = cov_m.sum(axis=0)

        exon_cov_summary[ex] = {
            "n_ref_codons": n_cod,
            "birds_present_in_msa": [BIRDS_6[bi] for bi, sp in enumerate(BIRDS_LEGACY) if sp in cod_m],
            "dual_aligner": {
                "all_7_birds_covered": int((tot_cnt_d == 7).sum()),
                "at_least_6_birds_covered": int((tot_cnt_d >= 6).sum()),
                "all_6_non_taegut_covered": int((cov_d[:6, :].sum(axis=0) == 6).sum()),
                "cross_clade_ge1_pal_ge1_neog": int(((pal_cnt_d >= 1) & (neog_cnt_d >= 1)).sum()),
                "cross_clade_ge2_pal_ge2_neog": int(((pal_cnt_d >= 2) & (neog_cnt_d >= 2)).sum()),
                "cross_clade_all3_pal_ge2_neog": int(((pal_cnt_d == 3) & (neog_cnt_d >= 2)).sum()),
                "cross_clade_all3_pal_ge3_neog": int(((pal_cnt_d == 3) & (neog_cnt_d >= 3)).sum()),
                "palaeo_only_covered_neog_0": int(((pal_cnt_d >= 1) & (neog_cnt_d == 0)).sum()),
                "neog_only_covered_palaeo_0": int(((pal_cnt_d == 0) & (neog_cnt_d >= 1)).sum()),
                "zero_birds_covered": int((tot_cnt_d == 0).sum()),
                "per_bird_covered_codons": {BIRDS_6[bi]: int(cov_d[bi, :].sum()) for bi in range(7)},
            },
            "mafft_only": {
                "all_7_birds_covered": int((tot_cnt_m == 7).sum()),
                "at_least_6_birds_covered": int((tot_cnt_m >= 6).sum()),
                "all_6_non_taegut_covered": int((cov_m[:6, :].sum(axis=0) == 6).sum()),
                "cross_clade_ge1_pal_ge1_neog": int(((pal_cnt_m >= 1) & (neog_cnt_m >= 1)).sum()),
                "cross_clade_ge2_pal_ge2_neog": int(((pal_cnt_m >= 2) & (neog_cnt_m >= 2)).sum()),
                "cross_clade_all3_pal_ge3_neog": int(((pal_cnt_m == 3) & (neog_cnt_m >= 3)).sum()),
            },
        }

    # Phase A2: Detection probability of a single genuine crown-Neornithes ancestral PTC
    # (present at the crown Neornithes root node prior to the Palaeognathae-Neognathae split)
    # evolving neutrally along empirical T1 branch lengths and observed through:
    # (i) unmasked 878 codons, (ii) empirical_mafft_mask, (iii) empirical_dual_aligner_mask
    assert cov_dual_enam8 is not None and cov_mafft_enam8 is not None
    # Load T1 empirical branch lengths from 10_audit_mutations_and_asr.py / treefile
    # From ENAM_exon8_aves T1 treefile:
    bl_t1 = {
        "t_pal": 0.0954433262,
        "t_subpal": 0.0252799710,
        "t_neo": 0.0566323799,
        "t_gal": 0.0680550918,
        "t_psi": 0.0706716885,
        "tips": {
            "Scam": 0.0587596682,
            "Arow": 0.0710557964,
            "Dnov": 0.0439344061,
            "Ggal": 0.1456338642,
            "Apla": 0.1235224258,
            "Shab": 0.0908253438,
            "Tgut": 0.2733774383,
        },
    }

    rng = np.random.default_rng(20261124)
    b2i = {"A": 0, "C": 1, "G": 2, "T": 3}
    stop_triplets = ["TAA", "TAG", "TGA"]
    # Simulate 878 codons x 100 replicates (= 87,800 ancestral stop trials) per stop codon & parameter regime
    n_cod = 878
    n_trials_per_cod = 100
    n_tot = n_cod * n_trials_per_cod

    def step_seq(seq: np.ndarray, cdf: np.ndarray) -> np.ndarray:
        u = rng.random(len(seq))
        c0, c1, c2 = cdf[seq, 0], cdf[seq, 1], cdf[seq, 2]
        return np.where(u < c0, 0, np.where(u < c1, 1, np.where(u < c2, 2, 3))).astype(np.int8)

    def is_stop_arr(seq: np.ndarray) -> np.ndarray:
        c0, c1, c2 = seq[0::3], seq[1::3], seq[2::3]
        return (c0 == 3) & (((c1 == 0) & ((c2 == 0) | (c2 == 2))) | ((c1 == 2) & (c2 == 0)))

    anc_stop_power: list[dict[str, Any]] = []
    for mask_name, cov_mat in [
        ("unmasked_full_878", np.ones((7, n_cod), dtype=bool)),
        ("empirical_mafft_mask", cov_mafft_enam8),
        ("empirical_dual_aligner_mask", cov_dual_enam8),
    ]:
        cov_tiled = np.tile(cov_mat, (1, n_trials_per_cod))
        for scale in [0.0, 0.65, 1.0, 1.25]:
            kappa = 4.0
            cdf_pal = k80_cdf(bl_t1["t_pal"] * scale, kappa)
            cdf_sub = k80_cdf(bl_t1["t_subpal"] * scale, kappa)
            cdf_neo = k80_cdf(bl_t1["t_neo"] * scale, kappa)
            cdf_gal = k80_cdf(bl_t1["t_gal"] * scale, kappa)
            cdf_psi = k80_cdf(bl_t1["t_psi"] * scale, kappa)
            cdf_tips = {sp: k80_cdf(bl_t1["tips"][sp] * scale, kappa) for sp in BIRDS_LEGACY}

            # Mix equally across TAA, TAG, TGA at crown Neornithes root
            root_nt = "".join(stop_triplets[i % 3] for i in range(n_tot))
            root_arr = np.array([b2i[c] for c in root_nt], dtype=np.int8)

            pal = step_seq(root_arr, cdf_pal)
            sub = step_seq(pal, cdf_sub)
            neo = step_seq(root_arr, cdf_neo)
            gal = step_seq(neo, cdf_gal)
            psi = step_seq(neo, cdf_psi)

            s_scam = is_stop_arr(step_seq(pal, cdf_tips["Scam"]))
            s_arow = is_stop_arr(step_seq(sub, cdf_tips["Arow"]))
            s_dnov = is_stop_arr(step_seq(sub, cdf_tips["Dnov"]))
            s_ggal = is_stop_arr(step_seq(gal, cdf_tips["Ggal"]))
            s_apla = is_stop_arr(step_seq(gal, cdf_tips["Apla"]))
            s_shab = is_stop_arr(step_seq(psi, cdf_tips["Shab"]))
            s_tgut = is_stop_arr(step_seq(psi, cdf_tips["Tgut"]))

            stops_mat = np.vstack([s_scam, s_arow, s_dnov, s_ggal, s_apla, s_shab, s_tgut])
            obs_stop = stops_mat & cov_tiled
            obs_sense = (~stops_mat) & cov_tiled

            n_pal_cov = cov_tiled[PALAEO_IDX, :].sum(axis=0)
            n_neog_cov = cov_tiled[NEOG_IDX, :].sum(axis=0)
            n_pal_stop = obs_stop[PALAEO_IDX, :].sum(axis=0)
            n_neog_stop = obs_stop[NEOG_IDX, :].sum(axis=0)
            n_tot_stop = obs_stop.sum(axis=0)
            n_tot_sense = obs_sense.sum(axis=0)

            # Criteria for detecting a crown-Neornithes synapomorphy:
            # 1. Strict 7/7 universal stop: all 7 birds covered and stop
            strict_7of7 = (n_tot_stop == 7)
            # 2. Strict 6/6 non-passerine stop (3/3 Palaeo + 3/3 non-Taegut Neog all stop, 0 sense anywhere)
            strict_6nonpass = (obs_stop[:6, :].sum(axis=0) == 6) & (n_tot_sense == 0)
            # 3. Multi-taxon uncontradicted cross-clade stop (>=2 Palaeo stop & >=2 Neog stop & 0 sense anywhere)
            uncontradicted_ge2_each = (n_pal_stop >= 2) & (n_neog_stop >= 2) & (n_tot_sense == 0)
            # 4. Relaxed uncontradicted cross-clade stop (lut_crown_neorn in Table S8: >=1 Palaeo stop & >=1 Neog stop & 0 sense)
            uncontradicted_ge1_each = (n_pal_stop >= 1) & (n_neog_stop >= 1) & (n_tot_sense == 0)
            # 5. Reversal-corrupted cross-clade polymorphism (>=1 Palaeo stop & >=1 Neog stop, BUT >=1 covered bird reverted to sense)
            reversal_corrupted_cross_clade = (n_pal_stop >= 1) & (n_neog_stop >= 1) & (n_tot_sense >= 1)
            # 6. Clade-truncated by missing coverage (all covered birds are stop, >=2 carriers, but 0 covered in Palaeo or 0 covered in Neog)
            clade_truncated_by_missing_cov = (n_tot_stop >= 2) & (n_tot_sense == 0) & ((n_pal_cov == 0) | (n_neog_cov == 0))

            anc_stop_power.append({
                "mask_mode": mask_name,
                "branch_scale": scale,
                "note": "scale=0.0 isolates pure coverage-mask effect (0 substitutions/reversals)" if scale == 0.0 else "neutral K80 stop-to-sense reversal + coverage mask",
                "p_strict_7of7_universal": round(float(np.mean(strict_7of7)), 4),
                "p_strict_6nonpass_uncontradicted": round(float(np.mean(strict_6nonpass)), 4),
                "p_uncontradicted_ge2_pal_ge2_neog": round(float(np.mean(uncontradicted_ge2_each)), 4),
                "p_uncontradicted_ge1_pal_ge1_neog_TableS8_def": round(float(np.mean(uncontradicted_ge1_each)), 4),
                "p_reversal_corrupted_cross_clade_appears_incongruent": round(float(np.mean(reversal_corrupted_cross_clade)), 4),
                "p_clade_truncated_appears_subclade_only": round(float(np.mean(clade_truncated_by_missing_cov)), 4),
            })

    # Phase C1: Audit exploratory empirical p-values in Table_S3c (504 cells)
    s3c_rows = list(csv.DictReader((SUPP / "Table_S3c_per_gene_seg_null_18sp.tsv").open(encoding="utf-8"), delimiter="\t"))
    marginal_cells = [
        {
            "species": r["species_code_6letter"],
            "dentition": r["dentition"],
            "query": r["query_name"],
            "seg_mode": r["seg_mode"],
            "obs_bits": float(r["obs_bits"]),
            "p95_1000": float(r["new_p95_N1000"]),
            "excess_1000": float(r["exc_N1000"]),
            "emp_p_1000": float(r["empirical_p_N1000"]),
        }
        for r in s3c_rows
        if 0.005 < float(r["empirical_p_N1000"]) < 0.10
    ]

    out_data = {
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "script": "scripts/enamel_relics/13_audit_final_neornithes_and_stats.py",
        "seed": 20261124,
        "exon_coverage_breakdown": exon_cov_summary,
        "ancestral_neornithes_stop_detection_power": anc_stop_power,
        "table_s3c_marginal_p_cells": marginal_cells,
    }
    OUT_JSON.write_text(json.dumps(out_data, indent=2), encoding="utf-8")
    print(json.dumps(out_data, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
