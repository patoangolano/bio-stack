"""Consolidate ALL audited v1.0.1 supplementary tables (Table_S1 through Table_S7 + species_ids.tsv)
with unambiguous 6-letter species codes, exact 1-based inclusive & 0-based BED indexing,
Phase B1 Exon 6 reciprocal audit, Phase B2 Table_S7 native model audit, Phase C 1,000-permutation
statistics + 20-aa block null + bootstrap 95% CIs, Phase D IQ-TREE 3 ML ASR + flanked-MSA splice borders,
and clean English Table_S6 with embedded numeric odaph_gap_bp (removing all legacy Portuguese strings).

Writes to runs/enamel_relics/supplementary/ and mirrors identically to submission_biorxiv/supplementary/,
regenerating MANIFEST_SHA256.tsv and appending provenance.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "runs" / "enamel_relics"
V101 = RUNS / "v101"
SUPP_RUNS = RUNS / "supplementary"
SUPP_SUB = ROOT / "submission_biorxiv" / "supplementary"
LEGACY = ROOT / "runs" / "atavismo"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    sp_map_rows = read_tsv(V101 / "species_ids.tsv")
    old2new = {r["legacy_4letter_code"]: r["unambiguous_6letter_code"] for r in sp_map_rows}

    # --- Table S1: add unambiguous 6-letter code + Phase B2 exact Gnomon audit notes ---
    s1_in = read_tsv(SUPP_RUNS / "Table_S1_assemblies_and_annotations.tsv")
    s7_in = read_tsv(V101 / "b2" / "Table_S7_native_model_audit_v101.tsv")
    s7_by_code = {r["locus"].split("_")[0]: r for r in s7_in}
    s1_out = []
    for r in s1_in:
        code = r.get("legacy_code") or r["species_code"]
        code6 = old2new.get(code, code)
        note = r["audit_note"]
        if code in s7_by_code:
            b2r = s7_by_code[code]
            note = (
                f"Gnomon LOW QUALITY PROTEIN ({b2r['transcript']}/{b2r['protein']}): "
                f"{b2r['n_gnomon_fs_compensations']} frameshift + {b2r['n_gnomon_ptc_compensations']} PTC compensated in native genome "
                f"[{b2r['revised_classification']}; see Table S7]"
            )
        row = {"species_code_6letter": code6, "legacy_code": code}
        for k, v in r.items():
            if k in ("species_code", "species_code_6letter", "legacy_code"):
                continue
            row[k] = v
        row["audit_note"] = note
        s1_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S1_assemblies_and_annotations.tsv", s1_out)

    # --- Table S2a: add 6-letter code + explicit 1-based inclusive / 0-based BED intergenic gap ---
    s2a_in = read_tsv(RUNS / "enamel_intervals.tsv")
    s2a_out = []
    for r in s2a_in:
        code = r["code"]
        code6 = old2new[code]
        strict_gap = ""
        inter_1b = ""
        inter_bed = ""
        if r["anchor_mode"] == "pair_NRG1_LPL":
            ns, ne = map(int, r["nrg1_coords"].split("(")[0].split("-"))
            ls, le = map(int, r["lpl_coords"].split("(")[0].split("-"))
            if ne < ls:
                lo_g, hi_g = ne + 1, ls - 1
            else:
                lo_g, hi_g = le + 1, ns - 1
            strict_gap = str(hi_g - lo_g + 1)
            inter_1b = f"{r['seq_acc']}:{lo_g}-{hi_g}"
            inter_bed = f"{r['seq_acc']}:{lo_g - 1}-{hi_g}"
        row = {"species_code_6letter": code6, "legacy_code": code}
        for k, v in r.items():
            if k == "code":
                continue
            row[k] = v
        row["strict_intergenic_bp"] = strict_gap
        row["intergenic_1based_inclusive"] = inter_1b
        row["intergenic_bed_0based"] = inter_bed
        s2a_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S2a_synteny_intervals_18sp.tsv", s2a_out)

    # --- Table S2b: add 6-letter code + explicit Phase B1 Ggal Psi-ENAM Exon 6 reciprocal audit row ---
    s2b_in = read_tsv(RUNS / "tables" / "whole_proteome_reciprocal.tsv")
    s2b_out = []
    for r in s2b_in:
        code = r["species_code"]
        row = {"species_code_6letter": old2new[code], "legacy_code": code}
        for k, v in r.items():
            if k == "species_code":
                continue
            row[k] = v
        s2b_out.append(row)
    # Append Phase B1 audited Chicken Psi-ENAM Exon 6 row so Table S2b explicitly records its negative reciprocal BLASTx!
    b1 = json.loads((V101 / "phaseAB" / "phase_A_B1_B3_audit.json").read_text(encoding="utf-8"))["phase_b1"]
    top_bx = b1["blastx_std_segyes_top3"]["Ggal_ENAM_ex6_blastn149"][0]
    s2b_out.append({
        "species_code_6letter": "Galgal",
        "legacy_code": "Ggal",
        "candidate_id": "Ggal__ENAM_Exon6_B1_audit",
        "query_name": "ENAM_Exon6_BLASTn_vs_Arow_and_Amis",
        "seq_acc": "NC_052572.1",
        "seq_name": "Z",
        "strand": "-",
        "n_clustered_hsps": "1",
        "total_hsps_ge28": "1",
        "top_tblastn_pid": "78.000(BLASTn_Arow)/72.593(BLASTn_Amis)",
        "top_tblastn_len": "150(nt_Arow)/135(nt_Amis)",
        "top_tblastn_qspan": "Exon6_codons_168-183(47bp_core_in_149bp_span)",
        "top_tblastn_evalue": "2.71e-27(Arow)/5.25e-11(Amis)",
        "top_tblastn_bits": "115.0(Arow)/60.8(Amis)",
        "abs_start": "55160529",
        "abs_end": "55160677",
        "span_bp": "149",
        "extracted_start": "55160570",
        "extracted_end": "55160616",
        "top_hsp_start": "55160570",
        "top_hsp_end": "55160616",
        "proteome_size": "48040",
        "reciprocal_confirmed_span": "False",
        "top_hsp_reciprocal_confirmed": "False",
        "top_proteome_acc": top_bx["sseqid"],
        "top_proteome_title": f"{top_bx['stitle']} [FAILED_RECIPROCAL: 47bp exon with TAA@cod4 + 1-nt del; Amis 48bp ctrl hits XP_019346377.1 37.7b E=3.6e-5]",
        "top_proteome_pid": str(top_bx["pident"]),
        "top_proteome_len": str(top_bx["length"]),
        "top_proteome_evalue": str(top_bx["evalue"]),
        "top_proteome_bits": str(top_bx["bitscore"]),
        "runnerup_nonortholog_acc": "XP_019340150.1(core47)",
        "runnerup_nonortholog_title": "dolichyl-diphosphooligosaccharide--protein glycosyltransferase subunit 2 (E=4.8, 23.1b)",
        "runnerup_nonortholog_evalue": "4.8",
        "runnerup_nonortholog_bits": "23.1",
        "top_hsp_proteome_acc": "XP_019340150.1",
        "top_hsp_proteome_evalue": "4.8",
        "top_hsp_proteome_bits": "23.1",
        "top_hsp_runnerup_acc": "NONE",
        "top_hsp_runnerup_bits": "",
    })
    write_tsv(SUPP_RUNS / "Table_S2b_whole_proteome_reciprocal_blastx.tsv", s2b_out)

    # --- Table S3a (unchanged composition, add source_species_6letter) ---
    s3a_in = read_tsv(RUNS / "tables" / "query_composition_seg.tsv")
    s3a_out = []
    for r in s3a_in:
        r2 = dict(r)
        r2["source_species_6letter"] = old2new.get(r["source_species"], r["source_species"])
        s3a_out.append(r2)
    write_tsv(SUPP_RUNS / "Table_S3a_query_composition_seg.tsv", s3a_out)

    # --- Table S3b & S3c: upgrade directly from Phase C 1,000-permutation + 20-aa block null outputs! ---
    s3b_in = read_tsv(V101 / "stats" / "retention_indices_1000perm_v101.tsv")
    s3b_out = []
    for r in s3b_in:
        code = r["species_code"]
        row = {"species_code_6letter": old2new[code], "legacy_code": code}
        for k, v in r.items():
            if k != "species_code":
                row[k] = v
        s3b_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S3b_retention_indices_18sp.tsv", s3b_out)

    s3c_in = read_tsv(V101 / "stats" / "perm_sensitivity_per_gene_v101.tsv")
    s3c_out = []
    for r in s3c_in:
        code = r["species_code"]
        row = {"species_code_6letter": old2new[code], "legacy_code": code}
        for k, v in r.items():
            if k != "species_code":
                row[k] = v
        s3c_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S3c_per_gene_seg_null_18sp.tsv", s3c_out)

    # --- Table S4: add 6-letter code + keep 25-kb QC note ---
    s4_in = read_tsv(RUNS / "tables" / "miniprot_level3_splice_audit.tsv")
    s4_out = []
    for r in s4_in:
        code = r["species_code"]
        row = {"species_code_6letter": old2new[code], "legacy_code": code}
        for k, v in r.items():
            if k != "species_code":
                row[k] = v
        if r["query_gene"] == "ENAM_Amis" and code == "Cpic":
            row["qc_note"] = "Spurious 58.4-kb upstream exon jump (NC_088796.1:129180576-129180957, reciprocal BLASTx = protein SON, E=0.010) filtered by 25-kb max-intron rule"
        else:
            row["qc_note"] = ""
        s4_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S4_miniprot_level3_and_splice_audit.tsv", s4_out)

    # --- Table S5a & S5b: upgrade directly from v1.0.2 / Phase D IQ-TREE 3 ML ASR + flanked-MSA splice borders ---
    v102_s5a = RUNS / "v102" / "Table_S5a_level4_shared_mutations_v102.tsv"
    s5a_src = v102_s5a if v102_s5a.exists() else (V101 / "phaseD" / "Table_S5a_level4_shared_mutations_v101.tsv")
    write_tsv(SUPP_RUNS / "Table_S5a_level4_shared_mutations.tsv", read_tsv(s5a_src))
    write_tsv(SUPP_RUNS / "Table_S5b_level4_splice_borders.tsv", read_tsv(V101 / "phaseD" / "Table_S5b_level4_splice_borders_v101.tsv"))

    # --- Table S6: Fix Portuguese leak + embed numeric odaph_gap_bp from s32_vaos.json + 6-letter codes! ---
    odaph_gaps = json.loads((LEGACY / "odam_odaph" / "s32_vaos.json").read_text(encoding="utf-8"))
    s6_in = read_tsv(SUPP_RUNS / "Table_S6_auxiliary_odaph_multiz_recomb_psiamel.tsv")
    s6_out = []
    for r in s6_in:
        code = r.get("legacy_code") or r["species_code"]
        code6 = old2new.get(code, code)
        ctrl = r["control_detail"].replace(
            "ODAM nao localizado em none especie COM dente; sensibilidade no ODAM do Xenopus:",
            "ODAM not recovered in any toothed sauropsid; Xenopus positive-control sensitivity:",
        ).replace(
            "ODAM nao localizado em nenhuma especie COM dente; sensibilidade no ODAM do Xenopus:",
            "ODAM not recovered in any toothed sauropsid; Xenopus positive-control sensitivity:",
        )
        gap_bp = ""
        if r["module"] == "ODAPH_ODAM" and r["gene"] == "ODAPH" and code in odaph_gaps:
            gap_bp = str(odaph_gaps[code]["bp"])
        row = {
            "module": r["module"],
            "gene": r["gene"],
            "species_code_6letter": code6,
            "legacy_code": code,
            "scientific_name": r["scientific_name"],
            "clade": r["clade"],
            "dentition": r["dentition"],
            "odaph_gap_bp": gap_bp,
            "metric_summary": r["metric_summary"],
            "evidence_level": r["evidence_level"],
            "coordinates": r["coordinates"],
            "control_detail": ctrl,
        }
        s6_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S6_auxiliary_odaph_multiz_recomb_psiamel.tsv", s6_out)

    # --- Table S7: NEW Phase B2 native RefSeq PROTEIN_CODING model audit! ---
    s7_out = []
    for r in s7_in:
        code = r["locus"].split("_")[0]
        row = {"species_code_6letter": old2new[code], "legacy_code": code}
        row.update(r)
        s7_out.append(row)
    write_tsv(SUPP_RUNS / "Table_S7_native_model_audit_v101.tsv", s7_out)

    # --- Also include species_ids.tsv as Table_S0 / reference key (strict LF via write_tsv) ---
    write_tsv(SUPP_RUNS / "Table_S0_species_identifiers.tsv", read_tsv(V101 / "species_ids.tsv"))

    # --- Also include v1.0.2 Table_S8 (K80 simulation sensitivity) & Table_S9 (T1/T2/T3 topology sensitivity) ---
    v102_s8 = RUNS / "v102" / "Table_S8_k80_simulation_sensitivity_v102.tsv"
    v102_s9 = RUNS / "v102" / "Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv"
    if v102_s8.exists():
        write_tsv(SUPP_RUNS / "Table_S8_k80_simulation_sensitivity_v102.tsv", read_tsv(v102_s8))
    if v102_s9.exists():
        write_tsv(SUPP_RUNS / "Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv", read_tsv(v102_s9))

    # --- Regenerate MANIFEST_SHA256.tsv and mirror identically to submission_biorxiv/supplementary/ ---
    manifest = []
    for p in sorted(SUPP_RUNS.glob("*.tsv")):
        if p.name == "MANIFEST_SHA256.tsv":
            continue
        manifest.append({
            "file": p.name,
            "bytes": p.stat().st_size,
            "sha256": sha256_file(p),
            "generated_utc": datetime.now(UTC).isoformat(),
        })
    write_tsv(SUPP_RUNS / "MANIFEST_SHA256.tsv", manifest)

    # Mirror cleanly to submission_biorxiv/supplementary/
    SUPP_SUB.mkdir(parents=True, exist_ok=True)
    for old in SUPP_SUB.glob("*.tsv"):
        old.unlink()
    for p in sorted(SUPP_RUNS.glob("*.tsv")):
        shutil.copy2(p, SUPP_SUB / p.name)

    print(f"Consolidated and mirrored {len(manifest)} v1.0.1 supplementary TSVs (+ MANIFEST_SHA256.tsv):")
    for m in manifest:
        print(f"  {m['file']:55s} {m['bytes']:>8,} B  {m['sha256'][:16]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
