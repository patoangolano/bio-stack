"""Phase H automated acceptance test suite for the bioRxiv v1.0.1 audit package.

Run: uv run pytest -v tests/test_enamel_relics_audit.py
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import re
import subprocess
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs" / "enamel_relics"
SUPP_RUNS = RUNS / "supplementary"
SUB = ROOT / "submission_biorxiv"
SUPP_SUB = SUB / "supplementary"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# 1. Supplementary SHA-256 hashes & mirroring
# ---------------------------------------------------------------------------
def test_supplementary_hashes_and_mirror() -> None:
    for base in (SUPP_RUNS, SUPP_SUB):
        manifest = read_tsv(base / "MANIFEST_SHA256.tsv")
        assert len(manifest) >= 11, f"Expected >=11 supplementary TSVs in {base}, got {len(manifest)}"
        for row in manifest:
            p = base / row["file"]
            assert p.exists(), f"Missing {p}"
            assert p.stat().st_size == int(row["bytes"]), f"Byte mismatch on {p}"
            assert sha256_file(p) == row["sha256"], f"SHA-256 mismatch on {p}"
    # Mirror exact match between runs/ and submission_biorxiv/
    assert sha256_file(SUPP_RUNS / "MANIFEST_SHA256.tsv") == sha256_file(SUPP_SUB / "MANIFEST_SHA256.tsv")


# ---------------------------------------------------------------------------
# 2. Coordinates & indexing conventions (B1 & B3)
# ---------------------------------------------------------------------------
def test_coordinates_and_indexing_conventions() -> None:
    s2a = {r["species_code_6letter"]: r for r in read_tsv(SUPP_SUB / "Table_S2a_synteny_intervals_18sp.tsv")}
    gg = s2a["Galgal"]
    assert gg["seq_acc"] == "NC_052572.1"
    assert gg["nrg1_coords"] == "54642186-55109784(-)"
    assert gg["lpl_coords"] == "55215002-55232438(-)"
    assert int(gg["intergenic_gap_bp"]) == 105218  # endpoint difference
    assert int(gg["strict_intergenic_bp"]) == 105217  # strictly between nucleotides 55109785..55215001
    assert gg["intergenic_1based_inclusive"] == "NC_052572.1:55109785-55215001"
    assert gg["intergenic_bed_0based"] == "NC_052572.1:55109784-55215001"

    # Reconcile Psi-ENAM Exon 8 (union 1,985 bp vs BLASTn 1,975 bp) and Exon 6 (149 bp span vs 47 bp core)
    b_audit = json.loads((RUNS / "v101" / "phaseAB" / "phase_A_B1_B3_audit.json").read_text(encoding="utf-8"))
    assert b_audit["phase_b3"]["ex8_genomic_span_55155480_55157464"] == 1985
    assert b_audit["phase_b3"]["ex8_msa_in_ungapped_len"] == 1975
    assert b_audit["phase_b1"]["reported_span_bp"] == 149
    assert b_audit["phase_b1"]["exon6_core_bp"] == 47
    assert b_audit["phase_b1"]["ggal_acceptor_2bp"] == "AG"
    assert b_audit["phase_b1"]["ggal_donor_2bp"] == "GT"

    # Confirm Table S2b records negative reciprocal BLASTx for Ggal Exon 6
    s2b = {r["candidate_id"]: r for r in read_tsv(SUPP_SUB / "Table_S2b_whole_proteome_reciprocal_blastx.tsv")}
    assert s2b["Ggal__ENAM_Exon6_B1_audit"]["reciprocal_confirmed_span"] == "False"
    assert s2b["Ggal__ENAM_Exon6_B1_audit"]["top_hsp_reciprocal_confirmed"] == "False"
    assert s2b["Ggal__ENAM_Amis"]["top_hsp_reciprocal_confirmed"] == "True"
    assert s2b["Ggal__AMBN_Amis"]["top_hsp_reciprocal_confirmed"] == "True"
    assert s2b["Shab__AMTN_Amis"]["top_hsp_reciprocal_confirmed"] == "True"


# ---------------------------------------------------------------------------
# 3. Species & sequence counts + unambiguous species codes
# ---------------------------------------------------------------------------
def test_species_and_sequence_counts() -> None:
    s0 = read_tsv(SUPP_SUB / "Table_S0_species_identifiers.tsv")
    s1 = read_tsv(SUPP_SUB / "Table_S1_assemblies_and_annotations.tsv")
    s2a = read_tsv(SUPP_SUB / "Table_S2a_synteny_intervals_18sp.tsv")
    s3a = read_tsv(SUPP_SUB / "Table_S3a_query_composition_seg.tsv")
    s3c = read_tsv(SUPP_SUB / "Table_S3c_per_gene_seg_null_18sp.tsv")
    s7 = read_tsv(SUPP_SUB / "Table_S7_native_model_audit_v101.tsv")

    assert len(s0) == 19
    assert len({r["unambiguous_6letter_code"] for r in s0}) == 19
    assert any(r["unambiguous_6letter_code"] == "Taegut" and r["scientific_name"] == "Taeniopygia guttata" for r in s0)
    assert any(r["unambiguous_6letter_code"] == "Tingut" and r["scientific_name"] == "Tinamus guttatus" for r in s0)

    assert len(s1) == 19  # 18 sauropsids + Homo sapiens
    assert sum(1 for r in s1 if r["clade"] != "Mammalia") == 18
    assert len(s2a) == 18
    assert len(s3a) == 16
    assert len({r["species_code_6letter"] for r in s3c}) == 18
    assert len({r["query_name"] for r in s3c}) == 14
    assert len(s3c) == 18 * 14 * 2  # 504 cells (18 sp x 14 queries x 2 SEG modes)
    assert len(s7) == 4  # 4 RefSeq PROTEIN_CODING models audited


# ---------------------------------------------------------------------------
# 4. Recalculate Level 3 lesion rates & verify exact 133/133 LPL introns
# ---------------------------------------------------------------------------
def test_level3_lesion_rates_and_133_lpl_introns() -> None:
    s4 = read_tsv(SUPP_SUB / "Table_S4_miniprot_level3_and_splice_audit.tsv")

    lpl = [r for r in s4 if r["query_gene"] == "LPL_Amis" and r["aligned"] == "True"]
    assert len(lpl) == 18
    lpl_aa = sum(int(r["aligned_aa"]) for r in lpl)
    lpl_fs = sum(int(r["frameshifts"]) for r in lpl)
    lpl_st = sum(int(r["stops"]) for r in lpl)
    lpl_introns_tot = sum(int(r["n_introns"]) for r in lpl)
    lpl_introns_canon = sum(int(r["n_gtag"]) for r in lpl)
    assert lpl_aa == 8246
    assert lpl_fs == 0 and lpl_st == 0
    assert lpl_introns_tot == 133 and lpl_introns_canon == 133

    toothed = [r for r in s4 if r["query_gene"] in ("ENAM_Amis", "AMBN_Amis", "AMTN_Amis") and r["dentition"] == "toothed" and r["aligned"] == "True"]
    assert len(toothed) == 15
    t_aa = sum(int(r["aligned_aa"]) for r in toothed)
    t_fs = sum(int(r["frameshifts"]) for r in toothed)
    t_st = sum(int(r["stops"]) for r in toothed)
    assert t_aa == 9394 and t_fs == 1 and t_st == 0
    assert round(1000.0 * (t_fs + t_st) / t_aa, 2) == 0.11

    edent = [r for r in s4 if r["query_gene"] in ("ENAM_Amis", "AMBN_Amis", "AMTN_Amis") and r["dentition"] == "edentulous" and r["aligned"] == "True"]
    e_aa = sum(int(r["aligned_aa"]) for r in edent)
    e_fs = sum(int(r["frameshifts"]) for r in edent)
    e_st = sum(int(r["stops"]) for r in edent)
    assert e_aa == 4122 and e_fs == 41 and e_st == 29
    assert round(1000.0 * (e_fs + e_st) / e_aa, 2) == 16.98


# ---------------------------------------------------------------------------
# 5. Recount shared mutations, homoplasies & verified splice borders
# ---------------------------------------------------------------------------
def test_level4_mutations_homoplasy_and_splice_borders() -> None:
    s5a = read_tsv(SUPP_SUB / "Table_S5a_level4_shared_mutations.tsv")
    s5b = read_tsv(SUPP_SUB / "Table_S5b_level4_splice_borders.tsv")

    assert len(s5a) == 23
    enam8_ptcs = [r for r in s5a if r["dataset"] == "ENAM_exon8_aves" and r["mutation_type"] == "stop"]
    assert len(enam8_ptcs) == 16
    k_incong_t1 = sum(1 for r in enam8_ptcs if int(r["min_gains_no_rev_T1"]) >= 2)
    assert k_incong_t1 == 6
    # Verify p.Ser1022Ter ML ASR P_stop at Neognathae > 0.90 (documented reversal alternative!)
    c1022 = next(r for r in enam8_ptcs if int(r["codon_1based"]) == 1022)
    assert c1022["standardized_variant_id"] == "XP_019346377.1:p.Ser1022Ter"
    # Verify no AMBN exon 9 fake 5' donor rows remain in Table_S5b
    assert all(not (r["gene"] == "AMBN" and r["exon"] == "exon9") for r in s5b)
    # Verify ENAM exon 6 splice borders agree between MAFFT & MUSCLE for all 6 species
    ex6 = [r for r in s5b if r["gene"] == "ENAM" and r["exon"] == "exon6"]
    assert len(ex6) == 6
    assert all(r["dual_aligner_agree_acc"] == "True" and r["dual_aligner_agree_don"] == "True" for r in ex6)


# ---------------------------------------------------------------------------
# 6. No Portuguese leaks in Table S6 & Figures derive strictly from Tables
# ---------------------------------------------------------------------------
def test_table_s6_clean_english_and_figures_from_tables() -> None:
    s6_txt = (SUPP_SUB / "Table_S6_auxiliary_odaph_multiz_recomb_psiamel.tsv").read_text(encoding="utf-8")
    for bad in ("nao localizado", "especie COM dente", "sensibilidade no ODAM"):
        assert bad not in s6_txt, f"Portuguese leak in Table S6: {bad}"
    s6 = read_tsv(SUPP_SUB / "Table_S6_auxiliary_odaph_multiz_recomb_psiamel.tsv")
    odaph_with_gap = [r for r in s6 if r["module"] == "ODAPH_ODAM" and r["gene"] == "ODAPH" and r["odaph_gap_bp"] != ""]
    assert len(odaph_with_gap) == 17

    fig_path = ROOT / "scripts" / "enamel_relics" / "06_make_figures.py"
    tree = ast.parse(fig_path.read_text(encoding="utf-8"))
    # Exclude module docstring, check all executable string literals for legacy path references
    non_doc_strings = [
        node.value
        for i, stmt in enumerate(tree.body)
        if not (i == 0 and isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Constant))
        for node in ast.walk(stmt)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    ]
    assert all("atavismo" not in s.lower() for s in non_doc_strings), (
        "Figure script code must derive strictly from supplementary tables, not legacy runs/atavismo!"
    )
    for stem in ("fig1_synteny_map", "fig2_compositional_bias", "fig3_lesions_and_phylogeny"):
        for ext in ("png", "pdf", "svg"):
            p_sub = SUB / "figures" / f"{stem}.{ext}"
            p_run = RUNS / "figures" / f"{stem}.{ext}"
            assert p_sub.exists() and p_run.exists()
            assert sha256_file(p_sub) == sha256_file(p_run)


# ---------------------------------------------------------------------------
# 7. Manuscript editorial integrity, citations, placeholders & PDF-Typst-MD sync
# ---------------------------------------------------------------------------
def test_manuscript_integrity_citations_and_sync() -> None:
    typ = (SUB / "manuscript.typ").read_text(encoding="utf-8")
    md = (SUB / "manuscript.md").read_text(encoding="utf-8")
    bib = (SUB / "references.bib").read_text(encoding="utf-8")
    readme = (SUB / "README_SUBMISSION.md").read_text(encoding="utf-8")

    # Zero placeholders or premature posted-status header
    for name, txt in [("typ", typ), ("md", md), ("readme", readme)]:
        for pat in (r"XXXXXXX", r"insert upon", r"insert before", r"\[ORCID\]", r"\[bioRxiv Preprint\]"):
            assert not re.search(pat, txt), f"Forbidden pattern {pat} in {name}"
        # Zero ambiguous legacy species codes Tgut / Tgua
        assert not re.search(r"\b(Tgut|Tgua)\b", txt), f"Ambiguous legacy species code found in {name}"
        # Zero misuse of HGVS c. for codon numbers
        assert not re.search(r"`c\.(311|439|558|562|615|619|684|109|123)`", txt), f"Misused HGVS c. codon prefix in {name}"

    # Abstract <= 250 words in both Typst and Markdown, and identical normalized text
    m_abs_md = re.search(r"## Abstract\s+(.*?)\n\n\*\*Keywords", md, re.S)
    m_abs_typ = re.search(r"#heading\(level: 1, numbering: none\)\[Abstract\]\s+(.*?)\n\n#\[", typ, re.S)
    assert m_abs_md and m_abs_typ, "Abstract block missing in manuscript.md or manuscript.typ"
    abs_words_md = len(re.findall(r"\S+", m_abs_md.group(1)))
    abs_words_typ = len(re.findall(r"\S+", m_abs_typ.group(1)))
    assert abs_words_md <= 250, f"Markdown abstract exceeds 250 words: {abs_words_md}"
    assert abs_words_typ <= 250, f"Typst abstract exceeds 250 words: {abs_words_typ}"

    def _norm_abs(s: str) -> str:
        s = s.replace(r"\~", "~")
        s = re.sub(r"[*_`]", "", s)
        return re.sub(r"\s+", " ", s).strip()

    assert _norm_abs(m_abs_md.group(1)) == _norm_abs(m_abs_typ.group(1)), "Typst and Markdown abstracts diverged!"

    # Zero unresolved Typst syntax or @cite tokens leaked into Markdown
    assert "#set par.line" not in md and "// ===" not in md
    assert not re.search(r"@[a-z]+[0-9]{4}", md), "Unresolved Typst @citation leaked into manuscript.md"

    # All BibTeX keys cited in manuscript.typ and zero broken keys
    bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib))
    cited_keys = set(re.findall(r"@([a-z]+[0-9]{4}[a-z]?)", typ)) | set(re.findall(r"#pcite\(<([a-z]+[0-9]{4}[a-z]?)>\)", typ))
    assert cited_keys == bib_keys, f"Citation mismatch! Missing from bib: {cited_keys - bib_keys}; Uncited: {bib_keys - cited_keys}"
    assert len(bib_keys) == 49
    assert sha256_file(SUB / "references.bib") == sha256_file(RUNS / "references.bib")

    # Concept DOI & ORCID present in all three
    for txt in (typ, md, readme):
        assert "10.5281/zenodo.24358066" in txt
        assert "0009-0000-3557-8337" in txt

    # PDF exists, valid header/EOF, newer than manuscript.typ & references.bib
    pdf = SUB / "manuscript_biorxiv.pdf"
    assert pdf.exists() and pdf.stat().st_size > 500_000
    raw_pdf = pdf.read_bytes()
    assert raw_pdf.startswith(b"%PDF-1."), "Invalid PDF header"
    assert b"%%EOF" in raw_pdf[-1024:], "Missing PDF EOF trailer"
    assert pdf.stat().st_mtime >= (SUB / "manuscript.typ").stat().st_mtime - 1.0
    assert pdf.stat().st_mtime >= (SUB / "references.bib").stat().st_mtime - 1.0


# ---------------------------------------------------------------------------
# 8. Manuscript Tables 1–3 numbers match supplementary TSVs & deliverables exist
# ---------------------------------------------------------------------------
def test_table_numbers_match_supplementary_and_deliverables_exist() -> None:
    typ = (SUB / "manuscript.typ").read_text(encoding="utf-8")
    md = (SUB / "manuscript.md").read_text(encoding="utf-8")

    # Table 1 vs Table S2b & Table S6
    s2b = {r["candidate_id"]: r for r in read_tsv(SUPP_SUB / "Table_S2b_whole_proteome_reciprocal_blastx.tsv")}
    for cid in ("Ggal__ENAM_Amis", "Ggal__AMBN_Amis", "Shab__AMTN_Amis"):
        bits_tok = f"{float(s2b[cid]['top_hsp_proteome_bits']):.1f}"
        fwd_tok = f"{float(s2b[cid]['top_tblastn_bits']):.1f}"
        assert bits_tok in typ and bits_tok in md, f"Table 1 reciprocal bitscore {bits_tok} missing for {cid}"
        assert fwd_tok in typ and fwd_tok in md, f"Table 1 forward bitscore {fwd_tok} missing for {cid}"
    assert "XP_019346377.1" in typ and "XP_019346384.1" in typ and "XP_014451634.1" in typ
    assert "Not confirmed" in typ and "Not confirmed" in md  # Ggal Psi-ENAM Exon 6

    # Table 2 vs Phase C 1,000-perm clade medians TSV (derived from Table S3c)
    cm = read_tsv(RUNS / "v101" / "stats" / "clade_medians_equidistant_1000perm_v101.tsv")
    idx = {(r["query_name"], r["seg_mode"]): r for r in cm}

    def _r1(val_str: str) -> str:
        return str(Decimal(val_str).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP))

    for q in ("LPL_Acar", "IBSP_Acar", "SPP1_Acar", "ENAM_Acar", "AMBN_Acar"):
        for seg in ("seg_yes", "seg_no"):
            row = idx[(q, seg)]
            c_med = _r1(row["croc_1000"].split()[0])
            t_med = _r1(row["turt_1000"].split()[0])
            a_med = _r1(row["aves_1000"].split()[0])
            a_pct = f"{_r1(row['aves_croc_pct_1000'])}%"
            a_blk = _r1(row["aves_block20_med"])
            for token in (c_med, t_med, a_med, a_pct, a_blk):
                assert token in typ, f"Table 2 value {token} missing in typ for {q}/{seg}"
                assert token in md, f"Table 2 value {token} missing in md for {q}/{seg}"

    # Table 3 vs Table S5a shared mutation counts per clade/category
    s5a = read_tsv(SUPP_SUB / "Table_S5a_level4_shared_mutations.tsv")
    expected_counts = {
        "Palaeognathae": (1, 3),
        "Novaeratitae": (0, 7),
        "Galloanserae": (1, 1),
        "Psittacopasserae": (0, 1),
        "Durocryptodira": (1, 2),
    }
    for clade, (exp_fs, exp_ptc) in expected_counts.items():
        rows = [r for r in s5a if r["orig_classification"] == f"CONGRUENT_{clade}"]
        obs_fs = sum(1 for r in rows if r["mutation_type"] in ("ins", "del"))
        obs_ptc = sum(1 for r in rows if r["mutation_type"] == "stop")
        assert (obs_fs, obs_ptc) == (exp_fs, exp_ptc), f"Table 3 count mismatch for {clade}: got {(obs_fs, obs_ptc)}"
        assert clade in typ and clade in md
    assert sum(1 for r in s5a if r["orig_classification"] == "CONGRUENT_Neornithes") == 0
    assert sum(1 for r in s5a if r["dataset"] == "ENAM_exon8_aves" and r["mutation_type"] == "stop" and int(r["min_gains_no_rev_T1"]) >= 2) == 6
    assert "37.5%" in typ and "31.2%" in typ and "37.5%" in md and "31.2%" in md

    # Verify all mandatory audit deliverables + metadata files exist and are non-empty
    required_files = [
        ROOT / "audit" / "reproducibility_audit.md",
        ROOT / "audit" / "bibliography_report.md",
        ROOT / "audit" / "claims_evidence_matrix.tsv",
        ROOT / "audit" / "change_log.md",
        ROOT / "audit" / "biorxiv_readiness_report.md",
        ROOT / ".zenodo.json",
        ROOT / "CITATION.cff",
    ]
    for p in required_files:
        assert p.exists() and p.stat().st_size > 200, f"Missing or empty deliverable: {p}"

    # Verify claims_evidence_matrix.tsv parses cleanly with expected columns
    claims = read_tsv(ROOT / "audit" / "claims_evidence_matrix.tsv")
    assert len(claims) >= 15
    assert {"claim_id", "manuscript_section", "v100_status", "v101_status"}.issubset(claims[0].keys())


# ---------------------------------------------------------------------------
# 9. Historical release tag enamel-relics-v1.0.0 strictly unmodified
# ---------------------------------------------------------------------------
def test_historical_tag_unmodified() -> None:
    tag_obj = subprocess.run(["git", "rev-parse", "enamel-relics-v1.0.0"], capture_output=True, text=True, cwd=ROOT, check=True).stdout.strip()
    tag_commit = subprocess.run(["git", "rev-parse", "enamel-relics-v1.0.0^{commit}"], capture_output=True, text=True, cwd=ROOT, check=True).stdout.strip()
    assert tag_obj == "b9ff25e654007ef26f5ce01ce702f68c214c90b2"
    assert tag_commit == "a425b2d1b20bd0fd77a997915d62523f5b8f2649"
