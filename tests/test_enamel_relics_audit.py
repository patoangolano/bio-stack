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
        manifest_path = base / "MANIFEST_SHA256.tsv"
        assert b"\r\n" not in manifest_path.read_bytes(), f"CRLF line endings found in {manifest_path}"
        manifest = read_tsv(manifest_path)
        assert len(manifest) == 14, f"Expected 14 supplementary TSVs (Table_S0..Table_S9) in {base}, got {len(manifest)}"
        for row in manifest:
            p = base / row["file"]
            assert p.exists(), f"Missing {p}"
            assert b"\r\n" not in p.read_bytes(), f"CRLF line endings found in {p}"
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
# 9. Historical release tags enamel-relics-v1.0.0, v1.0.1, and v1.0.2 unmodified
# ---------------------------------------------------------------------------
def test_historical_tags_unmodified() -> None:
    tag_v100_obj = subprocess.run(["git", "rev-parse", "enamel-relics-v1.0.0"], capture_output=True, text=True, cwd=ROOT, check=True).stdout.strip()
    tag_v100_commit = subprocess.run(["git", "rev-parse", "enamel-relics-v1.0.0^{commit}"], capture_output=True, text=True, cwd=ROOT, check=True).stdout.strip()
    assert tag_v100_obj == "b9ff25e654007ef26f5ce01ce702f68c214c90b2"
    assert tag_v100_commit == "a425b2d1b20bd0fd77a997915d62523f5b8f2649"

    tag_v101_commit = subprocess.run(["git", "rev-parse", "enamel-relics-v1.0.1^{commit}"], capture_output=True, text=True, cwd=ROOT, check=True).stdout.strip()
    assert tag_v101_commit == "13b34a6f77a7b7fdd56d95ee1eac4549d6be9370"

    tag_v102_commit = subprocess.run(["git", "rev-parse", "enamel-relics-v1.0.2^{commit}"], capture_output=True, text=True, cwd=ROOT, check=True).stdout.strip()
    assert tag_v102_commit == "db711560a41948ef5978ec146e1b1bbf7c39bd04"


# ---------------------------------------------------------------------------
# 10. v1.0.2 Phase A: Codon-to-amino-acid translation verification (Table S5a)
# ---------------------------------------------------------------------------
CODON_TO_AA3: dict[str, str] = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Ter", "TAG": "Ter",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Ter", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}
AA1_TO_AA3: dict[str, str] = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
}


def test_v102_codon_to_amino_acid_translations() -> None:
    # Load reference protein sequences
    faa_path = RUNS / "queries" / "all_real_queries.faa"
    seqs: dict[str, str] = {}
    cur_id = ""
    for line in faa_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            cur_id = line[1:].split()[0]
            seqs[cur_id] = ""
        elif cur_id:
            seqs[cur_id] += line.strip()

    ref_seqs = {
        "XP_019346377.1": seqs["ENAM_Amis"],
        "XP_019346386.1": seqs["AMBN_Amis"],
    }
    assert len(ref_seqs["XP_019346377.1"]) == 1078
    assert len(ref_seqs["XP_019346386.1"]) == 390

    s5a = read_tsv(SUPP_SUB / "Table_S5a_level4_shared_mutations.tsv")
    for r in s5a:
        prot = r["ref_protein"]
        cod_idx = int(r["codon_1based"])
        expected_ref_aa3 = AA1_TO_AA3[ref_seqs[prot][cod_idx - 1]]
        assert r["ref_aa_3letter"] == expected_ref_aa3, f"Ref AA mismatch at {r['standardized_variant_id']}"

        if r["mutation_type"] == "stop":
            assert r["standardized_variant_id"] == f"{prot}:p.{expected_ref_aa3}{cod_idx}Ter"
            assert CODON_TO_AA3[r["amis_codon"]] == expected_ref_aa3
            for c_stop in r["carrier_codons"].split(","):
                assert CODON_TO_AA3[c_stop] == "Ter", f"Carrier codon {c_stop} is not a stop codon"
            if r["neornithes_ml_codon"] not in ("", "-"):
                assert CODON_TO_AA3[r["neornithes_ml_codon"]] == r["neornithes_ml_aa_3letter"]

    # Specifically verify Codon 601 (Phe -> Tyr -> Ter) and Codon 1003 (Trp / site-independent ML Cys -> Ter)
    c601 = next(r for r in s5a if r["standardized_variant_id"] == "XP_019346377.1:p.Phe601Ter")
    assert c601["amis_codon"] == "TTT" and c601["ref_aa_3letter"] == "Phe"
    assert c601["neornithes_ml_codon"] == "TAT" and c601["neornithes_ml_aa_3letter"] == "Tyr"
    assert c601["carrier_codons"] == "TAA"

    c1003 = next(r for r in s5a if r["standardized_variant_id"] == "XP_019346377.1:p.Trp1003Ter")
    assert c1003["amis_codon"] == "TGG" and c1003["ref_aa_3letter"] == "Trp"
    assert c1003["neornithes_ml_codon"] == "TGC" and c1003["neornithes_ml_aa_3letter"] == "Cys"
    assert c1003["carrier_codons"] == "TGA"


# ---------------------------------------------------------------------------
# 11. v1.0.2 Phase B & C: Synapomorphy coverage tiers & topological sensitivity
# ---------------------------------------------------------------------------
def test_v102_synapomorphy_tiers_and_topological_sensitivity() -> None:
    s5a = read_tsv(SUPP_SUB / "Table_S5a_level4_shared_mutations.tsv")
    s9 = {r["topology_id"]: r for r in read_tsv(SUPP_SUB / "Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv")}

    # Every variant with missing sister-clade coverage must NOT be labeled OUTGROUP_AND_SISTER_VERIFIED_SYNAPOMORPHY
    for r in s5a:
        tier = r["synapomorphy_evidence_tier"]
        assert tier != "", f"Missing synapomorphy_evidence_tier on {r['standardized_variant_id']}"
        if r["covered_without_6letter"] in ("", "-"):
            assert "SISTER_UNRECOVERED" in tier, f"Variant {r['standardized_variant_id']} lacks sister coverage but tier is {tier}"
        if tier.startswith("OUTGROUP_AND_SISTER_VERIFIED_SYNAPOMORPHY"):
            assert r["covered_without_6letter"] not in ("", "-")

    # Verify ultra-sensitive Anapla AMBN Exon 4 alignment audit files exist
    ambn_dir = RUNS / "v102" / "ambn_exon4_coverage_audit"
    assert (ambn_dir / "mafft_apla.fna").exists() and (ambn_dir / "muscle_apla.fna").exists()

    # Verify T1, T2, T3 homoplasy counts in Table S5a match Table S9
    enam8_ptcs = [r for r in s5a if r["dataset"] == "ENAM_exon8_aves" and r["mutation_type"] == "stop"]
    assert sum(1 for r in enam8_ptcs if int(r["min_gains_no_rev_T1"]) >= 2) == 6
    assert sum(1 for r in enam8_ptcs if int(r["min_gains_no_rev_T2"]) >= 2) == 10
    assert sum(1 for r in enam8_ptcs if int(r["min_gains_no_rev_T3"]) >= 2) == 11

    assert int(s9["T1"]["k_incong_irreversible_parsimony"]) == 6 and float(s9["T1"]["pct_incong_irreversible_parsimony"]) == 37.5
    assert int(s9["T1"]["k_incong_reversible_ml_asr"]) == 5 and float(s9["T1"]["pct_incong_reversible_ml_asr"]) == 31.2
    assert int(s9["T2"]["k_incong_irreversible_parsimony"]) == 10 and float(s9["T2"]["pct_incong_irreversible_parsimony"]) == 62.5
    assert int(s9["T2"]["k_incong_reversible_ml_asr"]) == 9 and float(s9["T2"]["pct_incong_reversible_ml_asr"]) == 56.2
    assert int(s9["T3"]["k_incong_irreversible_parsimony"]) == 11 and float(s9["T3"]["pct_incong_irreversible_parsimony"]) == 68.8
    assert int(s9["T3"]["k_incong_reversible_ml_asr"]) == 10 and float(s9["T3"]["pct_incong_reversible_ml_asr"]) == 62.5

    # Verify formal IQ-TREE 3 KH/SH/AU p-values reject T2 and T3 in both aligners
    for topo in ("T2", "T3"):
        assert float(s9[topo]["mafft_p_SH"]) <= 0.0040
        assert float(s9[topo]["mafft_p_AU"]) <= 0.0018
        assert float(s9[topo]["muscle_p_SH"]) <= 0.0040
        assert float(s9[topo]["muscle_p_AU"]) <= 0.0018


# ---------------------------------------------------------------------------
# 12. v1.0.2 Phase D & Deliverables: Neutral K80 coverage-masked simulation
# ---------------------------------------------------------------------------
def test_v102_neutral_simulation_and_deliverables() -> None:
    s8 = read_tsv(SUPP_SUB / "Table_S8_k80_simulation_sensitivity_v102.tsv")
    assert len(s8) == 48

    # Unmasked full 878-codon baseline overpredicts shared PTCs (25-28 median, p <= 0.052)
    unmasked_t0 = next(
        r for r in s8
        if r["topology"] == "T1" and r["anc_seq_source"] == "Neornithes_ML"
        and r["mask_mode"] == "unmasked_full_878" and float(r["kappa"]) == 4.0
        and float(r["branch_scale"]) == 1.0 and float(r["t_stem"]) == 0.0
    )
    assert float(unmasked_t0["shared_ptc_median"]) == 25.0
    assert float(unmasked_t0["obs16_p_two_sided"]) == 0.052

    # Empirically masked dual-aligner baseline matches observed 16 shared PTCs (16.0 median, p = 1.000) and 37.5% incongruence
    masked_t0 = next(
        r for r in s8
        if r["topology"] == "T1" and r["anc_seq_source"] == "Neornithes_ML"
        and r["mask_mode"] == "empirical_dual_aligner_mask" and float(r["kappa"]) == 4.0
        and float(r["branch_scale"]) == 1.0 and float(r["t_stem"]) == 0.0
    )
    assert float(masked_t0["shared_ptc_median"]) == 16.0
    assert float(masked_t0["obs16_p_two_sided"]) == 1.0
    assert float(masked_t0["incong_rate_median"]) == 0.4
    assert float(masked_t0["obs_incong_p_two_sided"]) == 0.848
    assert float(masked_t0["cond_N16_incong_median"]) == 0.375
    assert float(masked_t0["crown_neornithes_syn_median"]) == 0.0

    # Verify v1.0.2 audit deliverables exist
    v102_dir = ROOT / "audit" / "v102"
    for fname in (
        "final_scientific_audit.md",
        "variant_nomenclature.md",
        "phylogenetic_robustness.md",
        "neutral_simulation_report.md",
        "claims_evidence_matrix.tsv",
        "test_results.md",
        "change_log.md",
        "biorxiv_readiness_report.md",
        "human_approval_checklist.md",
    ):
        p = v102_dir / fname
        assert p.exists() and p.stat().st_size > 200, f"Missing or empty v1.0.2 deliverable: {p}"

    v102_claims = read_tsv(v102_dir / "claims_evidence_matrix.tsv")
    assert len(v102_claims) >= 16
    assert {"claim_id", "finding_classification", "readiness_status"}.issubset(v102_claims[0].keys())


# ---------------------------------------------------------------------------
# 13. Pre-submission final Phase A-F: Neornithes detectability, K80 conditioning,
#     exploratory p-values, calibrated language, and audit/final/ deliverables
# ---------------------------------------------------------------------------
def test_pre_submission_final_neornithes_detectability_and_language() -> None:
    typ = (SUB / "manuscript.typ").read_text(encoding="utf-8")
    md = (SUB / "manuscript.md").read_text(encoding="utf-8")
    fig_code = (ROOT / "scripts" / "enamel_relics" / "06_make_figures.py").read_text(encoding="utf-8")

    # Verify neornithes_detectability_audit.json metrics
    det_path = RUNS / "v102" / "neornithes_detectability_audit.json"
    assert det_path.exists()
    det = json.loads(det_path.read_text(encoding="utf-8"))
    enam8_dual = det["exon_coverage_breakdown"]["ENAM_exon8_aves"]["dual_aligner"]
    assert enam8_dual["all_7_birds_covered"] == 0
    assert enam8_dual["all_6_non_taegut_covered"] == 306
    assert enam8_dual["cross_clade_ge1_pal_ge1_neog"] == 727
    assert enam8_dual["cross_clade_ge2_pal_ge2_neog"] == 675

    power_dual_1 = next(
        p for p in det["ancestral_neornithes_stop_detection_power"]
        if p["mask_mode"] == "empirical_dual_aligner_mask" and float(p["branch_scale"]) == 1.0
    )
    assert float(power_dual_1["p_strict_7of7_universal"]) == 0.0
    assert float(power_dual_1["p_uncontradicted_ge1_pal_ge1_neog_TableS8_def"]) == 0.1625
    assert float(power_dual_1["p_reversal_corrupted_cross_clade_appears_incongruent"]) == 0.4082

    # Forbidden overconclusive phrases must NOT appear in manuscript.typ, manuscript.md, or 06_make_figures.py
    forbidden_phrases = (
        "zero crown Neornithes synapomorphies",
        "0 congruent lesions",
        "Zero shared lesions span Palaeognathae + Neognathae",
        "ruling out generalized avian substitution-rate acceleration",
        "ruling out uniform avian rate acceleration",
        "confirming `Strhab` `ENAM_Acar`",
        "ENAM_ Exon 8 independently favors",
        "ENAM* Exon 8 independently favors",
        "unequivocally confirms",
    )
    for name, txt in (("manuscript.typ", typ), ("manuscript.md", md), ("06_make_figures.py", fig_code)):
        for phrase in forbidden_phrases:
            assert phrase not in txt, f"Overconclusive phrase '{phrase}' still present in {name}"

    # Required calibrated caveats must appear in both manuscript.typ and manuscript.md
    required_caveats = (
        "no crown-Neornithes-wide synapomorphy was detected",
        "16.3%",
        "40.8%",
        "unadjusted exploratory empirical",
        "conditional on the K80 model",
        "conditional descriptions of finite-site sampling within one locus",
    )
    for name, txt in (("manuscript.typ", typ), ("manuscript.md", md)):
        for req in required_caveats:
            assert req in txt, f"Required caveat '{req}' missing from {name}"

    assert "none detected (0/878" in fig_code

    # Verify all 10 mandatory deliverables in audit/final/ exist and are non-empty
    final_dir = ROOT / "audit" / "final"
    for fname in (
        "final_review_report.md",
        "neornithes_detection_limits.md",
        "k80_interpretation_review.md",
        "statistical_interpretation_review.md",
        "editorial_changes.md",
        "affiliation_and_submission.md",
        "claims_evidence_matrix.tsv",
        "test_results.md",
        "final_change_log.md",
        "submission_readiness.md",
    ):
        p = final_dir / fname
        assert p.exists() and p.stat().st_size > 200, f"Missing or empty final deliverable: {p}"

    final_claims = read_tsv(final_dir / "claims_evidence_matrix.tsv")
    assert len(final_claims) >= 14
    assert {"claim_id", "status", "verification_category"}.issubset(final_claims[0].keys())
