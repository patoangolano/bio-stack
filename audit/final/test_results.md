# Phase H Audit Report — Final Automated Acceptance Test Results (`test_results.md`)

**Status:** `PASSOU` (`13/13` automated `pytest` test functions passed; `0` `ruff` lint errors)  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final`  
**Test Suite:** [`tests/test_enamel_relics_audit.py`](../../tests/test_enamel_relics_audit.py)

---

## 1. Automated Test Matrix (`13/13` Passed)

| # | Test Function | Verification Scope | Status |
|---|---|---|---|
| 1 | `test_supplementary_hashes_and_mirror` | SHA-256 hashes, byte sizes, strict `LF` line endings, and exact mirroring across all 14 supplementary TSVs (`Table_S0`–`Table_S9`) and `MANIFEST_SHA256.tsv` | `PASSOU` |
| 2 | `test_coordinates_and_indexing_conventions` | 1-based inclusive and 0-based BED coordinates for chicken `GRCg7b` `NRG1–LPL` (`105,217 bp` strict vs. `105,218 bp` endpoint difference), `ΨENAM` Exon 8 (`1,985-bp` union), `ΨENAM` Exon 6 (`149-bp` span / `47-bp` core, negative reciprocal `BLASTx`), `ΨAMBN` Exon 3, and kakapo `ΨAMTN` | `PASSOU` |
| 3 | `test_species_and_sequence_counts` | 19 taxa (`18` sauropsids + *Homo sapiens*), unambiguous 6-letter species codes (`Taegut` vs. `Tingut`), 16 query profiles, `504` permutation cells (`18 × 14 × 2`), and 4 Gnomon `PROTEIN_CODING` models | `PASSOU` |
| 4 | `test_level3_lesion_rates_and_133_lpl_introns` | Level 3 `miniprot` lesion rates (`18/18` *LPL* controls = `0.00/1,000 aa`, `133/133` `GT-AG` introns; `15` toothed enamel alignments = `0.11/1,000 aa`; edentulous enamel pseudogenes = `16.98/1,000 aa`) | `PASSOU` |
| 5 | `test_level4_mutations_homoplasy_and_splice_borders` | 23 dual-aligner consensus shared mutations (`Table_S5a`), 16 avian *ENAM* Exon 8 PTCs (`6/16` incongruent on `T1`), and flanked-MSA splice borders (`Table_S5b`) | `PASSOU` |
| 6 | `test_table_s6_clean_english_and_figures_from_tables` | Clean English in `Table_S6`, strict derivation of `06_make_figures.py` from `supplementary/Table_S*.tsv`, and SHA-256 parity of all 9 figure files (`.png`, `.pdf`, `.svg`) | `PASSOU` |
| 7 | `test_manuscript_integrity_citations_and_sync` | Zero placeholders, zero ambiguous legacy codes, zero misused HGVS `c.` prefixes, Abstract word count (`248 <= 250` words) and verbatim Typst–Markdown synchronization, `49/49` BibTeX citation parity, Concept DOI & ORCID presence, and compiled PDF integrity | `PASSOU` |
| 8 | `test_table_numbers_match_supplementary_and_deliverables_exist` | Exact numerical correspondence between Manuscript Tables 1–3 and `Table_S2b`, `Table_S3c`, and `Table_S5a` | `PASSOU` |
| 9 | `test_historical_tags_unmodified` | Immutability of all three historical release tags: `enamel-relics-v1.0.0` (`b9ff25e6...`), `enamel-relics-v1.0.1` (`13b34a6f...`), and `enamel-relics-v1.0.2` (`db711560...`) | `PASSOU` |
| 10 | `test_v102_codon_to_amino_acid_translations` | Exact reference residue and standard genetic-code translation verification across all 23 rows of `Table_S5a` (`XP_019346377.1` and `XP_019346386.1`), including codons `601` (`TTT=Phe`, Neornithes ML `TAT=Tyr → TAA=Ter`) and `1003` (`TGG=Trp`, site-independent ML `TGC=Cys → TGA=Ter`) | `PASSOU` |
| 11 | `test_v102_synapomorphy_tiers_and_topological_sensitivity` | Sister-clade coverage tier consistency (`synapomorphy_evidence_tier`), ultra-sensitive `Anapla` *AMBN* Exon 4 alignment files, formal `IQ-TREE 3` RELL KH/SH/AU p-values, and `T1` (`6/16`, `5/16`), `T2` (`10/16`, `9/16`), `T3` (`11/16`, `10/16`) homoplasy counts (`Table_S9`) | `PASSOU` |
| 12 | `test_v102_neutral_simulation_and_deliverables` | 48-scenario K80 neutral simulation sensitivity grid (`Table_S8`): unmasked (`25.0` median, $p = 0.052$) vs. empirically coverage-masked (`16.0` median, $p = 1.000$; `40.0%` unconditional incongruence; `37.5%` conditional median at $N_{\text{shared}}=16$; `0.0` crown-Neornithes stops) | `PASSOU` |
| 13 | `test_pre_submission_final_neornithes_detectability_and_language` | **Phases A–H Final Pre-Submission Gate:** Verifies `neornithes_detectability_audit.json` (`0/878` 7-bird coverage, `306/878` 6 non-passerine, `727/878` $\ge 1$ per clade, `16.3%` neutral uncontradicted ancestral-stop detection power, `40.8%` secondary sense-reversal corruption), absence of all 9 forbidden overconclusive phrases across `manuscript.typ`, `manuscript.md`, and `06_make_figures.py`, presence of all 6 required methodological caveats (`no crown-Neornithes-wide synapomorphy was detected`, `16.3%`, `40.8%`, `unadjusted exploratory empirical`, `conditional on the K80 model`, `conditional descriptions of finite-site sampling within one locus`), and presence of all 10 deliverables in `audit/final/` | `PASSOU` |

*Note on Scope:* Internal automated consistency tests verify deterministic reproducibility, mathematical and translational invariants, hash integrity, and text–table–figure synchronization, and are distinct from out-of-scope external wet-lab or new genome-assembly validations (which are classified as `NÃO VERIFICADO` in `claims_evidence_matrix.tsv` and declared as Limitations in Section 3.5).

---

## 2. Real Execution Output

```text
$ uv run pytest -v tests/test_enamel_relics_audit.py
============================= test session starts =============================
platform win32 -- Python 3.11.18, pytest-9.1.1, pluggy-1.6.0 -- C:\Users\user\bio-stack\.venv\Scripts\python.exe
cachedir: .pytest_cache
rootdir: C:\Users\user\bio-stack
configfile: pyproject.toml
collecting ... collected 13 items

tests/test_enamel_relics_audit.py::test_supplementary_hashes_and_mirror PASSED [  7%]
tests/test_enamel_relics_audit.py::test_coordinates_and_indexing_conventions PASSED [ 15%]
tests/test_enamel_relics_audit.py::test_species_and_sequence_counts PASSED [ 23%]
tests/test_enamel_relics_audit.py::test_level3_lesion_rates_and_133_lpl_introns PASSED [ 30%]
tests/test_enamel_relics_audit.py::test_level4_mutations_homoplasy_and_splice_borders PASSED [ 38%]
tests/test_enamel_relics_audit.py::test_table_s6_clean_english_and_figures_from_tables PASSED [ 46%]
tests/test_enamel_relics_audit.py::test_manuscript_integrity_citations_and_sync PASSED [ 53%]
tests/test_enamel_relics_audit.py::test_table_numbers_match_supplementary_and_deliverables_exist PASSED [ 61%]
tests/test_enamel_relics_audit.py::test_historical_tags_unmodified PASSED [ 69%]
tests/test_enamel_relics_audit.py::test_v102_codon_to_amino_acid_translations PASSED [ 76%]
tests/test_enamel_relics_audit.py::test_v102_synapomorphy_tiers_and_topological_sensitivity PASSED [ 84%]
tests/test_enamel_relics_audit.py::test_v102_neutral_simulation_and_deliverables PASSED [ 92%]
tests/test_enamel_relics_audit.py::test_pre_submission_final_neornithes_detectability_and_language PASSED [100%]

============================= 13 passed in 0.30s ==============================

$ uv run ruff check scripts/enamel_relics/06_make_figures.py scripts/enamel_relics/11_consolidate_v101_tables.py scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py scripts/enamel_relics/13_audit_final_neornithes_and_stats.py tests/test_enamel_relics_audit.py
All checks passed!
```
