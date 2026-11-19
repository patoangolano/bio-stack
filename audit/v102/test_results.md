# Phase G Audit Report — Automated Acceptance Test Suite & Lint Verification (`v1.0.2`)

**Status:** `PASSOU` (`12/12` automated tests passed; `0` ruff lint errors)  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-final-v1.0.2`  
**Test Suite:** [`tests/test_enamel_relics_audit.py`](../../tests/test_enamel_relics_audit.py)

---

## 1. Automated Test Coverage (`tests/test_enamel_relics_audit.py`)

| # | Test Function | Phase & Verification Scope | Result |
|---|---|---|---|
| 1 | `test_supplementary_hashes_and_mirror` | Verifies all 14 consolidated supplementary TSVs (`Table_S0` through `Table_S9`) against `MANIFEST_SHA256.tsv` (exact byte sizes, SHA-256 hashes, strict `LF` line endings, and exact mirroring between `runs/enamel_relics/supplementary/` and `submission_biorxiv/supplementary/`) | `PASSOU` |
| 2 | `test_coordinates_and_indexing_conventions` | Verifies 1-based inclusive and 0-based BED coordinates for *Gallus gallus* `GRCg7b` `NRG1–LPL` intergenic interval (`105,217 bp` strict vs. `105,218 bp` endpoint difference), `ΨENAM` Exon 8 (`1,985-bp` union vs. `1,975-bp` `BLASTn` span), `ΨENAM` Exon 6 (`149-bp` flanked span vs. `47-bp` core and negative whole-proteome reciprocal `BLASTx`), `ΨAMBN` Exon 3, and kakapo `ΨAMTN` | `PASSOU` |
| 3 | `test_species_and_sequence_counts` | Verifies 19 taxa in `Table_S0`/`Table_S1` (18 sauropsids + *Homo sapiens*), unambiguous 6-letter species codes (`Taegut` vs. `Tingut`), 16 query profiles in `Table_S3a`, `504` permutation cells (`18 × 14 × 2`) in `Table_S3c`, and 4 Gnomon `PROTEIN_CODING` models in `Table_S7` | `PASSOU` |
| 4 | `test_level3_lesion_rates_and_133_lpl_introns` | Recomputes Level 3 `miniprot` lesion rates from `Table_S4`: `18/18` *LPL* positive controls (`8,246 aa`, `0` frameshifts, `0` stops, `133/133` canonical `GT-AG` introns), `15` toothed enamel alignments (`9,394 aa`, `0.11 lesions/1,000 aa`), and edentulous enamel pseudogenes (`4,122 aa`, `41` frameshifts + `29` stops = `16.98 lesions/1,000 aa`) | `PASSOU` |
| 5 | `test_level4_mutations_homoplasy_and_splice_borders` | Verifies 23 dual-aligner consensus shared mutations in `Table_S5a`, 16 avian *ENAM* Exon 8 PTCs (`6/16` incongruent on `T1`), and flanked-MSA splice borders in `Table_S5b` | `PASSOU` |
| 6 | `test_table_s6_clean_english_and_figures_from_tables` | Confirms clean English in `Table_S6`, strict derivation of `06_make_figures.py` from `supplementary/Table_S*.tsv`, and SHA-256 identity of all 9 figure files (`.png`, `.pdf`, `.svg`) across `runs/` and `submission_biorxiv/` | `PASSOU` |
| 7 | `test_manuscript_integrity_citations_and_sync` | Verifies zero placeholders, zero ambiguous legacy species codes, zero misused HGVS `c.` codon prefixes, Abstract word count (`247 <= 250` words) and verbatim synchronization between `manuscript.typ` and `manuscript.md`, `49/49` BibTeX citation consistency, Concept DOI `10.5281/zenodo.24358066` & ORCID `0009-0000-3557-8337`, and compiled `manuscript_biorxiv.pdf` validity | `PASSOU` |
| 8 | `test_table_numbers_match_supplementary_and_deliverables_exist` | Verifies every number in Manuscript Tables 1, 2, and 3 against `Table_S2b`, `Table_S3c`, and `Table_S5a`, plus existence of all baseline `v1.0.1` audit deliverables | `PASSOU` |
| 9 | `test_historical_tags_unmodified` | Verifies immutability of both historical Git tags: `enamel-relics-v1.0.0` (`b9ff25e654007ef26f5ce01ce702f68c214c90b2` → `a425b2d1b20bd0fd77a997915d62523f5b8f2649`) and `enamel-relics-v1.0.1` (`13b34a6f77a7b7fdd56d95ee1eac4549d6be9370`) | `PASSOU` |
| 10 | `test_v102_codon_to_amino_acid_translations` | **Phase A (`v1.0.2`):** Programmatically verifies every variant in `Table_S5a` against full-length reference proteins `XP_019346377.1` (*ENAM*, 1,078 aa) and `XP_019346386.1` (*AMBN*, 390 aa) and standard genetic-code translation of `amis_codon`, `neornithes_ml_codon`, and `carrier_codons` (explicitly asserting codon 601 `TTT=Phe`, Neornithes ML `TAT=Tyr → TAA=Ter` and codon 1003 `TGG=Trp`, site-independent ML `TGC=Cys → TGA=Ter`) | `PASSOU` |
| 11 | `test_v102_synapomorphy_tiers_and_topological_sensitivity` | **Phases B & C (`v1.0.2`):** Verifies `synapomorphy_evidence_tier` consistency with sister-clade coverage (`covered_without_6letter`), existence of ultra-sensitive `Anapla` *AMBN* Exon 4 alignment files, formal `IQ-TREE 3.1.4` KH/SH/AU topology test rejection of `T2`/`T3` ($p_{\text{SH}} \le 0.0040, p_{\text{AU}} \le 0.0018$), and exact parsimony/ML homoplasy counts across `T1` (`6/16`, `5/16`), `T2` (`10/16`, `9/16`), and `T3` (`11/16`, `10/16`) in `Table_S9` | `PASSOU` |
| 12 | `test_v102_neutral_simulation_and_deliverables` | **Phase D & Deliverables (`v1.0.2`):** Verifies 48-scenario K80 neutral simulation sensitivity grid (`Table_S8`): confirms unmasked 878-codon overprediction (`25.0` median, $p = 0.052$) vs. empirically coverage-masked match (`16.0` median, $p = 1.000$; `40.0%` unconditional incongruence, $p = 0.848$; `37.5%` conditional median at $N_{\text{shared}}=16$; `0.0` crown-Neornithes synapomorphies), and verifies all 9 required deliverables in `audit/v102/` | `PASSOU` |

---

## 2. Exact Execution Logs

### 2.1. Pytest Execution (`uv run pytest -v tests/test_enamel_relics_audit.py`)

```text
============================= test session starts =============================
platform win32 -- Python 3.11.18, pytest-9.1.1, pluggy-1.6.0 -- C:\Users\user\bio-stack\.venv\Scripts\python.exe
cachedir: .pytest_cache
rootdir: C:\Users\user\bio-stack
configfile: pyproject.toml
collecting ... collected 12 items

tests/test_enamel_relics_audit.py::test_v102_neutral_simulation_and_deliverables PASSED [  8%]
tests/test_enamel_relics_audit.py::test_supplementary_hashes_and_mirror PASSED [ 16%]
tests/test_enamel_relics_audit.py::test_species_and_sequence_counts PASSED [ 25%]
tests/test_enamel_relics_audit.py::test_level3_lesion_rates_and_133_lpl_introns PASSED [ 33%]
tests/test_enamel_relics_audit.py::test_level4_mutations_homoplasy_and_splice_borders PASSED [ 41%]
tests/test_enamel_relics_audit.py::test_historical_tags_unmodified PASSED [ 50%]
tests/test_enamel_relics_audit.py::test_table_numbers_match_supplementary_and_deliverables_exist PASSED [ 58%]
tests/test_enamel_relics_audit.py::test_manuscript_integrity_citations_and_sync PASSED [ 66%]
tests/test_enamel_relics_audit.py::test_table_s6_clean_english_and_figures_from_tables PASSED [ 75%]
tests/test_enamel_relics_audit.py::test_v102_synapomorphy_tiers_and_topological_sensitivity PASSED [ 83%]
tests/test_enamel_relics_audit.py::test_v102_codon_to_amino_acid_translations PASSED [ 91%]
tests/test_enamel_relics_audit.py::test_coordinates_and_indexing_conventions PASSED [100%]

============================= 12 passed in 0.29s ==============================
```

### 2.2. Ruff Lint Verification (`uv run ruff check`)

```text
$ uv run ruff check scripts/enamel_relics/06_make_figures.py scripts/enamel_relics/11_consolidate_v101_tables.py scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py tests/test_enamel_relics_audit.py
All checks passed!
```
