# Final Pre-Submission Review Report (`final_review_report.md`)

**Manuscript Title:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles*  
**Repository:** `https://github.com/patoangolano/bio-stack`  
**Audit Branch:** `audit/biorxiv-pre-submission-final` (branched from `enamel-relics-v1.0.2` / `59c1673`)  
**Historical Releases Preserved Untouched:**
- `enamel-relics-v1.0.0`: `b9ff25e654007ef26f5ce01ce702f68c214c90b2` $\to$ `a425b2d1b20bd0fd77a997915d62523f5b8f2649` ([`10.5281/zenodo.24358067`](https://doi.org/10.5281/zenodo.24358067))
- `enamel-relics-v1.0.1`: `13b34a6f77a7b7fdd56d95ee1eac4549d6be9370` ([`10.5281/zenodo.24360150`](https://doi.org/10.5281/zenodo.24360150))
- `enamel-relics-v1.0.2`: `db711560a41948ef5978ec146e1b1bbf7c39bd04` ([`10.5281/zenodo.24368656`](https://doi.org/10.5281/zenodo.24368656))
- **Zenodo Concept DOI:** [`10.5281/zenodo.24358066`](https://doi.org/10.5281/zenodo.24358066)

---

## 1. Executive Summary of Pre-Submission Corrections (`Phases A–H`)

| Phase | Objective | Classification | Detailed Deliverable |
|---|---|---|---|
| **Phase A** | Bound crown-Neornithes synapomorphy claims by empirical exon coverage (`0/878` 7-bird coverage; `34.9%` 6 non-passerine; `82.8%` $\ge 1$ per clade) and quantify ancestral-stop survival/detection power (`16.3%` uncontradicted cross-clade detection; `40.8%` secondary sense-reversal corruption) | `CORRIGIDO` | [`neornithes_detection_limits.md`](neornithes_detection_limits.md) |
| **Phase B** | Reframe coverage-masked K80 simulations (`Table_S8`) from causal proof to model- and mask-conditional statistical compatibility, documenting partial conditioning on the empirical coverage mask and `T2`/`T3` upper-tail behavior ($p = 0.072–0.244$) | `CORRIGIDO` | [`k80_interpretation_review.md`](k80_interpretation_review.md) |
| **Phase C** | Classify marginal empirical $p = 0.024$ (`Strhab` `ENAM_Acar`) across the 504-cell grid as unadjusted exploratory without multiplicity control; reinforce single-locus site-linkage caveat on Clopper–Pearson 95% CIs; separate external phylogenomic support for `T1` from within-locus RELL tests | `CORRIGIDO` | [`statistical_interpretation_review.md`](statistical_interpretation_review.md) |
| **Phase D** | Calibrate overly conclusive verbs (`"prove"`, `"ruling out"`, `"unequivocally"`, `"independently favors"`), delimit the single globular *LPL* control, and preserve graded orthology tiers (`Table_S2b`), Gnomon native-CDS categories (`Table_S7`), and variant nomenclature (`Table_S5a`) | `CORRIGIDO` / `PASSOU` | [`editorial_changes.md`](editorial_changes.md) |
| **Phase E** | Finalize 21-page PDF (`manuscript_biorxiv.pdf`, lines `1–679`), 248-word Abstract ($\le 250$), Figure 3B Crown Neornithes badge (`"none detected (0/878 7-bird cov.)"`), Table 3, and verbatim `manuscript.typ`/`manuscript.md` parity | `PASSOU` | [`editorial_changes.md`](editorial_changes.md) |
| **Phase F** | Document express openRxiv rules vs. prudential precautions, provide Option A and Option B headers, field-by-field bioRxiv portal guide, screener cover note, and 7-point human review checklist | `PENDENTE` (Human Decision) | [`affiliation_and_submission.md`](affiliation_and_submission.md) |
| **Phase G** | Preserve historical tags (`v1.0.0`, `v1.0.1`, `v1.0.2`) and Zenodo deposits untouched; compute diffs and SHA-256 checksums for candidate `v1.0.3` without premature release publication | `PASSOU` | [`final_change_log.md`](final_change_log.md) |
| **Phase H** | Expand and run automated `pytest` acceptance suite (`13/13` tests passed) and `ruff` lint checks (`0` errors) | `PASSOU` | [`test_results.md`](test_results.md) |

---

## 2. Complete Index of `audit/final/` Deliverables

1. [`final_review_report.md`](final_review_report.md) — Executive synthesis of the final pre-submission review.
2. [`neornithes_detection_limits.md`](neornithes_detection_limits.md) — Phase A audit of cross-clade coverage and ancestral-stop detection power.
3. [`k80_interpretation_review.md`](k80_interpretation_review.md) — Phase B epistemological and methodological review of K80 simulations.
4. [`statistical_interpretation_review.md`](statistical_interpretation_review.md) — Phase C review of exploratory *p*-values, site linkage, and topology tests.
5. [`editorial_changes.md`](editorial_changes.md) — Phases D & E log of calibrated language and 21-page PDF verification.
6. [`affiliation_and_submission.md`](affiliation_and_submission.md) — Phase F affiliation options, bioRxiv form guide, and human review checklist.
7. [`claims_evidence_matrix.tsv`](claims_evidence_matrix.tsv) — Tab-separated claims and requirements matrix (`R01`–`R14`).
8. [`test_results.md`](test_results.md) — Phase H automated `pytest` (`13/13`) and `ruff` execution results.
9. [`final_change_log.md`](final_change_log.md) — Phase G diff and SHA-256 checksum manifest relative to `enamel-relics-v1.0.2`.
10. [`submission_readiness.md`](submission_readiness.md) — Three-dimensional readiness decision (documentary, scientific, administrative).
