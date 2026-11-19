# Phase G Audit Report — Final Change Log, SHA-256 Checksums, and Version Control (`v1.0.2` → Candidate `v1.0.3`)

**Status:** `PASSOU`  
**Audit Date:** 2026-11-19  
**Baseline Release:** `enamel-relics-v1.0.2` (tag `db711560a41948ef5978ec146e1b1bbf7c39bd04` / commit `59c1673`, Zenodo Version DOI `10.5281/zenodo.24368656`, Concept DOI `10.5281/zenodo.24358066`)  
**Target Branch:** `audit/biorxiv-pre-submission-final`

---

## 1. Historical Releases & Zenodo Deposits Preserved Untouched

| Release Tag | Git Commit / Tag SHA | Zenodo Version DOI | Status |
|---|---|---|---|
| `enamel-relics-v1.0.0` | `b9ff25e654007ef26f5ce01ce702f68c214c90b2` $\to$ `a425b2d1b20bd0fd77a997915d62523f5b8f2649` | [`10.5281/zenodo.24358067`](https://doi.org/10.5281/zenodo.24358067) | **Intact (`PASSOU`)** |
| `enamel-relics-v1.0.1` | `13b34a6f77a7b7fdd56d95ee1eac4549d6be9370` | [`10.5281/zenodo.24360150`](https://doi.org/10.5281/zenodo.24360150) | **Intact (`PASSOU`)** |
| `enamel-relics-v1.0.2` | `db711560a41948ef5978ec146e1b1bbf7c39bd04` | [`10.5281/zenodo.24368656`](https://doi.org/10.5281/zenodo.24368656) | **Intact (`PASSOU`)** |
| **Zenodo Concept DOI** | Resolves automatically to latest published release | [`10.5281/zenodo.24358066`](https://doi.org/10.5281/zenodo.24358066) | **Intact (`PASSOU`)** |

Per strict version-control rules (Phase G):
- **No new GitHub release (`v1.0.3`), tag, or Zenodo deposit has been published automatically.**
- **No future version DOI has been invented.** The manuscript cites the permanent Zenodo Concept DOI (`10.5281/zenodo.24358066`) alongside the existing audited Version DOIs (`v1.0.2`: `10.5281/zenodo.24368656`; `v1.0.1`: `10.5281/zenodo.24360150`; `v1.0.0`: `10.5281/zenodo.24358067`).
- All 14 consolidated supplementary tables (`Table_S0` through `Table_S9`, `MANIFEST_SHA256.tsv` SHA-256 `20bc50020084802e9fcdc9ca6af8618d0faed74ff61873157c8a2a8371e9cd11`) remain 100% unchanged from `enamel-relics-v1.0.2`.

---

## 2. Summary of Differences Relative to `enamel-relics-v1.0.2` (`db71156` / `59c1673`)

1. **Added [`scripts/enamel_relics/13_audit_final_neornithes_and_stats.py`](../../scripts/enamel_relics/13_audit_final_neornithes_and_stats.py) and [`runs/enamel_relics/v102/neornithes_detectability_audit.json`](../../runs/enamel_relics/v102/neornithes_detectability_audit.json):**
   - Computes exact cross-clade codon coverage breakdowns across all Level 4 avian alignments (`ENAM_exon8_aves`, `AMBN_exon3_aves`, `AMBN_exon4_aves`, `AMBN_exon9_aves`).
   - Executes an 87,800-trial (`878` codons $\times$ `100` replicates, `seed 20261124`) conditional detection-power simulation for a genuine crown-Neornithes ancestral PTC under neutral K80 evolution (`scale = 0.00, 0.65, 1.00, 1.25`) across `unmasked_full_878`, `empirical_mafft_mask`, and `empirical_dual_aligner_mask`.
   - Audits all 35 marginal empirical permutation *p*-value cells ($0.005 < p < 0.10$) across the 504-cell benchmark (`Table_S3c`).
2. **Updated [`scripts/enamel_relics/06_make_figures.py`](../../scripts/enamel_relics/06_make_figures.py) & [`submission_biorxiv/figures/fig3_lesions_and_phylogeny.png`](../../submission_biorxiv/figures/fig3_lesions_and_phylogeny.png):**
   - Updated the Figure 3B Crown Neornithes annotation box from `"0 congruent lesions"` to `"none detected (0/878 7-bird cov.; T1: 5–6/16 incong.; T2/T3: 9–11/16)"`.
3. **Updated [`submission_biorxiv/manuscript.typ`](../../submission_biorxiv/manuscript.typ), [`submission_biorxiv/manuscript.md`](../../submission_biorxiv/manuscript.md), and recompiled [`submission_biorxiv/manuscript_biorxiv.pdf`](../../submission_biorxiv/manuscript_biorxiv.pdf) (21 pages):**
   - **Phase A:** Replaced absolute `"zero crown Neornithes synapomorphies"` claims across the Abstract, Section 2.5 (#6), Figure 3 caption, Table 3, Section 3.4, and Section 3.5 (#2) with coverage- and power-qualified formulations (`0/878` 7-bird coverage; `16.3%` neutral uncontradicted ancestral-stop detection power; `40.8%` conversion into apparent cross-clade homoplasies by secondary stop-to-sense substitutions).
   - **Phase B:** Framed K80 simulation results as statistical compatibility conditional on the K80 model, branch lengths, topology, and empirical coverage mask (`cov_dual`), explicitly noting the partial circularity of conditioning on the observed mask and the upper-tail position of `T2`/`T3` parsimony incongruence ($p = 0.072–0.244$).
   - **Phase C:** Labeled `Strhab` `ENAM_Acar` ($p = 0.024$) as an unadjusted exploratory score across the 504-cell benchmark without multiplicity control, reinforced the single-locus linked-site caveat on Clopper–Pearson 95% CIs, and distinguished external phylogenomic support for `T1` (Sackton et al., 2019) from within-locus RELL topology tests on `ENAM_exon8_aves`.
   - **Phase D:** Calibrated scope of the single flanking globular *LPL* control (`arguing against uniform locus-wide alignment failure across NRG1–LPL` while noting that a single control gene cannot exclude broader lineage- or gene-specific rate variation) and removed rhetorical intensifiers (`"unequivocally"`).
4. **Expanded [`tests/test_enamel_relics_audit.py`](../../tests/test_enamel_relics_audit.py) (`13/13` test functions passed)** and created the 10 mandatory reports in [`audit/final/`](./).

---

## 3. SHA-256 Checksums of Updated & Created Artifacts

| File Path | Size (Bytes) | SHA-256 Checksum |
|---|---:|---|
| `submission_biorxiv/manuscript_biorxiv.pdf` | `2,743,865` | `9812a7c31f6fecf2366d0ede62af4e6766037ec8959a74f9d98f7c519c7a239d` |
| `submission_biorxiv/manuscript.typ` | `93,143` | `e39f4c71bb2d1db5c8a8548dae3c3784b6c0fb18fee1fa2e0f44966f172a8fef` |
| `submission_biorxiv/manuscript.md` | `87,518` | `726f50cb5b19a180c65fcce40ccb4c9b34ca58395076d5bee1f8dc40cbf89c01` |
| `submission_biorxiv/figures/fig3_lesions_and_phylogeny.png` | `676,017` | `4aa8437d9763c35001fe9105cbb2a72925ddace5677c6166ca8ea564a0ee547d` |
| `submission_biorxiv/figures/fig3_lesions_and_phylogeny.pdf` | `68,394` | `53a77007cf6fb8d3bd6628636a0baf3023ec3a9e7031dcd5d4f7baeed09f086c` |
| `submission_biorxiv/figures/fig3_lesions_and_phylogeny.svg` | `279,558` | `4942cf058ba4aef6e46cd5131ede9149c023f1249ec6117b6447d21e58a1f974` |
| `submission_biorxiv/supplementary/MANIFEST_SHA256.tsv` | `2,046` | `20bc50020084802e9fcdc9ca6af8618d0faed74ff61873157c8a2a8371e9cd11` |
| `scripts/enamel_relics/06_make_figures.py` | `47,180` | `c7c4444745be5f26c05fbf986bc1d1cf8fa47f5bd9c591ede41c6434b1b689be` |
| `scripts/enamel_relics/13_audit_final_neornithes_and_stats.py` | `13,339` | `da78465f5e60b7b147584584e149a47ce87a5641f6bcb243e99c067ad3bf1c62` |
| `runs/enamel_relics/v102/neornithes_detectability_audit.json` | `19,641` | `1ac6441ebb6b91a65542e8d1c0c166a4682dc58479fbe137a1bc4d281f6f4796` |
| `tests/test_enamel_relics_audit.py` | `32,367` | `a1d15629a446ed9ed86e542a6de40e5711165f05fcc0a630d8471ed5b20673e6` |
