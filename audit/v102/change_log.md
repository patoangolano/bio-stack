# Change Log — From `enamel-relics-v1.0.1` (`13b34a6` / `d18be9d`) to Candidate `enamel-relics-v1.0.2` (`audit/biorxiv-final-v1.0.2`)

**Baseline Release:** `enamel-relics-v1.0.1` (tag `13b34a6f77a7b7fdd56d95ee1eac4549d6be9370`, Zenodo Version DOI `10.5281/zenodo.24360150`, followed on `main` by PDF DOI insertion commit `d18be9d9990236de0f5eb9ecd76565d1c4ed7bca`)  
**Target Branch:** `audit/biorxiv-final-v1.0.2`  
**Audit Date:** 2026-11-19  
**Historical Tags Preserved Untouched:** `enamel-relics-v1.0.0` (`b9ff25e654007ef26f5ce01ce702f68c214c90b2`), `enamel-relics-v1.0.1` (`13b34a6f77a7b7fdd56d95ee1eac4549d6be9370`)

---

## 1. Computational Scripts Added or Modified

1. **Added [`scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py`](../../scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py):**
   - **Phase A:** Programmatically verifies all 23 shared coding mutations in `Table_S5a` against *Alligator mississippiensis* reference proteins `XP_019346377.1` (*ENAM*, 1,078 aa) and `XP_019346386.1` (*AMBN*, 390 aa) and `IQ-TREE 3.1.4` (`GTR+G --ancestral`) ancestral state reconstructions; enriches `Table_S5a` with 12 new columns (`ref_aa_3letter`, `croc_consensus_aa_3letter`, `neornithes_ml_aa_3letter`, `neornithes_ml_p_best`, `neognathae_ml_p_stop`, `inferred_substitution_path`, `fitch_bin_steps_T2`, `fitch_bin_steps_T3`, `min_gains_no_rev_T2`, `min_gains_no_rev_T3`, `synapomorphy_evidence_tier`).
   - **Phase B:** Executes ultra-sensitive `blastn -task blastn -word_size 7 -evalue 10 -dust no` and `tblastn -word_size 2 -evalue 50 -seg no -comp_based_stats 0` across all 18 sauropsid `NRG1–LPL` intervals for *AMBN* Exon 4 (`runs/enamel_relics/v102/ambn_exon4_coverage_audit/`), recovering intact toothed squamate outgroups (`Anocar`, `Pytbiv`, `Podmur`) and a single diverged Neognathae relic in `Anapla` (`interval_Apla.fna:135,001–135,234`, `+`, 57.2 bits) and aligning it with `MAFFT L-INS-i` and `MUSCLE 5`.
   - **Phase C:** Runs formal 10,000-replicate RELL topology tests (`-zb 10000 -zw -au`: KH, SH, wKH, wSH, ELW, AU) and ML ancestral state reconstructions across topologies `T1`, `T2`, and `T3` in `IQ-TREE 3.1.4` on both MAFFT and MUSCLE alignments (`runs/enamel_relics/v102/iqtree_tests/`), producing `Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv`.
   - **Phase D:** Executes a 48-scenario × 1,000-replicate (`48,000` simulations, `numpy.random.default_rng(20261123)`) neutral K80 stop-gain sensitivity benchmark comparing `unmasked_full_878`, `empirical_mafft_mask`, and `empirical_dual_aligner_mask` across $\kappa \in \{2.0, 4.0, 6.0\}$, branch scales $0.65\times–1.25\times$, $t_{\text{stem}} \in \{0.00, 0.02, 0.05\}$, root sequences (`Allmis_ref` vs. `Neornithes_ML`), and topologies `T1`–`T3`, producing `Table_S8_k80_simulation_sensitivity_v102.tsv`.
2. **Updated [`scripts/enamel_relics/11_consolidate_v101_tables.py`](../../scripts/enamel_relics/11_consolidate_v101_tables.py):**
   - Incorporates `Table_S5a_level4_shared_mutations_v102.tsv`, `Table_S8_k80_simulation_sensitivity_v102.tsv`, and `Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv` into the consolidated supplementary package (`14` TSV tables + `MANIFEST_SHA256.tsv`), enforcing strict `LF` (`\n`) line endings across both `runs/enamel_relics/supplementary/` and `submission_biorxiv/supplementary/`.
3. **Updated [`scripts/enamel_relics/06_make_figures.py`](../../scripts/enamel_relics/06_make_figures.py):**
   - Updated Figure 3B panel title (`shared coding lesions`), crown Neornithes summary badge (`T1: 5–6/16 incong.; T2/T3: 9–11/16`), Palaeognathae badge (`AMBN Ex4 1 fs (16-nt del*)`), Durocryptodira badge (`Chemyd + Chrpic ...*`), and Psittacopasserae badge (`ref p.Trp1003Ter; avian TGG/TGC→TGA`).
   - Updated Figure 3C x-axis labels to show both reconstructed crown-Neornithes codon and 3-letter amino acid (`Cod 311\n[CAA: Gln]`, `Cod 439\n[TCA: Ser]`, `Cod 562\n[TGG: Trp]`, `Cod 615\n[TGG: Trp]`, `Cod 684\n[CAA: Gln]`, `Cod 1022\n[TCG: Ser]`).

---

## 2. Supplementary Tables & Figures Updated

- **Updated:** [`submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv`](../../submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv) (39 columns across 23 variants).
- **Added:** [`submission_biorxiv/supplementary/Table_S8_k80_simulation_sensitivity_v102.tsv`](../../submission_biorxiv/supplementary/Table_S8_k80_simulation_sensitivity_v102.tsv) (48 simulation scenarios × 27 metrics).
- **Added:** [`submission_biorxiv/supplementary/Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv`](../../submission_biorxiv/supplementary/Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv) (formal IQ-TREE 3 KH/SH/wSH/ELW/AU tests and `T1`–`T3` parsimony/ML homoplasy rates).
- **Updated:** [`submission_biorxiv/supplementary/MANIFEST_SHA256.tsv`](../../submission_biorxiv/supplementary/MANIFEST_SHA256.tsv) (14 TSV hashes + byte sizes, strict LF endings).
- **Re-rendered:** [`submission_biorxiv/figures/fig3_lesions_and_phylogeny.png`](../../submission_biorxiv/figures/fig3_lesions_and_phylogeny.png), `.pdf`, `.svg`.

---

## 3. Manuscript & Metadata Files Updated

- **Updated [`submission_biorxiv/manuscript.typ`](../../submission_biorxiv/manuscript.typ), [`submission_biorxiv/manuscript.md`](../../submission_biorxiv/manuscript.md), and recompiled [`submission_biorxiv/manuscript_biorxiv.pdf`](../../submission_biorxiv/manuscript_biorxiv.pdf) (21 pages):**
  - Abstract (247 words ≤ 250), Introduction (#4), Section 2.5, Figure 3 caption, Table 3, Section 3.1, Section 3.4, Section 3.5 (#2), Section 4.4, Section 4.5, Data and Code Availability, and Declaration of Generative AI and AI-Assisted Technologies.
- **Updated [`submission_biorxiv/README_SUBMISSION.md`](../../submission_biorxiv/README_SUBMISSION.md) and [`CITATION.cff`](../../CITATION.cff).**
- **Expanded [`tests/test_enamel_relics_audit.py`](../../tests/test_enamel_relics_audit.py)** with v1.0.2 tests for codon-to-AA translation consistency, outgroup/sister coverage tiers, `T1`/`T2`/`T3` homoplasy counts, coverage-masked vs. unmasked K80 simulations, 14-table SHA-256 manifest integrity, presence of all 9 `audit/v102/` deliverables, and immutability of tags `enamel-relics-v1.0.0` and `enamel-relics-v1.0.1`.
