# Human Approval Checklist for `enamel-relics-v1.0.2`

**Branch:** `audit/biorxiv-final-v1.0.2`  
**Baseline Releases (Preserved Untouched):** `enamel-relics-v1.0.1` (`10.5281/zenodo.24360150`) & `enamel-relics-v1.0.0` (`10.5281/zenodo.24358067`)  
**Concept DOI:** `10.5281/zenodo.24358066`

---

## 1. Mandatory Decisions Depending Exclusively on the Human Author

- [ ] **1.1. Institutional Affiliation & Pre-Submission Enquiries (`PENDENTE DE DECISÃO HUMANA`):**
  - Review the two draft emails in [`audit/v102/biorxiv_readiness_report.md`](biorxiv_readiness_report.md) (Section 6.1 to Centro Universitário Senac; Section 6.2 to `bioRxiv@openrxiv.org`).
  - Decide whether to submit with **Option A** (`Independent Researcher, São Paulo, Brazil` + factual student *Author Note*, currently compiled in `submission_biorxiv/manuscript_biorxiv.pdf`) after confirming eligibility with bioRxiv, OR update to **Option B** (`Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil` + non-endorsement disclaimer) **only after receiving written authorization from Senac**.
- [ ] **1.2. Verify Accuracy of the Revised Generative AI Declaration:**
  - Read the updated section *Declaration of Generative AI and AI-Assisted Technologies* in `submission_biorxiv/manuscript_biorxiv.pdf` (page 17) and confirm that it accurately reflects your personal contribution, direction, review, and full accountability without overstating manual alignment checks.
- [ ] **1.3. Scientific & Visual Sign-Off on `v1.0.2` Corrections:**
  - Confirm the four Phase A–D revisions in `submission_biorxiv/manuscript_biorxiv.pdf`:
    1. **Phase A (Variant Nomenclature):** Separation of *Alligator* reference protein coordinates (`XP_019346377.1:p.Phe601Ter`, `p.Trp1003Ter`) from reconstructed avian ancestral codons (`TTT [Phe] → TAT [Tyr, P=0.899] → TAA [Ter]` at codon 601; `TGG [Trp, P=0.168] / site-independent ML TGC [Cys, P=0.343] → TGA [Ter]` at codon 1003) in Section 2.5, Table 3, Figure 3B, and `Table_S5a`.
    2. **Phase B (Synapomorphy Evidence Tiers):** Qualification of `AMBN` Exon 4 16-nt deletion (Palaeognathae) and 1-nt insertion + `AG→AC` splice (Chemyd+Chrpic) as shared clade lesions where sister-clade coverage is unrecovered, and explicit separation of ML ASR inference from unreversed observation at `ENAM` codon 1022 (`p.Ser1022Ter`, `Anapla` uncovered).
    3. **Phase C (Topological Sensitivity `T1`/`T2`/`T3`):** Inclusion of `Table_S9` showing parsimony incongruence across `T1` (`6/16 = 37.5%`), `T2` (`10/16 = 62.5%`), and `T3` (`11/16 = 68.8%`) (reversible ML ASR: `31.2%`, `56.2%`, `62.5%`), formal IQ-TREE 3 KH/SH/wSH/ELW/AU rejection of `T2`/`T3` ($p < 0.01$), external phylogenomic grounding of `T1`, and single-locus linkage caveat on binomial intervals.
    4. **Phase D (K80 Neutral Simulation Coverage Mask & Sensitivity Grid):** Inclusion of `Table_S8` showing that applying the empirical per-species coverage mask reduces simulated shared PTCs from `26–29` (unmasked 878 codons) to `16.0–19.0` (matching `16` observed, two-sided $p = 0.54–1.00$) while maintaining `40.0–44.4%` incongruence and median `0` crown-Neornithes synapomorphies.
- [ ] **1.4. Pull Request Approval, Merge & `enamel-relics-v1.0.2` Release Minting:**
  - If approved, merge the PR (`audit/biorxiv-final-v1.0.2` → `main`) and publish GitHub Release `enamel-relics-v1.0.2`.
  - Once Zenodo automatically mints the new `v1.0.2` Version DOI under Concept DOI `10.5281/zenodo.24358066`, optionally record the new `v1.0.2` Version DOI in `Data and Code Availability` (currently citing Concept DOI `10.5281/zenodo.24358066` plus historical `v1.0.1` `10.5281/zenodo.24360150` and `v1.0.0` `10.5281/zenodo.24358067`).
- [ ] **1.5. Manual Submission on the bioRxiv Portal:**
  - Upload `submission_biorxiv/manuscript_biorxiv.pdf` and `submission_biorxiv/supplementary/` (`Table_S0`–`Table_S9` + `MANIFEST_SHA256.tsv`), paste the ≤250-word Abstract, select **Evolutionary Biology** (secondary: **Genomics**), license **CC BY 4.0**, and complete author accountability sign-off.
