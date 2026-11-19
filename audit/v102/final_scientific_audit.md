# Final Scientific Audit Report — Candidate Release `enamel-relics-v1.0.2`

**Manuscript Title:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles*  
**Repository:** `https://github.com/patoangolano/bio-stack`  
**Audit Branch:** `audit/biorxiv-final-v1.0.2` (branched from `enamel-relics-v1.0.1` / `d18be9d`)  
**Historical Releases Preserved Untouched:**
- `enamel-relics-v1.0.0` (tag `b9ff25e654007ef26f5ce01ce702f68c214c90b2`, Zenodo Version DOI [`10.5281/zenodo.24358067`](https://doi.org/10.5281/zenodo.24358067))
- `enamel-relics-v1.0.1` (tag `13b34a6f77a7b7fdd56d95ee1eac4549d6be9370`, Zenodo Version DOI [`10.5281/zenodo.24360150`](https://doi.org/10.5281/zenodo.24360150))
- Zenodo Concept DOI: [`10.5281/zenodo.24358066`](https://doi.org/10.5281/zenodo.24358066)

---

## 1. Executive Summary of Audit Classifications (`Phases A–H`)

| Phase | Scope | Classification | Readiness Status | Detailed Report |
|---|---|---|---|---|
| **Baseline (`v1.0.1`)** | 18-genome annotation audit (`RefSeq 106` + `Ensembl 116`), chicken chrZ relics (`ΨENAM` Ex8/Ex6, `ΨAMBN` Ex3), kakapo `ΨAMTN`, 4 Gnomon `LOW QUALITY PROTEIN` models (`Table_S7`), 1,000-permutation IDP benchmark (`Table_S3a–c`), `miniprot` coding lesions (`Table_S4`), 49/49 verified references | `CONFIRMADO` | `PASSOU` | [`audit/reproducibility_audit.md`](../reproducibility_audit.md), [`audit/statistical_review.md`](../statistical_review.md), [`audit/bibliography_report.md`](../bibliography_report.md) |
| **Phase A** | Molecular variant nomenclature & *Alligator* reference vs. reconstructed crown-Neornithes ancestral codons across all 23 shared variants (`codons 601, 1003, 311, 439, 559, 844, 1022`) | `CORRIGIDO` | `PASSOU` | [`audit/v102/variant_nomenclature.md`](variant_nomenclature.md) |
| **Phase B** | Synapomorphy validation & outgroup/sister-clade coverage grading (`AMBN` Exon 4 Palaeognathae 16-nt del + ultra-sensitive `Anapla` Exon 4 recovery; `Chemyd`+`Chrpic` lesions with `Gopevg`/`Pelsin` unrecovered; `ENAM` codon 1022 `Anapla` missing coverage + parsimony equiparsimony vs. ML ASR) | `CORRIGIDO` | `PASSOU` | [`audit/v102/phylogenetic_robustness.md`](phylogenetic_robustness.md) |
| **Phase C** | Palaeognathae topological sensitivity (`T1` Novaeratitae vs. `T2` `Strcam`+`Dronov` vs. `T3` `Strcam`+`Aptrow`) & 10,000-replicate RELL KH/SH/wSH/ELW/AU topology tests in `IQ-TREE 3.1.4` (`Table_S9`) | `CORRIGIDO` | `PASSOU` | [`audit/v102/phylogenetic_robustness.md`](phylogenetic_robustness.md) |
| **Phase D** | Forensic audit of neutral K80 stop-gain simulation (`10_audit_mutations_and_asr.py` unmasked 878-codon bug vs. `12_audit_v102_variants_phylo_sim.py` empirical dual-aligner per-species coverage mask across 48 scenarios × 1,000 replicates; `Table_S8`) | `CORRIGIDO` | `PASSOU` | [`audit/v102/neutral_simulation_report.md`](neutral_simulation_report.md) |
| **Phase E** | Scientific text precision, linked-site binomial CI caveat, Figure 3B–C & Table 3 alignment, Concept DOI vs. Version DOI clarity, and transparent AI Contribution Declaration | `CORRIGIDO` | `PASSOU` | [`audit/v102/change_log.md`](change_log.md) |
| **Phase F** | Author affiliation policy (Option A vs. Option B) and pre-submission human verification steps for openRxiv screening | `PENDENTE DE DECISÃO HUMANA` | `PENDENTE` | [`audit/v102/biorxiv_readiness_report.md`](biorxiv_readiness_report.md), [`audit/v102/human_approval_checklist.md`](human_approval_checklist.md) |
| **Phase G** | Automated acceptance test suite (`tests/test_enamel_relics_audit.py`) and SHA-256 manifest verification across all 14 supplementary tables (`Table_S0`–`Table_S9`) | `CORRIGIDO` | `PASSOU` | [`audit/v102/test_results.md`](test_results.md) |

---

## 2. Synthesis of Key Corrections in `v1.0.2`

1. **Reference Protein Coordinate vs. Crown-Avian Ancestral State (`Phase A`, `Table_S5a`):**
   - Explicitly separated the *Alligator mississippiensis* reference residue (`XP_019346377.1:p.Xxx###Ter`) from the reconstructed crown-Neornithes ancestral codon/amino acid across all 23 shared lesions.
   - At **`ENAM` codon 601 (`XP_019346377.1:p.Phe601Ter`)**, all three crocodylians carry `TTT` (**Phe**), crown Neornithes and Neognathae ML reconstruct `TAT` (**Tyr**, $P = 0.899$ in MAFFT / $0.870$ dual-aligner mean; retained as `TAT` in `Galgal`/`Strhab` and `TAC` in `Anapla`), and crown Palaeognathae underwent the 1-step transversion `TAT (Tyr) → TAA (Ter)` ($P_{\text{stop}} = 0.994$).
   - At **`ENAM` codon 1003 (`XP_019346377.1:p.Trp1003Ter`)**, all three crocodylians carry `TGG` (**Trp**), Palaeognathae carry `GGG`/`GGC` (**Gly**), `Galgal` carries `TAC` (**Tyr**), `Anapla` is uncovered (`---`), and Psittacopasserae (`Strhab` + `Taegut`) carry `TGA` (**Ter**, ML $P_{\text{stop}} = 0.895$). Under site-independent nucleotide `GTR+G` ML ASR, multiplying the three per-position base marginals (`T × G × C`) at crown Neornithes yields a composite argmax of `TGC` (**Cys**, $P = 0.343$) slightly ahead of `TGA` (**Ter**, $P = 0.291$) and `TGG` (**Trp**, $P = 0.168$), even though no extant bird carries `TGC` and 3-step nucleotide Fitch parsimony equally allows a direct `TGG (Trp) → TGA (Ter)` transition or an intermediate `TGC (Cys) → TGA (Ter)` transversion.

2. **Outgroup & Sister-Clade Coverage Grading (`Phase B`, `Table_S5a`):**
   - Added `synapomorphy_evidence_tier` to `Table_S5a` to distinguish **13 outgroup- and sister-verified synapomorphies** (`OUTGROUP_AND_SISTER_VERIFIED_SYNAPOMORPHY`) from **4 shared clade/subclade lesions whose sister lineages lack baseline exon recovery**:
     - Palaeognathae *AMBN* Exon 4 16-nt frameshift deletion at codon 109 (`SHARED_CLADE_LESION_SISTER_UNRECOVERED`): verified absent in all 3 crocodylians and 3 toothed squamates (`Anocar`, `Pytbiv`, `Podmur`). Baseline `dc-megablast` recovers zero Neognathae; ultra-sensitive `-word_size 7` `blastn` / `-word_size 2` `tblastn` recovers only `Anapla` (`57.2 bits`, $E = 9.25\times 10^{-10}$, canonical `AG..GT`), which preserves intact ancestral length at upstream `nt 61–76` but carries a shorter 10-nt deletion at codon 109.
     - Testudine shared lesions in `Chemyd` + `Chrpic` (*AMBN* Exon 4 1-nt insertion at codon 123, flanked-MSA `AG→AC` splice acceptor, and *ENAM* Exon 8 `p.Gln753Ter` & `p.Ala859Ter`; `SHARED_SUBCLADE_LESION_SISTER_UNRECOVERED`): unrecoverable in `Gopevg` (within Durocryptodira) and `Pelsin` (Trionychia) even under `-word_size 7` `blastn` and `-word_size 2` `tblastn`, and therefore reported as minimum inactivation bounds for the sampled `Chemyd` + `Chrpic` subclade (Chelonioidea + Testudinoidea).
   - At **`ENAM` codon 1022 (`p.Ser1022Ter`)**, `Anapla` is uncovered (`---`), `Galgal` and `Strhab` carry `TAG`, and `Taegut` carries sense `AAA` (**Lys**, 2 steps from `TAG`). Binary Fitch parsimony is strictly equiparsimonious (2 steps: two independent gains vs. one Neognathae gain + one `Taegut` reversal), whereas reversible `GTR+G` ML ASR assigns $P_{\text{stop}} = 0.926$ (`T1`) / $0.908$ (`T2`) / $0.909$ (`T3`) at crown Neognathae. Both interpretations (`6/16 = 37.5%` under irreversible parsimony vs. `5/16 = 31.2%` under reversible ML ASR) are presented transparently.

3. **Topological Sensitivity (`T1`–`T3`) & Formal IQ-TREE 3 Topology Tests (`Phase C`, `Table_S9`):**
   - Formal 10,000-replicate RELL topology tests (`-zb 10000 -zw -au`) in `IQ-TREE 3.1.4` confirm that *ENAM* Exon 8 independently favors `T1` (Novaeratitae) and rejects `T2` ($\Delta\ln L = -22.96$ MAFFT / $-22.10$ MUSCLE; $p_{\text{KH}} \le 0.0040, p_{\text{SH}} \le 0.0040, p_{\text{WSH}} \le 0.0045, c_{\text{ELW}} \le 0.0016, p_{\text{AU}} \le 0.0015$) and `T3` ($\Delta\ln L = -23.18$ / $-22.15$; $p_{\text{KH}} \le 0.0033, p_{\text{SH}} \le 0.0033, p_{\text{WSH}} \le 0.0096, c_{\text{ELW}} \le 0.0012, p_{\text{AU}} \le 0.0017$).
   - Across topologies, avian *ENAM* Exon 8 PTC incongruence is **`6/16` (`37.5%`, nominal linked-site 95% CI `15.2–64.6%`)** under irreversible parsimony and **`5/16` (`31.2%`, `11.0–58.7%`)** under reversible ML ASR on **`T1`**, rising to **`10/16` (`62.5%`) / `9/16` (`56.2%`) on `T2`** and **`11/16` (`68.8%`) / `10/16` (`62.5%`) on `T3`**.

4. **Forensic Resolution of the Neutral K80 Stop-Gain Simulation (`Phase D`, `Table_S8`):**
   - Identified that `10_audit_mutations_and_asr.py` (`v1.0.1`) had a comment stating it masked simulated birds to their empirical coverage but omitted the mask in code, simulating `100%` coverage (`878/878` codons in all 7 birds) and predicting `25.0–29.0` shared PTCs (observed `16` at two-sided $p = 0.004–0.052$).
   - In the empirical `ENAM_exon8_aves` alignment, dual-aligner codon coverage is `22.7%` in `Taegut` (`199/878` codons), `58.7–60.4%` in `Anapla`/`Galgal`, and `79.0–82.5%` in Palaeognathae/`Strhab` (`0/878` codons covered across all 7 birds).
   - Applying the exact empirical dual-aligner per-species coverage mask (`Table_S8`, seed `20261123`, 48 scenarios × 1,000 replicates) predicts a baseline median of **`16.0–19.0` shared PTCs** (95% range `10.0–27.0`; observed `16` at `48.8th–23.4th` percentile, two-sided empirical $p = 0.538–1.000$), **`40.0–44.4%` unconditional incongruence** (two-sided $p = 0.618–0.848$ vs. observed `37.5%`), **`37.5%` conditional median incongruence at $N_{\text{shared}} = 16$**, and **`0.0` crown-Neornithes PTC synapomorphies** ($P(N_{\text{syn}} = 0) = 50.3–94.6\%$).
