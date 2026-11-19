# Phase D Audit Report — Forensic Validation & Coverage-Masked Sensitivity Grid of the Neutral K80 Stop-Gain Simulation (`v1.0.2`)

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-final-v1.0.2`  
**Primary Script:** [`scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py`](../../scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py)  
**Consolidated Table:** [`submission_biorxiv/supplementary/Table_S8_k80_simulation_sensitivity_v102.tsv`](../../submission_biorxiv/supplementary/Table_S8_k80_simulation_sensitivity_v102.tsv)

---

## 1. Executive Summary & Root-Cause Diagnosis

In `enamel-relics-v1.0.1`, Section 2.5 (#6) and Section 3.4 reported that simulating 1,000 replicates of neutral K80 ($\kappa = 4$) evolution of the 878-codon (2,634-bp) *ENAM* Exon 8 along the empirical `IQ-TREE 3` 7-bird branch lengths (`seed 20261122`, `10_audit_mutations_and_asr.py`) predicted a median of **`26–29` shared PTCs** and **`48.1–51.8%` incongruence** (`95%` range `29.2–70.8%`), concluding that the observed data (`16` shared PTCs, `31.2–37.5%` incongruent) "match neutral paleogenomic expectations directly."

Forensic code inspection of [`scripts/enamel_relics/10_audit_mutations_and_asr.py`](../../scripts/enamel_relics/10_audit_mutations_and_asr.py) uncovered the exact reason why the `v1.0.1` simulation predicted `26–29` shared PTCs against the `16` observed empirically:

1. **Omitted Empirical Per-Species Coverage Mask in `10_audit_mutations_and_asr.py`:**
   - Line 436 of `10_audit_mutations_and_asr.py` contained the comment:
     `# Mask each simulated species to its ACTUAL empirical coverage mask in ENAM_exon8_aves!`
   - **However, lines 450–469 never applied the coverage mask** and instead evaluated every simulated bird across all `878/878` codons (`100%` coverage in all 7 birds) and tested congruence against full 7-bird clade sets (`Palaeognathae = {Scam, Arow, Dnov}`, `Neognathae = {Ggal, Apla, Shab, Tgut}`, etc.) without missing-taxon parsimony (`None` for uncovered tips).
2. **Actual Empirical Codon Coverage in `ENAM_exon8_aves`:**
   - Due to partial exon decay and terminal truncation across >100 Myr of evolution, the empirical fraction of the 878 *Alligator* Exon 8 codons covered by non-gap triplets in each bird is far below 100%:

| Species (6-letter ID) | MAFFT Coverage (`codons / 878`, %) | Dual-Aligner Consensus Coverage (`codons / 878`, %) |
|---|---|---|
| `Strcam` (*Struthio camelus*) | `706 / 878` (`80.4%`) | `694 / 878` (`79.0%`) |
| `Aptrow` (*Apteryx rowi*) | `732 / 878` (`83.4%`) | `724 / 878` (`82.5%`) |
| `Dronov` (*Dromaius novaehollandiae*) | `707 / 878` (`80.5%`) | `699 / 878` (`79.6%`) |
| `Galgal` (*Gallus gallus*) | `637 / 878` (`72.6%`) | `530 / 878` (`60.4%`) |
| `Anapla` (*Anas platyrhynchos*) | `597 / 878` (`68.0%`) | `515 / 878` (`58.7%`) |
| `Strhab` (*Strigops habroptila*) | `714 / 878` (`81.3%`) | `707 / 878` (`80.5%`) |
| `Taegut` (*Taeniopygia guttata*) | `223 / 878` (`25.4%`) | `199 / 878` (`22.7%`) |
| **Simultaneously covered in all 7 birds** | **`0 / 878` (`0.0%`)** | **`0 / 878` (`0.0%`)** |

Because `Taegut` covers only `22.7%` of Exon 8 and `Anapla`/`Galgal` cover `58.7–60.4%`, an unmasked simulation assuming `100%` coverage in all 7 birds artificially inflates the number of opportunities for two birds (especially `Strhab` + `Taegut` or `Galgal` + `Anapla`) to share a stop codon by ~1.5-fold!

---

## 2. Quantitative Comparison: Unmasked (`v1.0.1`) vs. Empirically Coverage-Masked (`v1.0.2`) Neutral Model (`Table_S8`)

We executed a 48-scenario sensitivity grid (`1,000` replicates per scenario = `48,000` simulations, `numpy.random.default_rng(20261123)`, logged in `runs/enamel_relics/v102/Table_S8_k80_simulation_sensitivity_v102.tsv`) varying:
- **Coverage mask (`mask_mode`):** `unmasked_full_878` (the `v1.0.1` unmasked 878-codon grid), `empirical_mafft_mask` (single-aligner MAFFT coverage), and `empirical_dual_aligner_mask` (requiring valid codon coverage in both MAFFT and MUSCLE, matching our empirical dual-aligner consensus filter), using exact missing-taxon Fitch parsimony (`None` for uncovered tips).
- **Root ancestral sequence (`anc_seq_source`):** `Allmis_ref` (`Alligator` Exon 8 sequence) vs. `Neornithes_ML` (`IQ-TREE 3` ML reconstructed crown-Neornithes ancestral sequence, with stops replaced by best sense codons).
- **Transition/transversion parameter ($\kappa$):** `2.0`, `4.0` (baseline), `6.0`.
- **Branch-length scaling factor (`branch_scale`):** `0.65×`, `0.80×`, `1.00×` (empirical IQ-TREE 3 branch lengths), `1.25×`.
- **Pre-crown Neornithes stem branch (`t_stem`):** `0.00`, `0.02`, `0.05` substitutions/site.
- **Palaeognathae topology (`topology`):** `T1` (Novaeratitae), `T2` (`Strcam` + `Dronov`), `T3` (`Strcam` + `Aptrow`).

### 2.1. Baseline Comparison ($\kappa = 4.0$, Branch Scale $1.0\times$, Topology `T1`)

| Root Sequence & Coverage Mask | $t_{\text{stem}}$ | Shared PTCs Median (IQR; 95% Range) | Observed `16` Percentile Rank | Two-Sided Empirical $p$ ($N_{\text{shared}}=16$) | Unconditional Incongruence Median (95% Range) | Two-Sided $p$ (vs. `37.5%`) | Conditional Incongruence Median at $N_{\text{shared}}=16$ ($n$ reps) | Conditional Incongruence Median at $N_{\text{shared}}\in[14,18]$ ($n$ reps) | Crown-Neornithes Synapomorphies Median ($P(N_{\text{syn}}=0)$) |
|---|---:|---|---:|---:|---|---:|---|---|---|
| `Allmis_ref` · `unmasked_full_878` | `0.00` | `25.0` (`22–28`; `17–35`) | `1.0%` | `0.032` | `47.8%` (`28.6–66.7%`) | `0.340` | `37.5%` ($n=11$) | `47.1%` ($n=56$) | `0.0` (`98.1%`) |
| `Allmis_ref` · `unmasked_full_878` | `0.02` | `27.0` (`24–30`; `18–37`) | `0.5%` | `0.014` | `50.0%` (`32.1–68.6%`) | `0.170` | `50.0%` ($n=3$) | `50.0%` ($n=24$) | `0.0` (`81.5%`) |
| `Allmis_ref` · `unmasked_full_878` | `0.05` | `29.0` (`26–33`; `20–39`) | `0.2%` | `0.004` | `52.2%` (`34.5–70.3%`) | `0.084` | — ($n=0$) | `51.5%` ($n=10$) | `0.0` (`62.5%`) |
| `Neornithes_ML` · `unmasked_full_878` | `0.00` | `25.0` (`22–28`; `16–34`) | `2.0%` | `0.052` | `50.0%` (`29.4–69.2%`) | `0.266` | `56.2%` ($n=11$) | `51.5%` ($n=68$) | `0.0` (`97.4%`) |
| `Neornithes_ML` · `unmasked_full_878` | `0.02` | `26.0` (`23–29`; `17–36`) | `1.5%` | `0.040` | `50.0%` (`32.0–69.7%`) | `0.156` | `53.1%` ($n=10$) | `50.0%` ($n=47$) | `0.0` (`78.2%`) |
| `Neornithes_ML` · `unmasked_full_878` | `0.05` | `28.0` (`25–32`; `20–38`) | `0.2%` | `0.008` | `53.6%` (`34.8–70.8%`) | `0.104` | `53.1%` ($n=4$) | `58.8%` ($n=13$) | `0.0` (`61.1%`) |
| `Neornithes_ML` · `empirical_mafft_mask` | `0.00` | `18.0` (`15–21`; `10–26`) | `32.1%` | `0.732` | `41.2%` (`16.7–64.7%`) | `0.782` | `43.8%` ($n=89$) | `41.2%` ($n=417$) | `0.0` (`93.0%`) |
| `Neornithes_ML` · `empirical_mafft_mask` | `0.02` | `19.0` (`16–22`; `12–27`) | `20.9%` | `0.504` | `43.5%` (`21.1–66.7%`) | `0.676` | `43.8%` ($n=87$) | `43.8%` ($n=370$) | `0.0` (`69.5%`) |
| **`Neornithes_ML` · `empirical_dual_aligner_mask`** | **`0.00`** | **`16.0` (`14–19`; `10–24`)** | **`48.8%`** | **`1.000`** | **`40.0%` (`16.7–65.0%`)** | **`0.848`** | **`37.5%` ($n=111$)** | **`40.0%` ($n=489$)** | **`0.0` (`94.6%`)** |
| **`Neornithes_ML` · `empirical_dual_aligner_mask`** | **`0.02`** | **`18.0` (`15–20.2`; `10–26`)** | **`33.8%`** | **`0.778`** | **`41.7%` (`17.4–66.7%`)** | **`0.762`** | **`43.8%` ($n=102$)** | **`42.9%` ($n=443$)** | **`0.0` (`73.0%`)** |
| **`Neornithes_ML` · `empirical_dual_aligner_mask`** | **`0.05`** | **`19.0` (`16–22`; `12–27`)** | **`23.4%`** | **`0.538`** | **`44.4%` (`21.4–68.4%`)** | **`0.618`** | **`43.8%` ($n=71$)** | **`44.4%` ($n=385$)** | **`0.0` (`50.3%`)** |
| `Allmis_ref` · `empirical_dual_aligner_mask` | `0.00` | `16.0` (`14–19`; `10–25`) | `45.2%` | `1.000` | `40.0%` (`16.7–64.7%`) | `0.860` | `43.8%` ($n=102$) | `40.0%` ($n=477$) | `0.0` (`94.7%`) |
| `Allmis_ref` · `empirical_dual_aligner_mask` | `0.02` | `18.0` (`15–21`; `11–26`) | `32.1%` | `0.730` | `41.2%` (`20.0–65.0%`) | `0.786` | `43.8%` ($n=88$) | `41.2%` ($n=439$) | `0.0` (`72.0%`) |

---

## 3. Sensitivity Across $\kappa$, Branch-Length Scaling, and Topologies `T1`–`T3`

1. **Effect of $\kappa \in \{2.0, 4.0, 6.0\}$ (at `empirical_dual_aligner_mask`, scale $1.0\times$, $t_{\text{stem}}=0.00–0.02$):**
   - $\kappa = 2.0$: median **`17.0–19.0`** shared PTCs ($p = 0.584–0.888$), incongruence median **`38.1–40.9%`** ($p = 0.788–0.984$).
   - $\kappa = 4.0$: median **`16.0–18.0`** shared PTCs ($p = 0.778–1.000$), incongruence median **`40.0–41.7%`** ($p = 0.762–0.848$).
   - $\kappa = 6.0$: median **`16.0–17.0`** shared PTCs ($p = 0.798–1.000$), incongruence median **`43.8–44.4%`** ($p = 0.570–0.642$).
2. **Effect of Branch-Length Scaling ($0.65\times$ to $1.25\times$, reflecting delayed subclade inactivation vs. accelerated pseudogene drift):**
   - $0.65\times$ branch scale (representing functional constraint during early crown divergence prior to subclade inactivation): predicts median **`11.0–13.0`** shared PTCs (`95%` range `6.0–21.0`, two-sided $p = 0.218–0.532$) and **`30.8–36.4%`** incongruence ($p = 0.702–0.930$).
   - $0.80\times$ branch scale: predicts median **`14.0–16.0`** shared PTCs (`95%` range `7.0–24.0`, two-sided $p = 0.564–1.000$) and **`35.0–40.0%`** incongruence ($p = 0.838–0.968$).
   - $1.25\times$ branch scale: predicts median **`19.0–21.0`** shared PTCs (`95%` range `11.0–30.0`, two-sided $p = 0.264–0.544$) and **`41.7–48.0%`** incongruence ($p = 0.320–0.726$).
3. **Effect of Topologies `T2` and `T3` (`empirical_dual_aligner_mask`, $\kappa = 4.0$, scale $0.80\times–1.00\times$):**
   - On `T2` (`Strcam` + `Dronov`), simulated neutral incongruence median is **`43.8–48.0%`** (`95%` range `16.7–72.2%`, two-sided $p = 0.162–0.244$ vs. observed `T2` incongruence `10/16 = 62.5%`).
   - On `T3` (`Strcam` + `Aptrow`), simulated neutral incongruence median is **`41.7–47.6%`** (`95%` range `15.4–72.7%`, two-sided $p = 0.072–0.092$ vs. observed `T3` incongruence `11/16 = 68.8%`).

---

## 4. Scientific Conclusion (`CORRIGIDO`)

Applying the actual empirical per-species codon coverage mask (`empirical_dual_aligner_mask`) completely resolves the discrepancy between the `v1.0.1` unmasked 878-codon simulation (`25.0–29.0` shared PTCs, two-sided $p = 0.004–0.052$) and the empirical dataset (`16` shared PTCs): under the coverage-masked neutral model, the observed **`16` shared PTCs** sit right at the center of the neutral distribution (**median `16.0` at $t_{\text{stem}}=0.00$, `48.8th` percentile, $p = 1.000$; `18.0` at $t_{\text{stem}}=0.02$, `33.8th` percentile, $p = 0.778$**), and the observed **`37.5%` incongruence** matches both the unconditional median (**`40.0–44.4%`**, $p = 0.618–0.848$) and the exact conditional median at $N_{\text{shared}} = 16$ (**`37.5%`**). Both the unmasked and coverage-masked results are transparently reported in Section 2.5 (#6), Section 3.4, Section 4.4, and `Table_S8`.
