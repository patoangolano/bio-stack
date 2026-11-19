# Phase A Audit Report — Crown-Neornithes Synapomorphy Detectability Limits, Exon Coverage Breakdown, and Ancestral-Stop Survival Power

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final` (baseline `enamel-relics-v1.0.2`, tag `db71156`, Zenodo Version DOI `10.5281/zenodo.24368656`)  
**Primary Script & Output:** [`scripts/enamel_relics/13_audit_final_neornithes_and_stats.py`](../../scripts/enamel_relics/13_audit_final_neornithes_and_stats.py) · [`runs/enamel_relics/v102/neornithes_detectability_audit.json`](../../runs/enamel_relics/v102/neornithes_detectability_audit.json)

---

## 1. Problem Statement (`Phase A`)

In `enamel-relics-v1.0.2`, the Abstract, Section 2.5 (#6), Figure 3B, Table 3, and Section 3.4 used shorthand formulations such as `"zero crown Neornithes synapomorphies"` and `"0 congruent lesions at crown Neornithes"`. However, in the empirical 7-bird `ENAM_exon8_aves` alignment (`878` *Alligator* reference codons), per-species codon coverage ranges from `22.7%` (`Taegut`, `199/878` codons) to `82.5%` (`Aptrow`, `724/878` codons), and **zero codons (`0/878 = 0.0%`) are simultaneously covered across all 7 birds**. Stating `"zero crown Neornithes synapomorphies"` without explicit qualification risks conflating **non-detection under incomplete coverage and secondary substitution** with **proven biological absence of crown-Neornithes inactivating mutations**.

---

## 2. A1 — Forensic Audit of How the Pipeline Evaluates Crown-Neornithes Mutations

Inspection of [`scripts/enamel_relics/04_run_level4_alignments.py`](../../scripts/enamel_relics/04_run_level4_alignments.py), [`scripts/enamel_relics/10_audit_mutations_and_asr.py`](../../scripts/enamel_relics/10_audit_mutations_and_asr.py), and [`scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py`](../../scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py) confirms that our pipeline **did not require 7/7 simultaneous bird coverage** to call or evaluate a shared mutation. Instead, it operated at two distinct levels:

1. **Empirical Variant Calling & Missing-Data-Tolerant Parsimony / ML ASR (`Table_S5a`):**
   - Any stop codon or frameshift indel shared by $\ge 2$ birds at the same reference coordinate in both `MAFFT L-INS-i` and `MUSCLE 5` was retained (`16` shared PTCs + `1` frameshift indel in `ENAM_exon8_aves`).
   - Uncovered species (`---`) were treated as missing data (`None`) in exact Fitch parsimony (`min_stop_gains_no_reversal`) and in `IQ-TREE 3.1.4` (`GTR+G --ancestral`) ML ancestral state reconstruction.
   - In the K80 simulation (`12_audit_v102_variants_phylo_sim.py`, line 869), a site was counted as an uncontradicted crown-Neornithes synapomorphy (`lut_crown_neorn`) whenever it had **$\ge 1$ Palaeognathae stop carrier AND $\ge 1$ Neognathae stop carrier AND `0` covered birds carrying a sense codon** (`(n_sense == 0) and pal_stop and neog_stop`).
2. **Why Zero Uncontradicted Crown-Neornithes Sites Were Detected Empirically:**
   - Across the 16 shared avian *ENAM* Exon 8 PTCs in `Table_S5a`, **four sites actually have stop carriers in both Palaeognathae and Neognathae**:
     - `Codon 439` (`XP_019346377.1:p.Ala439Ter`): `Aptrow + Dronov` (`TGA`) and `Anapla` (`TAG`), **but** covered `Strcam` (`AGA`), `Galgal` (`TCA`), and `Strhab` (`CTG`) carry sense codons (crown-Neornithes ML $P_{\text{stop}} = 0.2207$).
     - `Codon 562` (`p.Trp562Ter`): `Strcam` (`TGA`) and `Strhab` (`TGA`), **but** covered `Aptrow, Dronov, Galgal, Anapla` all retain ancestral `TGG` (crown-Neornithes ML $P_{\text{stop}} = 0.0314$).
     - `Codon 615` (`p.Trp615Ter`): `Dronov + Strcam` (`TGA`) and `Strhab` (`TAA`), **but** covered `Aptrow, Galgal, Anapla` retain ancestral `TGG` (crown-Neornithes ML $P_{\text{stop}} = 0.3035$).
     - `Codon 684` (`p.Gln684Ter`): `Aptrow + Dronov + Strcam` (`TAG`) and `Galgal` (`TAA`), **but** covered `Anapla` and `Strhab` retain ancestral `CAA` (crown-Neornithes ML $P_{\text{stop}} = 0.3655$).
   - Thus, zero sites in the empirical alignment satisfy `(n_sense == 0) and pal_stop and neog_stop`, and no site reaches crown-Neornithes ML $P_{\text{stop}} \ge 0.50$.

### 2.1. Cross-Clade Codon Coverage Breakdown Across All Level 4 Avian Exons (`13_audit_final_neornithes_and_stats.py`)

| Dataset | Reference Codons | Birds in Alignment | All 7 Birds (`7/7`) | All 6 Non-Passerine Birds (`6/6`) | $\ge 6$ of 7 Birds | $\ge 2$ Palaeo + $\ge 2$ Neog | $\ge 1$ Palaeo + $\ge 1$ Neog (`Table_S8` Def.) | Single-Clade Only or `0` Birds |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `ENAM_exon8_aves` (dual-aligner) | `878` | `7` (`Strcam, Aptrow, Dronov, Galgal, Anapla, Strhab, Taegut`) | **`0` (`0.0%`)** | `306` (`34.9%`) | `449` (`51.1%`) | `675` (`76.9%`) | `727` (`82.8%`) | `151` (`17.2%`) |
| `ENAM_exon8_aves` (MAFFT only) | `878` | `7` (`Strcam, Aptrow, Dronov, Galgal, Anapla, Strhab, Taegut`) | **`0` (`0.0%`)** | `348` (`39.6%`) | `514` (`58.5%`) | `747` (`85.1%`) | `799` (`91.0%`) | `79` (`9.0%`) |
| `AMBN_exon3_aves` (dual-aligner) | `39` | `5` (`Strcam, Aptrow, Dronov, Galgal, Anapla`) | **`0` (`0.0%`)** | `0` (`0.0%`) | `0` (`0.0%`) | `34` (`87.2%`) | `36` (`92.3%`) | `3` (`7.7%`) |
| `AMBN_exon4_aves` (dual-aligner) | `80` | `3` (`Strcam, Aptrow, Dronov`; `Anapla` only in `-word_size 7`) | **`0` (`0.0%`)** | `0` (`0.0%`) | `0` (`0.0%`) | **`0` (`0.0%`)** | **`0` (`0.0%`)** | `80` (`100.0%`) |
| `AMBN_exon9_aves` (dual-aligner) | `175` | `3` (`Strcam, Aptrow, Dronov`) | **`0` (`0.0%`)** | `0` (`0.0%`) | `0` (`0.0%`) | **`0` (`0.0%`)** | **`0` (`0.0%`)** | `175` (`100.0%`) |

### 2.2. Four Epistemological Categories Explicitly Distinguished in the Revised Manuscript

1. **No mutation detected (`Nenhuma mutação detectada`):** Zero uncontradicted crown-Neornithes coding lesions (`(n_sense == 0) and pal_stop and neog_stop`) were observed in the sampled alignments.
2. **Proven absence of shared mutations (`Ausência comprovada de mutações compartilhadas`):** **NOT supported** by our data, because (a) most of the *ENAM*/*AMBN*/*AMTN* coding sequence outside *ENAM* Exon 8 is physically deleted or unrecovered in Neognathae, (b) `0/878` codons of *ENAM* Exon 8 have 7-bird coverage (`17.2%` lack dual-aligner cross-clade coverage altogether), and (c) secondary substitutions over >100 Myr frequently revert ancestral stop codons to sense codons in one or more descendants.
3. **Absence of evidence due to insufficient coverage (`Ausência de evidência por cobertura insuficiente`):** Applies strictly to universal 7-bird coverage (`0/878` codons), to the `151/878` (`17.2%`) *ENAM* Exon 8 codons lacking dual-aligner cross-clade coverage, to baseline *AMBN* Exons 4 and 9 (`0%` Neognathae coverage), and to unrecovered exons (*ENAM* Exons 1–5, 7; *AMBN* Exons 1–2, 5–8).
4. **Probabilistic ancestral inference despite missing data (`Inferência ancestral probabilística apesar de dados faltantes`):** Evaluated via `IQ-TREE 3` (`GTR+G --ancestral`) across all 16 shared PTCs; while subclade nodes reach $P_{\text{stop}} = 0.880–0.999$ (and crown Neognathae reaches $P_{\text{stop}} = 0.926$ at codon 1022 where `Anapla` is missing and `Taegut = AAA`), crown Neornithes never exceeds $P_{\text{stop}} = 0.366$ (at codon 684).

---

## 3. A2 — Quantitative Detection Power for a Genuine Crown-Neornithes Ancestral Stop (`13_audit_final_neornithes_and_stats.py`)

To quantify the probability of detecting a genuine crown-Neornithes ancestral PTC if one had been fixed at the crown-Neornithes root node prior to the Palaeognathae–Neognathae divergence, we simulated `87,800` ancestral stop trials (`878` codons $\times$ `100` replicates, equally distributed across `TAA`, `TAG`, and `TGA`, `seed 20261124`) evolving neutrally under K80 ($\kappa = 4.0$) along empirical `T1` branch lengths at branch-length scales $s \in \{0.00, 0.65, 1.00, 1.25\}$ and observed through each coverage mask:

| Coverage Mask | Branch Scale ($s$) | Physical Regime | Strict `7/7` Universal Stop | Uncontradicted `6/6` Non-Passerine Stop | Uncontradicted $\ge 2$ Palaeo + $\ge 2$ Neog Stop | Uncontradicted $\ge 1$ Palaeo + $\ge 1$ Neog Stop (`Table_S8` Def.) | Corrupted by $\ge 1$ Sense Reversal (Appears Parsimony-Incongruent) | Truncated to 1 Clade by Missing Coverage |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| `unmasked_full_878` | `0.00` | Pure coverage effect (`0` substitutions) | `100.0%` | `100.0%` | `100.0%` | `100.0%` | `0.0%` | `0.0%` |
| `unmasked_full_878` | `0.65` | Reduced neutral drift + full coverage | `22.7%` | `22.7%` | `22.7%` | `22.7%` | `59.4%` | `0.0%` |
| `unmasked_full_878` | `1.00` | Empirical neutral drift + full coverage | `10.7%` | `10.7%` | `10.7%` | `10.7%` | `63.3%` | `0.0%` |
| `unmasked_full_878` | `1.25` | Accelerated drift + full coverage | `6.4%` | `6.4%` | `6.4%` | `6.4%` | `62.3%` | `0.0%` |
| `empirical_mafft_mask` | `0.00` | Pure MAFFT coverage mask (`0` subst.) | **`0.0%`** | `39.6%` | `85.1%` | `91.0%` | `0.0%` | `4.0%` |
| `empirical_mafft_mask` | `0.65` | Reduced drift + MAFFT mask | **`0.0%`** | `12.8%` | `28.1%` | `30.7%` | `41.3%` | `2.3%` |
| `empirical_mafft_mask` | `1.00` | Empirical drift + MAFFT mask | **`0.0%`** | `7.1%` | `16.0%` | `17.7%` | `45.2%` | `1.8%` |
| `empirical_dual_aligner_mask` | `0.00` | Pure dual-aligner mask (`0` subst.) | **`0.0%`** | `34.9%` | `76.9%` | `82.8%` | `0.0%` | `3.4%` |
| `empirical_dual_aligner_mask` | `0.65` | Reduced drift + dual-aligner mask | **`0.0%`** | `11.2%` | `25.4%` | `28.1%` | `37.4%` | `2.0%` |
| **`empirical_dual_aligner_mask`** | **`1.00`** | **Empirical drift + dual-aligner mask** | **`0.0%`** | **`6.2%`** | **`14.5%`** | **`16.3%`** | **`40.8%`** | **`1.5%`** |
| `empirical_dual_aligner_mask` | `1.25` | Accelerated drift + dual-aligner mask | **`0.0%`** | `4.2%` | `9.9%` | `11.3%` | `40.7%` | `1.3%` |

### Key Conclusions from the Power Audit:
1. **Strict 7-bird universal detection power is identically `0.0%`** under both empirical masks because `0/878` codons are simultaneously covered across all 7 birds.
2. **Even under the most permissive uncontradicted cross-clade criterion ($\ge 1$ Palaeognathae stop + $\ge 1$ Neognathae stop with `0` sense reversals anywhere), detection power for a single ancestral crown-Neornithes PTC under empirical neutral branch lengths (`scale = 1.0`) is only `16.3%`** (`14.5%` requiring $\ge 2$ carriers per clade; `6.2%` across all 6 non-passerine birds).
3. **In `40.8%` of trials at `scale = 1.0`, secondary stop-to-sense substitutions in one or more covered descendants convert a genuine ancestral crown-Neornithes stop into an apparently parsimony-incongruent cross-clade polymorphism** (indistinguishable from parallel gains at sites such as codons `439, 562, 615, 684`).
4. Therefore, **the absence of an uncontradicted crown-Neornithes synapomorphy in our 7-bird alignment cannot be used to infer the absence of crown-avian tooth loss or crown-avian pseudogenization**, and the status of crown-Neornithes-wide synapomorphies in this 8-bird RefSeq panel is explicitly declared **unresolved by the sampled data**.

---

## 4. A3 — Summary of Manuscript & Figure Corrections Applied

- **Abstract (`manuscript.typ` L119 / `manuscript.md` L14):** Replaced `"zero crown Neornithes synapomorphies, and 31.2–37.5% avian ENAM stop incongruence... matching coverage-masked neutral simulations"` with `"and 31.2–37.5% avian ENAM stop incongruence under topology T1 (56.2–68.8% under alternatives), consistent with coverage-masked K80 simulations; no crown-Neornithes-wide synapomorphy was detected, as incomplete coverage (0/878 seven-bird codons) limits inference."`
- **Results 2.5 (#6, `manuscript.typ` L278):** Retitled and rewritten to report the exact cross-clade coverage breakdown (`0/878 = 0.0%` 7-bird; `306/878 = 34.9%` 6 non-passerine; `727/878 = 82.8%` $\ge 1$ per clade; `675/878 = 76.9%` $\ge 2$ per clade) and the `16.3%` neutral uncontradicted ancestral-stop detection power (`40.8%` secondary reversal corruption).
- **Figure 3B & Caption (`06_make_figures.py` L778 / `manuscript.typ` L288):** Updated Figure 3B badge to `"Crown Neornithes: none detected (0/878 7-bird cov.; T1: 5–6/16 incong.; T2/T3: 9–11/16)"` and caption to `"no uncontradicted crown-Neornithes lesion detected [0/878 codons simultaneously covered across all 7 birds]"`.
- **Table 3 (`manuscript.typ` L303):** Updated row to `Crown Neornithes (uncontradicted) | Palaeognathae + Neognathae | 0 | 0 | None detected (unresolved: 0/878 codons covered in all 7 birds; 16.3% neutral ancestral-stop detection power)`.
- **Discussion 3.4 & Limitations 3.5 (#2):** Explicitly state that non-detection of an uncontradicted crown-Neornithes synapomorphy does not prove absence of crown-avian inactivation and remains compatible with the single crown-avian tooth loss (~116 Ma) inferred from 48-bird frameshift indels by Meredith et al. (2014).
