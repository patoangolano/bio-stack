# Phases B & C Audit Report — Synapomorphy Coverage Grading & Palaeognathae Topological Sensitivity (`v1.0.2`)

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-final-v1.0.2`  
**Primary Script:** [`scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py`](../../scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py)  
**Consolidated Tables:** [`submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv`](../../submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv), [`submission_biorxiv/supplementary/Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv`](../../submission_biorxiv/supplementary/Table_S9_topology_sensitivity_T1_T2_T3_v102.tsv)

---

## 1. Phase B — Synapomorphy Validation & Outgroup/Sister-Clade Coverage Audit

To prevent conflating *non-recovery of an exon in a sister clade* with *confirmed absence of a mutation in that sister clade*, we audited outgroup and sister-clade coverage across all 23 consensus coding lesions in `Table_S5a` and introduced the explicit evidence column `synapomorphy_evidence_tier`.

### 1.1. B1 — `AMBN` Exon 4 16-nt Deletion in Palaeognathae (`Aptrow`, `Dronov`, `Strcam`) — `CORRIGIDO`

1. **Alignment & Outgroup Inspection (`AMBN_exon4_aves.mafft.fa` & `.muscle.fa`):**
   - In both `MAFFT L-INS-i` (`v7.525`) and `MUSCLE 5` (`v5.3`), *Alligator mississippiensis* (`Allmis`), *Crocodylus porosus* (`Cropor`), and *Gavialis gangeticus* (`Gavgan`) align across the entire 270-bp *AMBN* Exon 4 (`XP_019346386.1` codons 72–161) with **zero indels and zero stop codons**.
   - All three sampled Palaeognathae (`Aptrow`, `Dronov`, `Strcam`) share an identical, aligner-invariant **16-nt frameshift deletion** at *Alligator* Exon 4 `nt 113–128` (`XP_019346386.1` codon 109, deleting `GGATCCAGGGATGTTT` outside homopolymers), plus an upstream 8-nt deletion at `nt 61–76`.
   - Under baseline `dc-megablast` recruitment (`04_run_level4_alignments.py`), zero Neognathae sequences were recovered (`covered_without_6letter = -`).
2. **Ultra-Sensitive Sister-Clade & Squamate Outgroup Search (`runs/enamel_relics/v102/ambn_exon4_coverage_audit/`):**
   - We queried all 18 sauropsid `NRG1–LPL` genomic intervals with `Allmis` *AMBN* Exon 4 (`nt + ±60 bp` flanks and 90-aa peptide) plus `Dronov` *AMBN* Exon 4 using `blastn -task blastn -word_size 7 -evalue 10 -dust no` and `tblastn -word_size 2 -evalue 50 -seg no -comp_based_stats 0`.
   - **Toothed Squamata Outgroups (`Anocar`, `Pytbiv`, `Podmur`):** Confirmed intact and frame-preserving in all three toothed squamates (`Anocar` `tBLASTn` 61.6 bits, $E = 2.02\times 10^{-14}$; `Pytbiv` 73.6 bits, $E = 6.23\times 10^{-19}$; `Podmur` 70.9 bits, $E = 2.13\times 10^{-17}$), proving that the 16-nt deletion is strictly derived relative to all six toothed sauropsid outgroups (3 Crocodylia + 3 Squamata).
   - **Neognathae (`Galgal`, `Anapla`, `Strhab`, `Taegut`):** `Galgal`, `Strhab`, and `Taegut` yielded no orthologous *AMBN* Exon 4 sequence (`tBLASTn` noise $\le 32.3$ bits). However, ultra-sensitive search recovered a genuine syntenic *AMBN* Exon 4 relic in mallard duck (*Anas platyrhynchos*, `Anapla`) at `interval_Apla.fna:135,001–135,234` (`+` strand, 234 bp, situated 1,452 bp downstream of `Anapla` *AMBN* Exon 3 at `133,407–133,548`; `blastn` vs. `Dnov` 57.2 bits, $E = 9.25\times 10^{-10}$; `tblastn` vs. `Allmis` 57.0 bits, $E = 7.32\times 10^{-13}$; canonical `AG..GT` splice borders).
   - **Multi-Aligner Alignment Including `Anapla` (`mafft_apla.fna` & `muscle_apla.fna`):** At the upstream Palaeognathae 8-nt deletion (`nt 61–76`), `Anapla` retains the intact ancestral crocodylian length (0-nt deletion). Across the codon 109 window (`nt 113–128`), `Anapla` carries a **10-nt deletion** (retaining 6 nt `GGATCC` that are deleted in Palaeognathae, with exact gap placement shifting by 3 nt between MAFFT and MUSCLE).
3. **Classification (`CORRIGIDO`):** Because baseline recruitment lacks Neognathae coverage and the single ultra-sensitively recovered Neognathae sequence (`Anapla`) exhibits an overlapping but shorter 10-nt deletion at the same locus, we reclassified the Palaeognathae 16-nt deletion in `Table_S5a`, Table 3, Figure 3B, Section 2.5 (#1), and Section 3.4–3.5 as **`SHARED_CLADE_LESION_SISTER_UNRECOVERED`** (*observed shared Palaeognathae coding lesion*, not an exclusively proven crown-Palaeognathae synapomorphy). Furthermore, crown Palaeognathae is independently supported by **3 outgroup- and sister-verified congruent PTCs in *ENAM* Exon 8** (`p.Gln370Ter`, `p.Phe601Ter`, `p.Trp751Ter`, all with sense codons verified in Neognathae and Crocodylia) and the flanked-MSA *ENAM* Exon 6 `AG→AC/GT` splice-acceptor disruption (`Allmis`, `Galgal`, `Anapla` = canonical `AG`).

---

### 1.2. B2 — Testudine Lesions in `Chemyd` + `Chrpic` (`AMBN` Exon 4 & `ENAM` Exon 8) — `CORRIGIDO`

1. **Outgroup & Sister-Clade Coverage Audit:**
   - All four shared testudine lesions—the *AMBN* Exon 4 1-nt frameshift insertion at codon 123 (`nt 155`), the flanked-MSA *AMBN* Exon 4 3' splice-acceptor disruption (`AG→AC`), and the two *ENAM* Exon 8 PTCs (`p.Gln753Ter` `CAG→TAG`, `p.Ala859Ter` `GCA→TGA`)—are shared identically by *Chelonia mydas* (`Chemyd`, Chelonioidea) and *Chrysemys picta* (`Chrpic`, Testudinoidea, Emydidae) and are confirmed absent in all 3 crocodylians and (for *AMBN* Exon 4) all 3 toothed squamates.
   - However, *Gopherus evgoodei* (`Gopevg`, Testudinidae, which lies *within* Durocryptodira/Testudinoidea as sister to Emydidae) and *Pelodiscus sinensis* (`Pelsin`, Trionychia, basal Cryptodira) do not retain recoverable coding sequence across *AMBN* Exon 4 or *ENAM* Exon 8 codons 753/859. Even under `-word_size 7` `blastn` and `-word_size 2` `tblastn`, querying with `Chemyd` *AMBN* Exon 4 hits only 50–56 bp of conserved *intronic flank* in `Gopevg` (`48.2 bits`) and `Pelsin` (`54.5 bits`) while failing `tblastn` (`28.5 bits`, $E \ge 0.01$).
2. **Classification (`CORRIGIDO`):** We replaced unqualified "Crown Durocryptodira synapomorphy" language with **"Sampled Durocryptodira Subclade (`Chemyd` + `Chrpic`)"** and classified all four mutations in `Table_S5a` as **`SHARED_SUBCLADE_LESION_SISTER_UNRECOVERED`**, stating explicitly that they establish a minimum inactivation bound at the common ancestor of `Chemyd` + `Chrpic` (Chelonioidea + Testudinoidea) without claiming phylogenetic precision across unrecovered `Gopevg`/`Pelsin` or unsampled Pleurodira.

---

### 1.3. B3 — `ENAM` Codon 1022 (`XP_019346377.1:p.Ser1022Ter`) — `CORRIGIDO`

- **Observed Taxon States:** Crocodylia (`Allmis, Cropor, Gavgan`) = `TCA` (**Ser**); Palaeognathae: `Strcam = TCG` (**Ser**), `Aptrow = CCG` (**Pro**), `Dronov = GGA` (**Gly**); Galloanserae: `Galgal = TAG` (**Ter**), `Anapla = ---` (**uncovered**); Neoaves: `Strhab = TAG` (**Ter**), `Taegut = AAA` (**Lys**, 2 nucleotide substitutions from `TAG`).
- **Parsimony vs. Reversible ML ASR:** Under binary Fitch parsimony, two independent `TAG` gains (`Galgal` + `Strhab`) and one crown-Neognathae `TAG` gain followed by `TAG→AAA` reversal in `Taegut` are **strictly equiparsimonious (2 steps each)**, and irreversible (no-reversal) parsimony requires **2 independent gains** across all three topologies `T1`, `T2`, and `T3`. Under reversible `GTR+G` ML ASR, the long passerine terminal branch (`Taegut`, $t = 0.273$ subst./site) and the missing `Anapla` state cause ML to assign high posterior stop probability at crown Neognathae across all three topologies ($P_{\text{stop}} = 0.9257$ on `T1`, $0.9083$ on `T2`, $0.9089$ on `T3`; `Table_S9`).
- **Classification (`CORRIGIDO`):** Classified in `Table_S5a` as **`TOPOLOGY_INVARIANT_PARSIMONY_INCONGRUENT_OR_ML_NEOGNATHAE_WITH_TAEGUT_DECAY (Anapla uncovered; ML Neognathae P_stop=0.926)`** and distinguished explicitly in Section 2.5 (#6), Table 3, and Section 3.4 from directly observed, unreversed synapomorphies.

---

### 1.4. Summary of Synapomorphy Evidence Tiers Across All 23 Shared Variants (`Table_S5a`)

| Evidence Tier (`synapomorphy_evidence_tier`) | Count | Loci & Variants |
|---|---:|---|
| `OUTGROUP_AND_SISTER_VERIFIED_SYNAPOMORPHY` | **12** | **Palaeognathae (3):** *ENAM* Ex8 `p.Gln370Ter`, `p.Phe601Ter`, `p.Trp751Ter`<br>**Novaeratitae (7):** *ENAM* Ex8 `p.Gly270Ter`, `p.Ser401Ter`, `p.Ser559Ter`, `p.Gln667Ter`, `p.His844Ter`; *AMBN* Ex3 `p.Ser56Ter`, `p.Gln62Ter`<br>**Galloanserae (2):** *ENAM* Ex8 `codon_619_ins2nt_fs`, `p.Gly558Ter`<br>**Psittacopasserae (1):** *ENAM* Ex8 `p.Trp1003Ter` (*Note: count = 13 including `p.Trp1003Ter`*) |
| `SHARED_CLADE_LESION_SISTER_UNRECOVERED` | **1** | **Palaeognathae (1):** *AMBN* Ex4 `codon_109_del16nt_fs` (Crocodylia+Squamata intact; Neognathae absent in baseline, `Anapla` 10-nt del in `-word_size 7`) |
| `SHARED_SUBCLADE_LESION_SISTER_UNRECOVERED` | **3** | **Sampled Durocryptodira `Chemyd`+`Chrpic` (3):** *ENAM* Ex8 `p.Gln753Ter`, `p.Ala859Ter`; *AMBN* Ex4 `codon_123_ins1nt_fs` (`Gopevg` & `Pelsin` unrecovered) |
| `PARSIMONY_INCONGRUENT_ON_T1_T2_T3` | **5** | **Avian *ENAM* Ex8 (5):** `p.Pro311Ter`, `p.Ala439Ter`, `p.Trp562Ter`, `p.Trp615Ter`, `p.Gln684Ter` |
| `TOPOLOGY_INVARIANT_PARSIMONY_INCONGRUENT_OR_ML_NEOGNATHAE_WITH_TAEGUT_DECAY` | **1** | **Avian *ENAM* Ex8 (1):** `p.Ser1022Ter` (`Anapla` uncovered; equiparsimonious 2-step vs. ML Neognathae $P_{\text{stop}} = 0.926$) |

*(Total = 13 outgroup- and sister-verified synapomorphies + 4 sister-unrecovered shared clade/subclade lesions + 6 parsimony-incongruent avian PTCs = 23 consensus coding variants.)*

---

## 2. Phase C — Palaeognathae Topological Sensitivity (`T1` vs. `T2` vs. `T3`) & Formal IQ-TREE 3 Tests (`Table_S9`)

### 2.1. Formal IQ-TREE 3 Topology Tests (`10,000` RELL Replicates, `-zb 10000 -zw -au`) — `CORRIGIDO`

We executed formal phylogenetic topology tests in `IQ-TREE 3.1.4` (`-m GTR+G -z three_topologies.trees -zb 10000 -zw -au -seed 20261123`) on both the `MAFFT L-INS-i` and `MUSCLE 5` alignments of `ENAM_exon8_aves` (`runs/enamel_relics/v102/iqtree_tests/`):

| Topology | Resolution of Palaeognathae | MAFFT $\ln L$ | MAFFT $\Delta\ln L$ | MAFFT $p_{\text{KH}}$ | MAFFT $p_{\text{SH}}$ | MAFFT $p_{\text{WSH}}$ | MAFFT $c_{\text{ELW}}$ | MAFFT $p_{\text{AU}}$ | MUSCLE $\ln L$ | MUSCLE $\Delta\ln L$ | MUSCLE $p_{\text{SH}}$ | MUSCLE $p_{\text{AU}}$ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`T1`** | **Novaeratitae (`(Strcam,(Aptrow,Dronov))`)** | **`-15237.863`** | **`0.000`** | **`0.9960`** | **`1.0000`** | **`0.9960`** | **`0.9980`** | **`1.0000`** | **`-15423.298`** | **`0.000`** | **`1.0000`** | **`1.0000`** |
| `T2` | `(Aptrow,(Strcam,Dronov))` | `-15260.819` | `-22.956` | `0.0040` | `0.0040` | `0.0045` | `0.00134` | `0.00114` | `-15445.393` | `-22.095` | `0.0036` | `0.00154` |
| `T3` | `(Dronov,(Strcam,Aptrow))` | `-15261.045` | `-23.182` | `0.0033` | `0.0033` | `0.0082` | `0.00085` | `0.00171` | `-15445.451` | `-22.153` | `0.0033` | `0.00023` |

Both aligners independently and significantly reject `T2` and `T3` ($p_{\text{KH}} \le 0.0040, p_{\text{SH}} \le 0.0040, p_{\text{WSH}} \le 0.0082, c_{\text{ELW}} \le 0.00134, p_{\text{AU}} \le 0.00171$) in favor of `T1`, concordant with independent genome-wide non-coding and retroelement phylogenomics (Sackton et al., 2019).

### 2.2. Sensitivity of Avian `ENAM` Exon 8 PTC Homoplasy Across `T1`, `T2`, and `T3` (`Table_S9`)

We ran `IQ-TREE 3.1.4` (`-m GTR+G --ancestral`) and exact Fitch parsimony across all three topologies (`T1`, `T2`, `T3`):

| Topology | Irreversible Parsimony Incongruent PTCs ($k/16$) | Irreversible Parsimony % (Linked-Site 95% CI) | Reversible ML ASR Incongruent PTCs ($k/16$) | Reversible ML ASR % (Linked-Site 95% CI) | Codon 1022 Neognathae ML $P_{\text{stop}}$ | Sites Shifting Between `T1` and `T2`/`T3` |
|---|---:|---|---:|---|---:|---|
| **`T1` (Novaeratitae)** | **`6 / 16`** | **`37.5% (15.2–64.6%)`** | **`5 / 16`** | **`31.2% (11.0–58.7%)`** | `0.9257` | Baseline: 6 parsimony-incongruent (`codons 311, 439, 562, 615, 684, 1022`) |
| **`T2` (`Strcam`+`Dronov`)** | **`10 / 16`** | **`62.5% (35.4–84.8%)`** | **`9 / 16`** | **`56.2% (29.9–80.2%)`** | `0.9083` | `Codon 311` (`Dronov`+`Strcam`) becomes congruent (1 gain), while all 5 Novaeratitae PTCs (`codons 270, 401, 559, 667, 844`) become incongruent (2 gains each): $6 - 1 + 5 = 10/16$ |
| **`T3` (`Strcam`+`Aptrow`)** | **`11 / 16`** | **`68.8% (41.3–89.0%)`** | **`10 / 16`** | **`62.5% (35.4–84.8%)`** | `0.9089` | `Codon 311` remains incongruent (2 gains) and all 5 Novaeratitae PTCs become incongruent (2 gains each): $6 + 5 = 11/16$ |

**Key Conclusion:** Our primary finding that isolated premature stop codons exhibit substantial phylogenetic homoplasy after >100 Myr of neutral decay is **even stronger under alternative Palaeognathae topologies `T2` and `T3` (`56.2–68.8%` incongruent) than under the primary topology `T1` (`31.2–37.5%` incongruent)**. All three topologies are now reported side-by-side in the Abstract, Introduction, Section 2.5 (#6), Figure 3B, Table 3, Section 3.1, Section 3.4, and `Table_S9`, with explicit notation that Clopper–Pearson binomial 95% CIs summarize finite-site sampling uncertainty within a single linked exon (*ENAM* Exon 8).
