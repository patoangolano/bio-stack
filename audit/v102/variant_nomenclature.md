# Phase A Audit Report — Molecular Variant Nomenclature & Reference vs. Avian Ancestral Codons (`v1.0.2`)

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-final-v1.0.2` (baseline `enamel-relics-v1.0.1`, commit `13b34a6` / `d18be9d`)  
**Primary Script:** [`scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py`](../../scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py)  
**Consolidated Table:** [`submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv`](../../submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv)

---

## 1. Executive Summary

In `enamel-relics-v1.0.1`, all 23 Level 4 dual-aligner consensus coding mutations in `Table_S5a` were indexed using *Alligator mississippiensis* RefSeq reference protein coordinates (`XP_019346377.1:p.Xxx###Ter` for *ENAM*, 1,078 aa; `XP_019346386.1:p.Xxx###Ter` for *AMBN*, 390 aa). However, shorthand phrasing in Section 2.5 and Table 3 of the `v1.0.1` manuscript juxtaposed the *Alligator* reference amino-acid designation directly against the reconstructed avian codon transition without explicitly stating that a non-synonymous sense substitution occurred along the ~240–250 Myr crocodylian–avian stem prior to stop-codon formation:

- **`ENAM` codon 601 (`XP_019346377.1:p.Phe601Ter`):** Written as ``p.Phe601Ter` avian `TAT→TAA``, without noting that *Alligator* carries `TTT` (**Phe**) whereas crown Neornithes reconstructed `TAT` (**Tyr**, $P = 0.8992$ in MAFFT / $0.8702$ average across MAFFT+MUSCLE) prior to the Palaeognathae `TAT (Tyr) → TAA (Ter)` transversion.
- **`ENAM` codon 1003 (`XP_019346377.1:p.Trp1003Ter`):** Written as ``p.Trp1003Ter`, avian `TGC→TGA``, without noting that *Alligator* carries `TGG` (**Trp**) and that `TGC` (**Cys**, $P = 0.3433$) is an artifact of multiplying independent nucleotide-site ML base marginals (`T × G × C`) at crown Neornithes (`TGC` $P = 0.3433$ vs. `TGA` $P = 0.2907$ and `TGG` $P = 0.1678$), even though **no extant bird carries `TGC`**.

We audited all 23 shared coding variants in `Table_S5a` against the exact *Alligator* coding sequences and protein references, verified every codon-to-amino-acid translation programmatically, added six explicit columns to `Table_S5a` (`ref_aa_3letter`, `croc_consensus_aa_3letter`, `neornithes_ml_aa_3letter`, `neornithes_ml_p_best`, `neognathae_ml_p_stop`, `inferred_substitution_path`), updated Figure 3B–C and Table 3, revised Sections 2.5 and 3.4 of the manuscript, and added automated translation assertions to [`tests/test_enamel_relics_audit.py`](../../tests/test_enamel_relics_audit.py).

---

## 2. Detailed Forensic Analysis of `ENAM` Codon 601 and Codon 1003

### 2.1. A1 — `ENAM` Codon 601 (`XP_019346377.1:p.Phe601Ter`, Exon 8 `nt 1201–1203`) — `CORRIGIDO`

| Parameter | Empirical Value / Reconstructed State |
|---|---|
| Reference protein & 1-based residue | `XP_019346377.1` residue `601` = **`F` (`Phe`, Phenylalanine)** |
| *Alligator mississippiensis* Exon 8 coordinate | `nt 1201–1203` (1-based within 2,634-bp Exon 8; complete 8-exon CDS = 3,237 bp) |
| Crocodylia observed codons (`Allmis`, `Cropor`, `Gavgan`) | `TTT`, `TTT`, `TTT` (**Phe** in all 3 crocodylian outgroups; 0 lesions) |
| Palaeognathae observed codons (`Aptrow`, `Dronov`, `Strcam`) | `TAA`, `TAA`, `TAA` (**Ter** in all 3 sampled palaeognaths) |
| Neognathae observed codons (`Galgal`, `Anapla`, `Strhab`, `Taegut`) | `Galgal = TAT` (**Tyr**), `Anapla = TAC` (**Tyr**), `Strhab = TAT` (**Tyr**), `Taegut = ---` (uncovered) |
| Crown Neornithes (`Node1`, `IQ-TREE 3.1.4 GTR+G`) ML reconstruction | **`TAT` (`Tyr`, Tyrosine)**: $P(\text{TAT}) = 0.8992$ (MAFFT), $0.8412$ (MUSCLE), dual-aligner mean $P = 0.8702$; $P_{\text{stop}} = 0.0480$ |
| Crown Neognathae (`Node5`) ML reconstruction | **`TAT` (`Tyr`)**: $P(\text{TAT}) = 0.9461$ (MAFFT), $P_{\text{stop}} = 0.0044$ |
| Crown Palaeognathae (`Node2`) ML reconstruction | **`TAA` (`Ter`)**: $P_{\text{stop}} = 0.9944$ (`T1` MAFFT), $0.9998$ (Novaeratitae `Node3`) |
| Evolutionary pathway | $\text{TTT (Phe, Crocodylia)} \xrightarrow{\text{stem T}\to\text{A}} \text{TAT (Tyr, crown Neornithes)} \xrightarrow{\text{Palaeognathae T}\to\text{A}} \text{TAA (Ter, Palaeognathae)}$ |

**Resolution:** Because all variant identifiers in `Table_S5a` are explicitly anchored to the *Alligator* reference protein (`XP_019346377.1`), `XP_019346377.1:p.Phe601Ter` is the mathematically exact reference-anchored identifier (changing it to `p.Tyr601Ter` under `XP_019346377.1:` would falsify residue 601 of `XP_019346377.1`, which is `Phe`). Instead, `Table_S5a`, Table 3, and Section 2.5 now explicitly state both the *Alligator* reference codon/residue (`TTT` [Phe]) and the reconstructed crown-Neornithes ancestral codon/residue (`TAT` [Tyr], $P = 0.870–0.899$) that underwent the 1-step transversion `TAT (Tyr) → TAA (Ter)` in crown Palaeognathae.

---

### 2.2. A2 — `ENAM` Codon 1003 (`XP_019346377.1:p.Trp1003Ter`, Exon 8 `nt 2407–2409`) — `CORRIGIDO`

| Parameter | Empirical Value / Reconstructed State |
|---|---|
| Reference protein & 1-based residue | `XP_019346377.1` residue `1003` = **`W` (`Trp`, Tryptophan)** |
| *Alligator mississippiensis* Exon 8 coordinate | `nt 2407–2409` (1-based within Exon 8) |
| Crocodylia observed codons (`Allmis`, `Cropor`, `Gavgan`) | `TGG`, `TGG`, `TGG` (**Trp** in all 3 crocodylians) |
| Palaeognathae observed codons (`Strcam`, `Aptrow`, `Dronov`) | `Strcam = GGG` (**Gly**), `Aptrow = GGC` (**Gly**), `Dronov = GGC` (**Gly**) |
| Galloanserae observed codons (`Galgal`, `Anapla`) | `Galgal = TAC` (**Tyr**), `Anapla = ---` (uncovered) |
| Neoaves / Psittacopasserae observed codons (`Strhab`, `Taegut`) | `Strhab = TGA` (**Ter**), `Taegut = TGA` (**Ter**) |
| Per-nucleotide ML marginals at crown Neornithes (`Node1`, MAFFT) | Position 1: `p(T) = 0.6961, p(G) = 0.2899`; Position 2: `p(G) = 0.9034, p(A) = 0.0890`; Position 3: `p(C) = 0.5459, p(A) = 0.4624, p(G) = 0.2669` |
| Composite site-independent codon probabilities at crown Neornithes | `TGC` (**Cys**): $P = 0.3433$ (argmax); `TGA` (**Ter**): $P = 0.2907$ (total $P_{\text{stop}} = 0.2941$); `TGG` (**Trp**): $P = 0.1678$; `GGC` (**Gly**): $P = 0.1430$ |
| Crown Neognathae (`Node5`) & Psittacopasserae (`Node7`) ML $P_{\text{stop}}$ | Neognathae (`Node5`): $P_{\text{stop}} = 0.4498$ (`TGA` $P = 0.4366$, `TAC` $P = 0.3129$); Psittacopasserae (`Node7`): $P_{\text{stop}} = 0.8947$ (`TGA`) |
| 3-step nucleotide Fitch parsimony across the 3 codon positions | Pos 1: `{G, T}` (1 step); Pos 2: `G` (`G→A` in `Galgal`, 1 step); Pos 3: `{A, C, G}` (2 steps). Total = 4 nucleotide substitutions; `TGG (Trp) → TGA (Ter)` (1-step transition) and `TGC (Cys) → TGA (Ter)` (1-step transversion) are equally parsimonious |

**Resolution:** Under nucleotide-level `GTR+G` ML ASR, each codon position is reconstructed independently: because position 1 favors `T` (shared by Crocodylia `TGG`, `Galgal` `TAC`, and `Strhab`/`Taegut` `TGA`), position 2 favors `G` (shared by Crocodylia `TGG`, Palaeognathae `GGG`/`GGC`, and `Strhab`/`Taegut` `TGA`), and position 3 slightly favors `C` (`Aptrow`/`Dronov` `GGC` + `Galgal` `TAC`) over `A` (`Strhab`/`Taegut` `TGA`) and `G` (Crocodylia `TGG` + `Strcam` `GGG`), the product of the three independent marginals peaks at `TGC` (**Cys**, $P = 0.3433$), even though **zero extant birds carry `TGC`**. Calling the mutation simply ``p.Trp1003Ter`, avian `TGC→TGA`` obscured both the `TGG` (**Trp**) reference state and the site-independent ML reconstruction uncertainty. We corrected `Table_S5a`, Table 3, Figure 3B (`ref p.Trp1003Ter; avian TGG/TGC→TGA`), Section 2.5 (#4), and Section 3.4 to report this explicitly.

---

## 3. Complete Audit of All 23 Shared Coding Variants in `Table_S5a` (`CORRIGIDO`)

All 23 dual-aligner (`MAFFT L-INS-i` + `MUSCLE 5`) shared variants were audited for exact reference amino-acid identity, crocodylian outgroup codon consensus, crown-Neornithes ML ancestral codon and amino acid, and inferred substitution pathway:

| # | Dataset | Coordinate (`standardized_variant_id`) | *Alligator* Ref (`amis_codon` / AA) | Crocodylia Outgroups | Neornithes ML (`codon` / AA / $P_{\text{best}}$) | Carrier Codon(s) | Carriers (6-letter) | Inferred Substitution Pathway |
|---|---|---|---|---|---|---|---|---|
| 1 | `ENAM_exon8_aves` | `XP_019346377.1:p.Gly270Ter` | `GGA` (`Gly`) | `GGA` (`Gly`) | `GGA` (`Gly`, $0.590$) | `TGA` | `Aptrow,Dronov` | `GGA(Gly)->TGA(Ter)` |
| 2 | `ENAM_exon8_aves` | `XP_019346377.1:p.Pro311Ter` | `CCG` (`Pro`) | `CCA/CCG` (`Pro`) | `CAA` (`Gln`, $0.542$) | `TGA` | `Dronov,Strcam` | `CCG(Pro,Amis)->CAA(Gln,Neorn;P=0.542)->CGA(Arg,Palaeo;P=0.541)/TGA(Ter;P=0.453)` |
| 3 | `ENAM_exon8_aves` | `XP_019346377.1:p.Gln370Ter` | `CAA` (`Gln`) | `CAA` (`Gln`) | `CAA` (`Gln`, $0.683$) | `TGA` | `Aptrow,Dronov,Strcam` | `CAA(Gln)->TGA(Ter)` |
| 4 | `ENAM_exon8_aves` | `XP_019346377.1:p.Ser401Ter` | `TCC` (`Ser`) | `TCC` (`Ser`) | `TCC` (`Ser`, $0.885$) | `TGA` | `Aptrow,Dronov` | `TCC(Ser)->TGA(Ter)` |
| 5 | `ENAM_exon8_aves` | `XP_019346377.1:p.Ala439Ter` | `GCA` (`Ala`) | `CCA/GCA` (`Pro/Ala`) | `TCA` (`Ser`, $0.320$) | `TAG,TGA` | `Anapla,Aptrow,Dronov` | `GCA(Ala,Amis)->TCA(Ser,Neorn;P=0.320)->TAG,TGA(Ter)` |
| 6 | `ENAM_exon8_aves` | `XP_019346377.1:p.Gly558Ter` | `GGA` (`Gly`) | `GGA` (`Gly`) | `GGA` (`Gly`, $0.983$) | `TGA` | `Anapla,Galgal` | `GGA(Gly)->TGA(Ter)` |
| 7 | `ENAM_exon8_aves` | `XP_019346377.1:p.Ser559Ter` | `AGC` (`Ser`) | `AGC` (`Ser`) | `TGT` (`Cys`, $0.787$) | `TGA` | `Aptrow,Dronov` | `AGC(Ser,Amis)->TGT(Cys,Neorn;P=0.787)->TGA(Ter)` |
| 8 | `ENAM_exon8_aves` | `XP_019346377.1:p.Trp562Ter` | `TGG` (`Trp`) | `TGG` (`Trp`) | `TGG` (`Trp`, $0.965$) | `TGA` | `Strcam,Strhab` | `TGG(Trp)->TGA(Ter)` |
| 9 | `ENAM_exon8_aves` | `XP_019346377.1:p.Phe601Ter` | `TTT` (`Phe`) | `TTT` (`Phe`) | `TAT` (`Tyr`, $0.870$) | `TAA` | `Aptrow,Dronov,Strcam` | `TTT(Phe,Amis)->TAT(Tyr,Neorn;P=0.870)->TAA(Ter)` |
| 10 | `ENAM_exon8_aves` | `XP_019346377.1:p.Trp615Ter` | `TGG` (`Trp`) | `TGG` (`Trp`) | `TGG` (`Trp`, $0.695$) | `TAA,TGA` | `Dronov,Strcam,Strhab` | `TGG(Trp)->TAA,TGA(Ter)` |
| 11 | `ENAM_exon8_aves` | `XP_019346377.1:codon_619_ins2nt_fs` | Codon 619 (`Pro`) | In-frame | `+3 nt` in-frame in Aves | `+2 nt` (net 1-nt del) | `Anapla,Galgal` | `codon_619(Pro,Amis)_ins2nt_fs` |
| 12 | `ENAM_exon8_aves` | `XP_019346377.1:p.Gln667Ter` | `CAG` (`Gln`) | `CAG` (`Gln`) | `CAG` (`Gln`, $0.941$) | `TAG` | `Aptrow,Dronov` | `CAG(Gln)->TAG(Ter)` |
| 13 | `ENAM_exon8_aves` | `XP_019346377.1:p.Gln684Ter` | `CAA` (`Gln`) | `CAA` (`Gln`) | `CAA` (`Gln`, $0.548$) | `TAA,TAG` | `Aptrow,Dronov,Galgal,Strcam` | `CAA(Gln)->TAA,TAG(Ter)` |
| 14 | `ENAM_exon8_aves` | `XP_019346377.1:p.Trp751Ter` | `TGG` (`Trp`) | `TGG` (`Trp`) | `TGG` (`Trp`, $0.776$) | `TGA` | `Aptrow,Dronov,Strcam` | `TGG(Trp)->TGA(Ter)` |
| 15 | `ENAM_exon8_aves` | `XP_019346377.1:p.His844Ter` | `CAC` (`His`) | `CAC` (`His`) | `CAA` (`Gln`, $0.886$) | `TGA` | `Aptrow,Dronov` | `CAC(His,Amis)->CAA(Gln,Neorn;P=0.886)->TGA(Ter)` |
| 16 | `ENAM_exon8_aves` | `XP_019346377.1:p.Trp1003Ter` | `TGG` (`Trp`) | `TGG` (`Trp`) | `TGC` (`Cys`, $0.343$) | `TGA` | `Strhab,Taegut` | `TGG(Trp,Amis;Neorn_P=0.168)/TGC(Cys,Neorn_ML_argmax_P=0.343;no_extant_TGC)->TGA(Ter,Psittaco_P=0.892)` |
| 17 | `ENAM_exon8_aves` | `XP_019346377.1:p.Ser1022Ter` | `TCA` (`Ser`) | `TCA` (`Ser`) | `TCG` (`Ser`, $0.652$) | `TAG` | `Galgal,Strhab` | `TCA(Ser,Amis)->TCG(Ser,Neorn;P=0.652)->TAG(Ter)` |
| 18 | `ENAM_exon8_turtles` | `XP_019346377.1:p.Gln753Ter` | `CAG` (`Gln`) | `CAG` (`Gln`) | — | `TAG` | `Chemyd,Chrpic` | `CAG(Gln,Amis)->TAG(Ter)` |
| 19 | `ENAM_exon8_turtles` | `XP_019346377.1:p.Ala859Ter` | `GCA` (`Ala`) | `ACA/GCA` (`Thr/Ala`) | — | `TGA` | `Chemyd,Chrpic` | `GCA(Ala,Amis)->TGA(Ter)` |
| 20 | `AMBN_exon3_aves` | `XP_019346386.1:p.Ser56Ter` | `TCA` (`Ser`) | `TCA` (`Ser`) | — | `TAA` | `Aptrow,Dronov` | `TCA(Ser,Amis)->TAA(Ter)` |
| 21 | `AMBN_exon3_aves` | `XP_019346386.1:p.Gln62Ter` | `CAA` (`Gln`) | `CAA/CGA` (`Gln/Arg`) | — | `TGA` | `Aptrow,Dronov` | `CAA(Gln,Amis)->TGA(Ter)` |
| 22 | `AMBN_exon4_aves` | `XP_019346386.1:codon_109_del16nt_fs` | Codon 109 (`Pro`) | Intact (0 indels) | — | `16-nt del` | `Aptrow,Dronov,Strcam` | `codon_109(Pro,Amis)_del16nt_fs` |
| 23 | `AMBN_exon4_turtles` | `XP_019346386.1:codon_123_ins1nt_fs` | Codon 123 (`Pro`) | Intact (0 indels) | — | `1-nt ins` | `Chemyd,Chrpic` | `codon_123(Pro,Amis)_ins1nt_fs` |

---

## 4. Summary of Files Updated

1. [`scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py`](../../scripts/enamel_relics/12_audit_v102_variants_phylo_sim.py): Programmatic verification of all 23 variants against `XP_019346377.1` and `XP_019346386.1` and generation of `Table_S5a_level4_shared_mutations_v102.tsv`.
2. [`submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv`](../../submission_biorxiv/supplementary/Table_S5a_level4_shared_mutations.tsv): Added 12 v1.0.2 columns including `ref_aa_3letter`, `croc_consensus_aa_3letter`, `neornithes_ml_aa_3letter`, `neornithes_ml_p_best`, `neognathae_ml_p_stop`, and `inferred_substitution_path`.
3. [`scripts/enamel_relics/06_make_figures.py`](../../scripts/enamel_relics/06_make_figures.py) & [`submission_biorxiv/figures/fig3_lesions_and_phylogeny.png`](../../submission_biorxiv/figures/fig3_lesions_and_phylogeny.png): Updated Figure 3B (`ref p.Trp1003Ter; avian TGG/TGC→TGA`) and Figure 3C x-axis labels (`Cod 311\n[CAA: Gln]`, `Cod 439\n[TCA: Ser]`, `Cod 562\n[TGG: Trp]`, `Cod 615\n[TGG: Trp]`, `Cod 684\n[CAA: Gln]`, `Cod 1022\n[TCG: Ser]`).
4. [`submission_biorxiv/manuscript.typ`](../../submission_biorxiv/manuscript.typ), [`submission_biorxiv/manuscript.md`](../../submission_biorxiv/manuscript.md), and [`submission_biorxiv/manuscript_biorxiv.pdf`](../../submission_biorxiv/manuscript_biorxiv.pdf): Updated Section 2.5, Table 3, Figure 3 caption, and Section 3.4.
5. [`tests/test_enamel_relics_audit.py`](../../tests/test_enamel_relics_audit.py): Added automated unit tests asserting standard genetic-code translation consistency across all 23 rows of `Table_S5a`.
