# Phases D & E Audit Report — Calibration of Scientific Claims & Editorial Finalization of the 21-Page PDF

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final`  
**Manuscript Files:** [`submission_biorxiv/manuscript.typ`](../../submission_biorxiv/manuscript.typ) · [`submission_biorxiv/manuscript.md`](../../submission_biorxiv/manuscript.md) · [`submission_biorxiv/manuscript_biorxiv.pdf`](../../submission_biorxiv/manuscript_biorxiv.pdf)

---

## 1. Phase D — Line-by-Line Calibration of Overly Conclusive Expressions

We audited every occurrence of `"prove"`, `"proving"`, `"demonstrate"`, `"confirm"`, `"rule out"`, `"ruling out"`, `"eliminate"`, `"unequivocally"`, `"absence"`, and `"independently favors"` across the manuscript, calibrating each statement to match the exact reach of its underlying empirical evidence:

| Section & Topic | Prior Phrasing (`v1.0.2`) | Calibrated Phrasing (`audit/biorxiv-pre-submission-final`) | Scientific Rationale |
|---|---|---|---|
| **1. Introduction (#1) & 3.1 (#1)** — Dual-database annotation audit | `"proving that all nine odontogenic genes remain completely unannotated"` / `"dual-database proof"` / `"forensic native-CDS proof"` | `"showing that all nine odontogenic genes remain completely unannotated"` / `"dual-database verification"` / `"forensic native-CDS verification"` | Database parsing verifies annotation status in `RefSeq 106` and `Ensembl 116` without requiring mathematical `"proof"` terminology. |
| **1. Introduction (#4)** — Indels/splice sites vs. PTCs | `"demonstrating why multi-nucleotide frameshift indels and canonical splice-junction disruptions outperform isolated PTCs"` | `"showing why multi-nucleotide frameshift indels and canonical splice-junction disruptions are less prone to homoplasy than isolated PTCs"` | Precisely specifies the empirical property evaluated (lower homoplasy/reversal ambiguity). |
| **2.2. Reciprocal `BLASTx`** — Ratite *AMBN* top-HSP | `"top-HSP exon reciprocal BLASTx unequivocally confirms Alligator ameloblastin as the proteome-wide best hit"` | `"top-HSP exon reciprocal BLASTx identifies Alligator ameloblastin as the proteome-wide best hit"` | Reports the quantitative proteome-wide best hit (`51.2` and `50.8 bits`, $E \le 9.5\times 10^{-9}$) without rhetorical intensifiers. |
| **2.2 (#2). Chicken *ΨAMEL* (`EU340348`)** | `"proving that absence of tBLASTn signal cannot rule out a retained pseudogene relic"` | `"showing that absence of tBLASTn signal does not exclude a retained pseudogene relic"` | Replaces `"proving ... cannot rule out"` with precise methodological language. |
| **2.3 (`C1`). `Strhab` `ENAM_Acar` permutation score** | `"while confirming Strhab ENAM_Acar at 0.80 bits, empirical p = 0.024"` | `"while leaving Strhab ENAM_Acar at a marginal 0.80-bit excess, unadjusted exploratory empirical p = 0.024 across the 504-cell grid without multiplicity correction; orthology of Strhab ENAM is established independently by whole-proteome reciprocal BLASTx in Table S2b and 866-aa miniprot alignment in Table S4 rather than by this exploratory score"` | Prevents an unadjusted exploratory $p = 0.024$ in a 504-test family from being read as confirmatory proof. |
| **2.3. Equidistant *Anolis* query retention index** | `"Switching to equidistant Anolis carolinensis queries (~281 Ma...) eliminates this overlap cleanly:"` | `"Switching to equidistant Anolis carolinensis queries (~281 Ma...) resolves this query-distance overlap in the sampled panel:"` | Limits the statement to the 18-species sampled panel. |
| **2.3 & 3.3 (`D1`). Flanking globular *LPL* control** | `"This directly demonstrates clade-specific rate heterogeneity despite equidistant divergence—and rules out generalized avian substitution-rate acceleration or locus-wide alignment breakdown across NRG1–LPL."` / `"ruling out uniform avian rate acceleration across NRG1–LPL."` | `"This indicates clade- and protein-class-specific rate heterogeneity despite equidistant divergence and argues against uniform locus-wide alignment breakdown across the NRG1–LPL interval, although a single conserved globular control cannot exclude broader lineage- or gene-specific evolutionary rate variation."` / `"arguing against uniform locus-wide alignment failure across NRG1–LPL."` | A single conserved globular enzyme (*LPL*) shows that the syntenic neighborhood aligns well and does not suffer locus-wide collapse in birds, but one control gene cannot rule out all forms of avian substitution-rate variation. |
| **2.4. Level 3 `miniprot` *LPL* control & *C. picta* `protein SON` artifact** | `"ruling out systematic miniprot alignment artifacts or assembly-quality confounding"` / `"which whole-proteome reciprocal BLASTx proved was protein SON"` | `"indicating that our miniprot parameters do not generate spurious lesions on conserved globular coding sequence within the exact same genomic intervals"` / `"which whole-proteome reciprocal BLASTx matched to protein SON (E = 0.010) rather than enamelin"` | Matches the scope of the positive control and reciprocal `BLASTx` result. |
| **2.5 (`C3`). `ENAM` Exon 8 RELL topology tests** | `"Alongside its external phylogenomic support, ENAM Exon 8 independently favors T1... formally rejecting T2... and T3"` | `"Within this single-locus alignment, ENAM Exon 8 exhibits a higher log-likelihood under T1... than under T2... or T3...; however, because these RELL tests are computed on the same pseudogene alignment used to score PTC congruence—and because incomplete lineage sorting is well documented at the base of Palaeognathae (Sackton et al., 2019)—we anchor T1 on external phylogenomic evidence and evaluate PTC incongruence across all three topologies."` | Avoids treating a single-exon likelihood test on the focal alignment as an independent species-tree validation. |

---

## 2. Preservation of Validated Scientific Distinctions (`D2`, `D3`, `D4`)

1. **D2 — Graded Orthology Evidence Tiers (`CONFIRMADO`):**
   - Chicken *ΨENAM* Exon 8 (`NC_052572.1:55,155,480–55,157,464`, `−`, 1,985-bp union): whole-proteome reciprocal `BLASTx`-confirmed (`XP_019346377.1`, `41.6 bits`, $E = 2.7\times 10^{-5}$).
   - Chicken *ΨENAM* Exon 6 (`55,160,529–55,160,677`, `−`, 149-bp flanked span; 47-bp exon core `55,160,570–55,160,616`): nucleotide `BLASTn` ($E = 5.3\times 10^{-11}$ vs. *Alligator*; $E = 2.7\times 10^{-27}$ vs. *A. rowi*), canonical `AG...GT` splice borders, and synteny confirmed; **not** reciprocal `BLASTx`-confirmed ($E = 0.29$, truncated at codon 4 by `CAA→TAA` + 1-nt deletion).
   - Chicken *ΨAMBN* Exon 3 (`55,168,855–55,168,992`, `−`, 138 bp): whole-proteome reciprocal `BLASTx`-confirmed (`XP_019346384.1`, `45.8 bits`, $E = 5.1\times 10^{-7}$).
   - Kakapo *ΨAMTN* (`NC_044302.2:51,216,319–51,216,444`, `−`, 126 bp): whole-proteome reciprocal `BLASTx`-confirmed (`XP_014451634.1`, `38.1 bits`, $E = 1.9\times 10^{-4}$).
2. **D3 — Gnomon `LOW QUALITY PROTEIN` `PROTEIN_CODING` Models (`CONFIRMADO`, `Table_S7`):**
   - Strictly preserves the distinction between the **three chromosome-level confirmed annotation-compensated pseudogenes** (*C. picta* *ENAM*, *C. mydas* *AMBN*, *S. camelus* *AMELX*) and the **one scaffold-level single-assembly lesion** (*A. rowi* *AMELX*, `NW_020450766.1`, contig N50 129.4 kb, codon-177 `TAA` masked to `X`, pending a second conspecific assembly).
3. **D4 — Molecular Variant Nomenclature (`CONFIRMADO`, `Table_S5a`):**
   - Strictly preserves the separation between *Alligator* reference residues (`XP_019346377.1:p.Phe601Ter` [`TTT` Phe], `p.Trp1003Ter` [`TGG` Trp], `p.Ser1022Ter` [`TCA` Ser]) and reconstructed crown-Neornithes ancestral states (`TAT` [Tyr, $P = 0.870–0.899$] $\to$ Palaeognathae `TAA` at codon 601; site-independent Neornithes ML argmax `TGC` [Cys, $P = 0.343$] vs. `TGG` [$P = 0.168$] $\to$ Psittacopasserae `TGA` at codon 1003; Neornithes `TCG` [Ser, $P = 0.652$] $\to$ equiparsimonious 2-step `Galgal`+`Strhab` `TAG` vs. reversible ML Neognathae `TAG` $P_{\text{stop}} = 0.926$ with `Taegut` `AAA` and `Anapla` uncovered at codon 1022).

---

## 3. Phase E — Final Editorial Verification of the 21-Page PDF

- **E1 (Abstract):** Verified at exactly **248 words** ($\le 250$-word bioRxiv limit), identical normalized text between `manuscript.typ` and `manuscript.md`, clear separation between empirical findings and coverage/simulation limitations.
- **E2 (Figures 1–3 & Tables 1–3):** Visually inspected via 150-DPI full-page rasterization of all modified pages (`Page 1`, `Pages 10–19`) and direct inspection of `fig3_lesions_and_phylogeny.png` (`320 DPI`, Okabe-Ito palette). Figure 3B now explicitly displays `"Crown Neornithes: none detected (0/878 7-bird cov.; T1: 5–6/16 incong.; T2/T3: 9–11/16)"`.
- **E3 (Bibliography):** All `49/49` BibTeX entries in `submission_biorxiv/references.bib` match cited keys `1:1` in `manuscript.typ` (`0` uncited, `0` missing, `0` retractions).
- **E4 (Compiled PDF):** Recompiled with `Typst 0.15.1` (`21` pages, continuous line numbers `1–679`, searchable vector text and embedded fonts, zero broken links or placeholders).
