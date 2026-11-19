# Change Log: `enamel-relics-v1.0.0` → `audit/biorxiv-revision-v1.0.1`

**Repository:** `https://github.com/patoangolano/bio-stack`  
**Historical Baseline Tag (Untouched):** `enamel-relics-v1.0.0` (tag object `b9ff25e654007ef26f5ce01ce702f68c214c90b2` → commit `a425b2d1b20bd0fd77a997915d62523f5b8f2649`)  
**Branch Base (`main`):** `99427639b778e341ad6bcbef32fe1675ce251616` (post-release PDF-only DOI update)  
**Revision Branch:** `audit/biorxiv-revision-v1.0.1`  
**Audit Date:** 2026-11-19

---

## 1. Scientific & Genomic Corrections (Phases B1–B3)

### 1.1. Chicken *ΨENAM* Exon 6 Evidence-Tier Correction (Phase B1)
- **Issue in `v1.0.0`:** Abstract, Introduction, and Figure 1A implied whole-proteome reciprocal `BLASTx` confirmation for all three chicken chrZ relics, even though Table 1 only reported `BLASTn` for *ΨENAM* Exon 6 and `Table_S2b` omitted it.
- **Empirical Audit (`scripts/enamel_relics/07_audit_phaseA_B1_B3.py`):** Ran whole-proteome `BLASTx` against all 48,040 *Alligator mississippiensis* proteins (`-seg yes`, `-seg no`, and ultra-sensitive `-word_size 2 -evalue 1000 -comp_based_stats 0`) on the 149-bp `BLASTn` span (`NC_052572.1:55,160,529–55,160,677`, `−`), the exact 47-bp coding core (`55,160,570–55,160,616`, `−`), `±30 bp`, and `±100 bp` windows, plus an intact 48-bp *Alligator* Exon 6 positive control.
- **Result & Revision:** Chicken *ΨENAM* Exon 6 returned **0 enamelin hits** in whole-proteome `BLASTx` (top non-specific hit: `XP_019343727.1`, 28.5 bits, *E* = 0.29), whereas the intact 48-bp *Alligator* control passed (`XP_019346377.1`, 37.7 bits, *E* = 3.6×10⁻⁵). Direct `BLASTn` (`-word_size 7`) vs. *Alligator* Exon 6 (±60 bp) confirmed orthology at 72.6% identity over 135 bp (60.8 bits, *E* = 5.3×10⁻¹¹ vs. runner-up 30.1 bits, *E* = 0.089), alongside intra-avian `BLASTn` vs. *A. rowi* (115 bits, *E* = 2.7×10⁻²⁷), single-copy status (1 hit in `dc-megablast` and `blastn`), canonical `AG..GT` splice borders, and the exact in-frame `CAA→TAA` stop at codon 171 + 1-nt deletion (`GGA→-GA`).
- **Updated Files:** `Table_S2b`, `manuscript.typ`, `manuscript.md`, Table 1, Figure 1B, `audit/reproducibility_audit.md`.

### 1.2. Native-CDS Forensic Dissection of Four RefSeq `PROTEIN_CODING` Models (Phase B2)
- **Issue in `v1.0.0`:** Lumped *C. picta* *ENAM*, *C. mydas* *AMBN*, *S. camelus* *AMELX*, and *A. rowi* *AMELX* together as "demonstrably pseudogenized" false-positive RefSeq models even though `miniprot` had only aligned the two turtle loci and neither ratite *AMELX* native CDS had been inspected.
- **Empirical Audit (`scripts/enamel_relics/08_audit_protein_coding_models.py`, `08b`, `08c`):** Retrieved GenBank/GenPept records, native genomic CDS exons, intron dinucleotides, ±2-kb `N`-gaps, flanking synteny, and whole-proteome `BLASTp` orthology for all four models.
- **Result & Revision:** All four are genuine orthologs explicitly flagged `LOW QUALITY PROTEIN` by NCBI Gnomon with synthetic base patches and zero nearby `N`-gaps. Reclassified with forensic precision in new **`Table_S7_native_model_audit_v101.tsv`**, Abstract, Section 2.1, Section 3.2, and Section 3.5:
  1. `Chrpic` *ENAM* (`XM_065550579.1`): 761-bp native CDS, `+1N` synthetic insertion at cod 33 (non-homopolymer) → **`CONFIRMED_ANNOTATION_COMPENSATED_PSEUDOGENE`**.
  2. `Chemyd` *AMBN* (`XM_037903130.1`): 1,201-bp native CDS, `+2N` at cod 246 + masked native genomic `TAA` PTC at cod 44 + skipped `GT→AT` Intron 8 donor → **`CONFIRMED_ANNOTATION_COMPENSATED_PSEUDOGENE`**.
  3. `Strcam` *AMELX* (`XM_068931113.1`): 572-bp native CDS, `+1N` at cod 95 + 3 non-homopolymer native PTCs (cod 72/92/106) + non-canonical Intron 1 (`CC..CC`) → **`CONFIRMED_ANNOTATION_COMPENSATED_PSEUDOGENE`**.
  4. `Aptrow` *AMELX* (`XM_026079913.1`, scaffold `NW_020450766.1`, contig N50 129.4 kb, 0 second conspecific assemblies on NCBI): 618-bp native CDS, 1 masked native `TAA` PTC at cod 177 (near 4-A run) + non-canonical Intron 2 (`CA..CC`) → conservatively reclassified as **`ANNOTATION_COMPENSATED_LESION_SINGLE_ASSEMBLY`**.

### 1.3. Coordinate & Indexing Convention Reconciliation (Phase B3)
- **Issue in `v1.0.0`:** Mixed endpoint difference (`105,218 bp`) with strict intergenic base count (`105,217 bp`), showed a `138-kb` plotted window in Figure 1B without explaining the flanking genes included, and cited `1,975 bp` alongside coordinates `55,155,480–55,157,464` (`1,985 bp`).
- **Revision (`07_audit_phaseA_B1_B3.py`, `11_consolidate_v101_tables.py`):** Standardized all genomic coordinates to **1-based inclusive** throughout text, Table 1, and Figure 1, added explicit `strict_intergenic_bp`, `intergenic_1based_inclusive`, and `intergenic_bed_0based` columns to `Table_S2a`, explicitly documented the `138-kb` plotted coordinate window (`55.098–55.236 Mb`) in Figure 1B, and reconciled *ΨENAM* Exon 8 (`1,646-bp` `tBLASTn` core `55,155,480–55,157,125` ∪ `1,975-bp` intra-avian `BLASTn` span `55,155,490–55,157,464` = `1,985-bp` genomic union).

---

## 2. Statistical Re-Evaluation (Phase C)

- **Issue in `v1.0.0`:** Estimated 95th-percentile null thresholds from only `N = 20` composition-preserving permutations (high quantile variance; extreme `25.87-bit` collapse on heavily `SEG`-masked `DSPP_Amis`), lacked uncertainty intervals, lacked local block-preserving controls, and did not explicitly separate equal divergence time from clade-specific substitution rates.
- **Upgrade (`scripts/enamel_relics/09_permutation_sensitivity.py`):** Executed `N = 20, 200, 500, 1,000` global permutations (seed `20261119`), `N = 200` local 20-aa block-preserving permutations (seed `20261120`), and 2,000-replicate bootstrap 95% CIs (seed `20261121`) across all `504` species × query × `SEG` cells (reproducing `504/504` observed bitscores identically).
- **Updated Quantitative Findings:**
  - Median bootstrap 95% CI width on $p_{95}$ narrowed **3.4-fold** (`2.70 → 0.80 bits`; `DSPP_Amis` `25.87 → 0.70 bits`).
  - Twelve marginal noise-floor calls (`0.26–1.42 bits`) flipped cleanly (notably `Pelsin` and `Tingut` `ENAM_Amis` dropped to `0.00` excess).
  - Equidistant *Anolis* medians updated to `N = 1,000` across Table 2, Figure 2, Abstract, and text: `LPL` Aves/Croc = `112.6%` (`135.6` vs. `120.4`; top-HSP aa identity `88.7%` Aves vs. `80.0%` Croc vs. `73.9%` Turt, proving rate heterogeneity despite equal ~281 Ma divergence); `IBSP` `74.5%` (`-seg yes`) / `48.5%` (`-seg no`, `51.5%` drop); `SPP1` `50.0%` (`16.8`, range `7.3–23.9 bits`, overlapping avian `ENAM` `0.0–12.7` and `AMBN` `0.0–15.4` bits; robust under 20-aa block shuffles).
  - Quantified query length ($\rho = +0.651/+0.702$), entropy ($\rho = +0.769$), and `SEG` fraction ($r = -0.795/-0.670$) scaling of null baselines.
  - Explicitly avoided invalid pooled-species inferential *p*-values across *n* = 2 evolutionary origins of edentulism (Aves, Testudines).

---

## 3. Mutation, Splice-Border, Nomenclature & Phylogenetic ASR Audit (Phase D)

- **Script (`scripts/enamel_relics/10_audit_mutations_and_asr.py`):**
  1. **Unambiguous Species Identifiers:** Created **`Table_S0_species_identifiers.tsv`** replacing ambiguous `Tgut` / `Tgua` with unambiguous 6-letter codes (`Taegut` = *Taeniopygia guttata*; `Tingut` = *Tinamus guttatus*; etc.) across all tables, figures, and manuscript text.
  2. **Standardized Variant Nomenclature:** Replaced incorrect HGVS `c.` codon prefixes (`c.311`, `c.619`, etc.) with `XP_019346377.1:p.Xxx###Ter` / `codon_###` and `XP_019346386.1:codon_###`.
  3. **IQ-TREE 3.1.4 ML Ancestral State Reconstruction (`GTR+G --ancestral`) & Topology Tests:** Tested all three Palaeognathae resolutions (`T1` Novaeratitae $\ln L = -15,237.86$ strongly preferred over `T2` $\Delta\ln L = -22.66$ and `T3` $\Delta\ln L = -22.86$). Corrected overstatements at:
     - `p.Pro311Ter` (`Allmis` `CCG` → Neornithes `CAA` → Palaeognathae `CGA` CpG codon $P = 0.541$ vs. `TGA` $P = 0.453$: genuine CpG context within Palaeognathae, ambiguous between convergence and `Aptrow` reversal).
     - `p.Ser1022Ter` (reversible ML ASR reconstructs crown Neognathae `TAG` at $P_{\text{stop}} = 0.926$ with secondary passerine decay to `AAA` in `Taegut`, giving `5/16 = 31.2%` [95% CI `11.0–58.7%`] under ML ASR vs. `6/16 = 37.5%` [95% CI `15.2–64.6%`] under irreversible parsimony).
     - `p.Gln684Ter` (`TAA↔TAG`) and `p.Trp615Ter` (`TAA↔TGA`): explicitly noted 1-transition drift caveat and grounded independent-origin inference in phylogeny + ML ASR rather than stop-triplet differences alone.
  4. **Neutral Stop-Gain Simulation (`1,000` replicates, K80 $\kappa = 4$, seed `20261122`):** Demonstrated expected median `26–29` shared PTCs, `48.1–51.8%` incongruence (95% range `29.2–70.8%`), and median `0` crown-Neornithes PTCs.
  5. **Flanked-MSA Splice-Border Re-Audit (`Table_S5b`):** Replaced raw local-BLASTn endpoint slicing with true ±60–80 bp intronic flanks aligned in both MAFFT and MUSCLE, removing invalid terminal *AMBN* Exon 9 5' donor rows and resolving *AMBN* Exon 3 (`AG..GT` in 4 birds; chicken-derived `AC..GT`) and *AMBN* Exon 4 (Durocryptodira `AC..GC`).

---

## 4. Bibliographic Integrity (Phase E)

- **Audit (`audit/bibliography_audit.py`, `audit/phaseE_complete_check.py`, `audit/bibliography_report.md`):**
  - Verified **46/46** original entries via Crossref (`HTTP 200`), Europe PMC (`44/46`), and Retraction Watch (`0` retractions/concerns).
  - Corrected missing consortium author (`and {Zoonomia Consortium}`) in `kirilenko2023`.
  - Added **3** newly verified primary references supporting v1.0.1 analyses: `wong2026` (IQ-TREE 3, *MBE* 43(5):msag117), `sackton2019` (*Science* 364:74–78), and `clopper1934` (*Biometrika* 26(4):404–413).
  - Updated `submission_biorxiv/references.bib` and `runs/enamel_relics/references.bib` to **49/49 verified entries** with 100% 1-to-1 text-bibliography correspondence.

---

## 5. Policy, Affiliation, Zenodo Metadata & Editorial Fixes (Phases A, F, G, H)

- **bioRxiv / openRxiv Readiness (`audit/biorxiv_readiness_report.md`):** Archived primary openRxiv policies (20 Oct 2026 AI & Accountability policy; 29 May 2026 Trust Signals), assessed affiliation screening risk, drafted conservative Option A (compiled into `manuscript.typ`/`manuscript.md` with factual Senac EAD student *Author Note*, zero unverified institutional claims) + Option B + formal inquiry emails to Senac (PT) and `bioRxiv@openrxiv.org` (EN).
- **Zenodo & Citation Metadata:** Diagnosed why tag `enamel-relics-v1.0.0` (`10.5281/zenodo.24358067`) carried provisional `XXXXXXX` placeholders and GitHub username `patoangolano` as creator (`audit/reproducibility_audit.md`). Preserved tag `enamel-relics-v1.0.0` untouched, cited Zenodo **Concept DOI `10.5281/zenodo.24358066`** consistently across `manuscript.typ`, `manuscript.md`, and `README_SUBMISSION.md`, and added root `.zenodo.json` and `CITATION.cff` crediting `Lima, Matheus` (`ORCID 0009-0000-3557-8337`).
- **Supplementary Tables & Figures (`11_consolidate_v101_tables.py`, `06_make_figures.py`):**
  - Consolidated `Table_S0`–`Table_S7` (12 TSVs + `MANIFEST_SHA256.tsv`), purged all Portuguese strings from `Table_S6`, embedded numeric `odaph_gap_bp`, and mirrored identically between `runs/enamel_relics/supplementary/` and `submission_biorxiv/supplementary/`.
  - Decoupled `06_make_figures.py` completely from legacy `runs/atavismo/` (100% driven by `supplementary/Table_S*.tsv`), enlarged labels, fixed overlap artifacts, and re-rendered `.png`/`.pdf`/`.svg` across all three figures.
- **Manuscript Typography & Compilation:** Fixed running-header spacing, Typst `#super[...]\;` semicolon escapes, unescaped `(~1.16 GB)` tilde, `_LPL_-normalized` italics, double-parenthesis `#pcite` formatting, and floating figure placement (`placement: auto`) → clean 19-page PDF (`manuscript_biorxiv.pdf`, 245-word Abstract, continuous lines 1–557), mirrored verbatim into pure Markdown `manuscript.md`.
- **Automated Acceptance Test Suite (`tests/test_enamel_relics_audit.py`):** Added automated checks verifying SHA-256 manifests, coordinate/indexing invariants, species/sequence counts, Level 3 lesion rates and 133/133 `LPL` introns, Level 4 counts and splice borders, Table 1–3 numeric consistency with supplementary TSVs, citation completeness, zero placeholders/ambiguous codes, PDF–Typst–Markdown sync, and strict immutability of tag `enamel-relics-v1.0.0`.
