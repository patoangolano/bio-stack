# Phase A–D Technical, Genomic, and Reproducibility Audit Report

**Audit Date:** 2026-11-19  
**Audit Branch:** `audit/biorxiv-revision-v1.0.1`  
**Reference Archival Release (Unmodified):** `enamel-relics-v1.0.0` (annotated tag `b9ff25e654007ef26f5ce01ce702f68c214c90b2` → commit `a425b2d1b20bd0fd77a997915d62523f5b8f2649`)  
**Updated `main` Baseline Commit:** `99427639b778e341ad6bcbef32fe1675ce251616`  
**Zenodo Deposit Audited:** `https://doi.org/10.5281/zenodo.24358067` (Version DOI) / `https://doi.org/10.5281/zenodo.24358066` (Concept DOI)

---

## 1. Phase A — Repository Structure, Integrity, Release Comparison, and Zenodo Verification

### 1.1. SHA-256 & File Integrity Check (`v1.0.0` Baseline)
- **Supplementary Tables (`runs/enamel_relics/supplementary/` & `submission_biorxiv/supplementary/`):** **10/10** original TSV files matched `MANIFEST_SHA256.tsv` byte-for-byte (**PASSOU**).
- **Publication Figures (`runs/enamel_relics/figures/` vs. `submission_biorxiv/figures/`):** **9/9** files (`.png`, `.pdf`, `.svg` for Figures 1–3) matched SHA-256 identically (**PASSOU**).
- **Tool Versions Verified Live in WSL2 (`Ubuntu-26.04`, kernel `6.18.58.1-microsoft-standard-WSL2`):**
  - `tblastn` / `blastx` / `blastn`: `2.17.0+` (exact match to `tool_versions.txt`)
  - `segmasker`: `1.0.0` (exact match)
  - `diamond`: `2.2.8` (exact match)
  - `miniprot`: `0.18-r281` (exact match)
  - `MAFFT`: `v7.525 (2024/Mar/13)` (exact match)
  - `MUSCLE`: `5.3.linux64` (exact match)
  - `IQ-TREE 3`: `3.1.4` (live verified; used in Phase D ML ASR)
  - `samtools`: `1.24`, `seqkit`: `v2.14.0`, `NCBI Datasets CLI`: `18.38.0`, `Typst`: `0.15.1 (6734eafd)` (exact match)

### 1.2. Discrepancies Found in `v1.0.0` Reproducibility & Provenance
1. **Why `enamel-relics-v1.0.0` Still Displays Placeholder DOI `10.5281/zenodo.XXXXXXX`:**
   - Tag `enamel-relics-v1.0.0` points immutably to commit `a425b2d` (created `2026-11-19T02:56Z`), whose PDF (`SHA-256: f768fb366ce079653cdd61168b13c1e06a6832b5ef36189b60816c1d8131c092`), `manuscript.typ` (`a3505e74...`), `manuscript.md` (`45611011...`), and `README_SUBMISSION.md` contained `[Zenodo DOI: 10.5281/zenodo.XXXXXXX — insert upon GitHub release]`.
   - Publishing the GitHub Release (`2026-11-19T03:18:05Z`) triggered the Zenodo webhook (`2026-11-19T03:18:18Z`), which archived `patoangolano/bio-stack-enamel-relics-v1.0.0.zip` (`83,096,375` bytes, `MD5: 4943fe7d17410ae4abc570679933ba46`) and minted Version DOI `10.5281/zenodo.24358067` (Concept DOI `10.5281/zenodo.24358066`). Direct download and inspection of the Zenodo ZIP confirmed its embedded `manuscript_biorxiv.pdf` matches `f768fb366ce07965...` (placeholder DOI).
   - Commit `9942763` on `main` inserted `10.5281/zenodo.24358067` and recompiled `manuscript_biorxiv.pdf` (`SHA-256: 9aec8d70f68e...`) *after* the release, creating a self-referential version mismatch (the `main` PDF cited the `v1.0.0` version DOI whose own archive predated that commit).
   - **Zenodo & DataCite Metadata Verification (`audit/evidence/zenodo_record_24358067.json` & `datacite_24358067.json`):** DOI resolves via `https://doi.org/10.5281/zenodo.24358067` (`HTTP 302 → 302 → 200`), state `findable`, license `cc-by-4.0`, linked to `https://github.com/patoangolano/bio-stack/tree/enamel-relics-v1.0.0`. However, creator metadata is `[{"name": "patoangolano", "affiliation": null}]` with empty ORCID because `.zenodo.json` / `CITATION.cff` were missing.
   - **Fix in `v1.0.1`:** Added `.zenodo.json` and `CITATION.cff` crediting `Lima, Matheus` (`ORCID: 0009-0000-3557-8337`), preserved tag `enamel-relics-v1.0.0` untouched, and cited the **Concept DOI `10.5281/zenodo.24358066`** as the primary manuscript archival link (automatically resolving to the latest release).
2. **`tblastn` `-max_target_seqs` Methods Discrepancy:** `manuscript.typ` line 336 stated `-max_target_seqs 100`, whereas `02_run_alignments_wsl.sh` executed `-max_target_seqs 50` (functionally identical against single-sequence subject FASTAs, but fixed in Methods to match the exact CLI flag).
3. **Offline vs. Online Reconstruction Scope & Legacy Directory Dependency:**
   - Tables `S1–S5` and Figures `2–3` could be rebuilt offline from tracked files (`runs/enamel_relics/{tables,level4,queries,audit_metadata.json}`), whereas raw genomic intervals, proteomes, BLAST databases, and GFFs (`1.16 GB` total in `refdata/`, `blastdb/`, `intervals_*`, `tblastn_*`, `miniprot/`, `reciprocal/*.fna|tsv`) are intentionally `.gitignore`-excluded and require re-running Steps `01`–`02` against NCBI/Ensembl to reconstruct from scratch.
   - **Found & Fixed:** `05_consolidate_tables.py` and `06_make_figures.py` depended on `runs/atavismo/odam_odaph/s32_vaos.json` for Figure 1C, and `05_consolidate_tables.py` had a string-replacement order bug leaving Portuguese (`ODAM nao localizado em none especie COM dente; sensibilidade no ODAM do Xenopus:`) in all 18 `ODAM` rows of `Table_S6`. Both fixed in `11_consolidate_v101_tables.py` and the updated figure pipeline (Figure 1C now reads `odaph_gap_bp` directly from `Table_S6`).
4. **Provenance Log Gaps:** `provenance_run.jsonl` logged 12 events (`01`, `03`, `04`, `06`, Typst compiles) but lacked entries for `01b`, `02`, `05`, and commit `9942763`. All `v1.0.1` audit steps write dedicated provenance logs under `runs/enamel_relics/v101/provenance_*.jsonl`.
5. **WSL Distro Documentation Drift:** Root `AGENTS.md` documents `Ubuntu-24.04`, whereas the active runtime distro used by the pipeline is `Ubuntu-26.04` (correctly documented in `scripts/enamel_relics/README.md` and `manuscript.typ`).

---

## 2. Phase B — Genomic Evidence Corrections

### 2.1. Phase B1: *Gallus gallus* `GRCg7b` $\Psi$*ENAM* Exon 6 (`NC_052572.1:55,160,529–55,160,677`, `−`)
- **Script & Output:** `scripts/enamel_relics/07_audit_phaseA_B1_B3.py` → `runs/enamel_relics/v101/phaseAB/phase_A_B1_B3_audit.json`.
- **Critical Correction (FALHOU in `v1.0.0`, CORRIGIDO in `v1.0.1`):**
  - The `v1.0.0` Abstract, Introduction, and Discussion grouped $\Psi$*ENAM* Exon 6 with $\Psi$*ENAM* Exon 8 and $\Psi$*AMBN* Exon 3 as `"whole-proteome reciprocal-confirmed"`.
  - Direct empirical testing of chicken $\Psi$*ENAM* Exon 6 (`149-bp` BLASTn span `NC_052572.1:55,160,529–55,160,677`, `47-bp` exon core `NC_052572.1:55,160,570–55,160,616`, plus `±30 bp` and `±100 bp` windows) via `blastx` against all `48,040` *Alligator mississippiensis* proteins **FAILED** whole-proteome reciprocity under every tested parameter regime:
    - Standard (`-seg yes`, `-evalue 10`): top hit for 149-bp span is `XP_019343727.1` (neurochondrin, `28.5 bits`, $E = 0.29$); top hit for 47-bp core is `XP_019340150.1` (`23.1 bits`, $E = 4.8$); **0 enamelin hits**.
    - Unmasked (`-seg no`, `-evalue 10`): identical (`0` enamelin hits).
    - Ultra-sensitive (`-word_size 2 -seg no -evalue 1000 -comp_based_stats 0`): still **0 enamelin top hits** (while the intact 48-bp *Alligator* Exon 6 positive control returns `XP_019346377.1` enamelin at `37.7 bits`, $E = 3.61 \times 10^{-5}$ under standard `-seg yes`, and even intact-frame *Apteryx rowi* Exon 6 only ranks 3rd at `26.6 bits`, $E = 2.3$).
  - **Biophysical reason:** *Alligator* *ENAM* Exon 6 is only 48 bp (16 aa, `XP_019346377.1` codons `168–183`, phase 0), and chicken Exon 6 carries an in-frame `TAA` stop at the 4th codon (codon 171: avian ancestral `CAA → TAA`, shared convergently with *Struthio camelus*) immediately followed by a 1-nt frameshift deletion (`GGA → -GA`), leaving only 3 intact codons before truncation.
  - **Confirmed Evidence Level for $\Psi$*ENAM* Exon 6 (PASSOU):**
    1. **Direct *Alligator* nucleotide `BLASTn` (`-word_size 7`):** `Amis` Exon 6 (`±60 bp`, 168 bp) recovers the exact chicken locus as the sole high-scoring hit in the 690-kb `NRG1–LPL` window (`72.59%` identity over `135 bp`, `60.8 bits`, $E = 5.25 \times 10^{-11}$; runner-up `30.1 bits`, $E = 0.089$; bare 48-bp exon core alone hits at `35.6 bits`, $E = 4.23 \times 10^{-4}$).
    2. **Intra-avian `BLASTn` (`-word_size 7`):** *A. rowi* Exon 6 (`±60 bp`) hits at `78.00%` identity over `150` alignment columns (`149 bp` genomic span, `115 bits`, $E = 2.71 \times 10^{-27}$).
    3. **Single-copy in `NRG1–LPL`:** Exactly 1 hit in both `dc-megablast` and `blastn` (`100%`, `149 bp`, `269 bits`, $E = 4.87 \times 10^{-74}$).
    4. **Syntenic collinearity & splice borders:** Exact `−` strand, situated `3,064 bp` strictly upstream of $\Psi$*ENAM* Exon 8 (149-bp span) / `3,105 bp` (47-bp exon core `NC_052572.1:55,160,570–55,160,616`), retaining canonical `AG..GT` splice dinucleotides in both MAFFT and MUSCLE, with zero RefSeq 106 or Ensembl 116 overlap.
  - **Manuscript Fix:** Downgraded $\Psi$*ENAM* Exon 6 everywhere (Abstract, Intro, Results 2.2, Table 1, Figure 1B, Discussion 3.1) to **nucleotide- and synteny-confirmed** (explicitly noting its failure of whole-proteome `BLASTx` reciprocity due to short 47-bp length + early codon-4 stop + 1-nt frameshift).

### 2.2. Phase B2: Forensic Audit of the Four RefSeq `PROTEIN_CODING` Enamel Models
- **Scripts & Outputs:** `scripts/enamel_relics/08_audit_protein_coding_models.py`, `08b_b2_native_exons_and_replication.py`, `08c_b2_exact_ptc_and_tableS7.py` → `runs/enamel_relics/v101/b2/Table_S7_native_model_audit_v101.tsv`.
- **Key Discovery (PASSOU with nuanced reclassification):** Auditing the native GenBank/Gnomon models directly on the genomic sequence (without relying on inter-species alignment) proved that **all 4 models (`4/4`) are explicitly tagged `LOW QUALITY PROTEIN` by NCBI Gnomon**, which artificially inserted synthetic `N` bases and/or masked native genomic stop codons to `X` in `XP_` so it could emit `PROTEIN_CODING` records:

| Locus | Assembly & RefSeq Model | Orthology (Synteny + `BLASTp` vs. *Alligator*) | Native Genomic Lesions Patched by Gnomon | Annotated Introns | Assembly Quality / Gaps (`±2 kb`) | Revised Classification |
|---|---|---|---|---|---|---|
| **Chrpic *ENAM*** | `GCF_011386835.1` (`XM_065550579.1` / `XP_065406651.1`, 253 aa, 2 exons) | Flanked by `LPL`; `XP_019346375.1` `202 b` ($E=3.4\times 10^{-59}$) | **1 frameshift (`+1N` inserted)** at `NC_088796.1:129,121,140` (cod 33, non-homopolymer); native CDS = `761 bp`; `miniprot` confirms 12 fs + 4 PTCs across 713 aa | `1/1` canonical (`GT..AG`) | Chromosome-level (contig N50 `13.4 Mb`); `0` N-gaps | **CONFIRMED_ANNOTATION_COMPENSATED_PSEUDOGENE** |
| **Chemyd *AMBN*** | `GCF_015237465.2` (`XM_037903130.1` / `XP_037759058.1`, 400 aa, 6 exons) | Flanked by `NRG1`, `LPL`, `PSD3`; `XP_019346385.1` `252 b` ($E=4.7\times 10^{-80}$) | **1 frameshift (`+2N` inserted)** at `NC_051245.2:128,896,655` (cod 246, in 4-C run) + **1 native PTC (`TAA→X`)** at `128,902,406–128,902,408` (cod 44); native CDS = `1,201 bp`; `miniprot` confirms `GT→AT` Intron 8 donor | `5/5` annotated canonical (`GT/GC..AG`); `miniprot` shows skipped `GT→AT` donor | Chromosome-level (contig N50 `39.4 Mb`); `0` N-gaps | **CONFIRMED_ANNOTATION_COMPENSATED_PSEUDOGENE** |
| **Strcam *AMELX*** | `GCF_040807025.1` (`XM_068931113.1` / `XP_068787214.1`, 190 aa, 3 exons) | Adjacent to `ARHGAP6` (chr1); `XP_014463267.1` `93.6 b` (`-seg no`, $E=1.3\times 10^{-23}$) | **1 non-homopolymer frameshift (`+1N`)** at `NC_090942.1:141,844,154` (cod 95) + **3 non-homopolymer native PTCs**: `TAG` cod 72 (`141,844,221–223`), `TAA` cod 92 (`141,844,161–163`), `TGA` cod 106 (`141,844,120–122`); native CDS = `572 bp` | `1/2` canonical; **Intron 1 non-canonical `CC..CC` (849 bp)** | Chromosome-level (contig N50 `28.8 Mb`); `0` N-gaps | **CONFIRMED_ANNOTATION_COMPENSATED_PSEUDOGENE** |
| **Aptrow *AMELX*** | `GCF_003343035.1` (`XM_026079913.1` / `XP_025935698.1`, 205 aa, 3 exons) | Adjacent to `ARHGAP6` (`NW_020450766.1`); `XP_014463267.1` `121 b` ($E=4.0\times 10^{-34}$) | **0 frameshifts** (native CDS = `618 bp`); **1 native PTC (`TAA→X`)** at `NW_020450766.1:13,439–13,441` (cod 177, adjacent to `AAAA` run) | `1/2` canonical; **Intron 2 non-canonical `CA..CC` (848 bp)** | Scaffold-level (`aptRow1`, contig N50 `129.4 kb`); `0` N-gaps; **0 second conspecific assemblies** on NCBI | **ANNOTATION_COMPENSATED_LESION_SINGLE_ASSEMBLY** |

- **Manuscript Precision Fix:** Distinguish clearly that `3/4` chromosome-level loci (`Chrpic ENAM`, `Chemyd AMBN`, `Strcam AMELX`) carry multiple confirmed native frameshifts/PTCs/non-canonical introns compensated by Gnomon (`LOW QUALITY PROTEIN`), whereas `Aptrow AMELX` (`aptRow1` scaffold-level) has 0 annotated frameshifts and carries 1 Gnomon-compensated PTC (`TAA` at codon 177) plus a non-canonical `CA..CC` intron without an independent conspecific assembly to exclude rare sequencing error (though its `LOW QUALITY PROTEIN` status and non-canonical splice junction are native RefSeq features, not inter-species alignment artifacts).

### 2.3. Phase B3: Exact Coordinate & Indexing Reconciliation
- **Verified Indexing Convention (explicitly declared in Methods, Table 1, Table S2a, and Figure 1):**
  - Genomic intervals reported as `start–end` are **1-based inclusive** (span $= \text{end} - \text{start} + 1$); corresponding **0-based half-open BED** intervals are `[start - 1, end)`.
- **Four Reconciled Differences:**
  1. **`GRCg7b` *NRG1–LPL* Intergenic Gap:** RefSeq *NRG1* ends at `55,109,784` and *LPL* starts at `55,215,002` (1-based inclusive). The endpoint difference `55,215,002 - 55,109,784 = 105,218 bp` (reported in `Table_S2a` column `intergenic_gap_bp`), whereas the **strictly intergenic nucleotides** (`NC_052572.1:55,109,785–55,215,001`, BED `55109784–55215001`) span **`105,217 bp`** (both columns now reported explicitly in `Table_S2a`).
  2. **Figure 1B Plotted Span (`138 kb` vs. `105.2 kb`):** Figure 1B plots `NC_052572.1:55.098–55.236 Mb` (`138,000 bp`) to display `11.8 kb` of the 5' *NRG1* flank + the full `105,217-bp` intergenic gap + the full `17,437-bp` *LPL* locus (`55,215,002–55,232,438`) + `3.6 kb` right margin (`PSD3` edge). Caption explicitly clarifies both spans.
  3. **$\Psi$*ENAM* Exon 8 (`1,985 bp` Genomic Union vs. `1,975 bp` `BLASTn` Subject Span vs. `2,074` Alignment Columns):**
     - `tBLASTn` core (`Table_S2b`): `NC_052572.1:55,155,480–55,157,125 (−)` (`1,646 bp`, BED `55155479–55157125`).
     - Intra-avian `BLASTn` vs. kakapo (`ENAM_exon8_aves`): `NC_052572.1:55,155,490–55,157,464 (−)` (`1,975 bp` genomic span, BED `55155489–55157464`, `2,074` gapped alignment columns, `63.4%` identity, `317 bits`, $E = 7.32 \times 10^{-87}$).
     - **Union span:** `NC_052572.1:55,155,480–55,157,464 (−)` = **`1,985 bp`** (1-based inclusive, BED `55155479–55157464`), differing by `10 bp` at the 3' (`55,155,480` vs. `55,155,490`) edge. `v1.0.0` erroneously printed `1,975 bp` next to `55,155,480–55,157,464`; fully reconciled in `v1.0.1`!
  4. **$\Psi$*ENAM* Exon 6 (`149 bp` `BLASTn` Span vs. `47 bp` Exon Core) & Distances:**
     - `BLASTn` hit span: `NC_052572.1:55,160,529–55,160,677 (−)` (`149 bp`, BED `55160528–55160677`), `3,064 bp` strictly upstream of Exon 8 (`55,157,464`) and `1,108 bp` strictly downstream of Ensembl lncRNA `ENSGALG00010014557` (`55,161,786–55,165,166`).
     - Exact homologous Exon 6 coding core: `NC_052572.1:55,160,570–55,160,616 (−)` (`47 bp`, BED `55160569–55160616`), `3,105 bp` strictly upstream of Exon 8 and `1,169 bp` from `ENSGALG00010014557`.

---

## 3. Phase D — Mutation Calling, IQ-TREE 3 ML ASR, and Splice-Border Audit

- **Script & Outputs:** `scripts/enamel_relics/10_audit_mutations_and_asr.py` → `runs/enamel_relics/v101/phaseD/{Table_S5a_level4_shared_mutations_v101.tsv, Table_S5b_level4_splice_borders_v101.tsv, phaseD_audit_summary.json}`.
- **Species Identifier Disambiguation:** Replaced ambiguous 4-letter codes (`Tgut` vs. `Tgua`) with unambiguous 6-letter codes (`Taegut` = *Taeniopygia guttata*, `Tingut` = *Tinamus guttatus*; complete 19-species map in `Table_S0_species_identifiers.tsv`).
- **Standardized Variant Nomenclature:** Replaced incorrect HGVS `c.` prefix (which denotes cDNA nucleotide coordinates, not codons) with explicit protein-anchored notation (`XP_019346377.1:p.Pro311Ter` / codon 311, `XP_019346386.1:codon_109_del16nt_fs`).
- **Topology & ML ASR (`IQ-TREE 3.1.4`, `GTR+G --ancestral` on avian *ENAM* Exon 8):**
  - *ENAM* Exon 8 itself strongly prefers the phylogenomic consensus topology **T1** (Novaeratitae = *Apteryx* + *Dromaius*, *Struthio* basal; $\ln L = -15,237.86$) over **T2** (*Struthio*+*Dromaius*, $\ln L = -15,260.52$, $\Delta\ln L = -22.66$) and **T3** (*Struthio*+*Apteryx*, $\ln L = -15,260.72$, $\Delta\ln L = -22.86$).
  - Under **T1**, the 23 dual-aligner consensus shared mutations (16 avian *ENAM* Exon 8 PTCs + 1 Galloanserae frameshift + 2 turtle *ENAM* PTCs + 2 avian *AMBN* Exon 3 PTCs + 1 Palaeognathae *AMBN* Exon 4 16-nt deletion + 1 Durocryptodira *AMBN* Exon 4 1-nt insertion) reproduce **23/23** exactly across both MAFFT L-INS-i and MUSCLE 5.
  - **Key Mechanistic & Evolutionary Corrections to `v1.0.0`:**
    1. **ENAM codon 311 (`p.Pro311Ter`):** `v1.0.0` labeled `CCG→TGA` (Alligator reference, 3 nt steps) as "CpG deamination." IQ-TREE 3 ML ASR resolves the actual trajectory: Neornithes `CAA` ($P=0.542$) → Palaeognathae **`CGA`** ($P=0.541$) vs. **`TGA`** ($P=0.453$). Thus, within Palaeognathae, `CGA→TGA` *is* a genuine 1-step CpG transition, and ML ASR explicitly shows ambiguity between parallel `CGA→TGA` transitions in `Strcam`/`Dronov` versus a single Palaeognathae `CGA→TGA` stop followed by `TGA→CGA` reversal in `Aptrow`.
    2. **ENAM codon 1022 (`p.Ser1022Ter`, `Galgal` + `Strhab` = `TAG`, `Anapla` missing, `Taegut` = `AAA`):** Under irreversible parsimony (`v1.0.0`), codon 1022 counts among the 6 incongruent PTCs (`6/16 = 37.5%`, exact Clopper–Pearson 95% CI `15.2–64.6%`). Under reversible ML ASR (`GTR+G`), however, codon 1022 is reconstructed as a **putative crown Neognathae synapomorphy** (`TAG`, $P_{\text{stop}}(\text{Neognathae}) = 0.926$) with secondary decay to `AAA` on the long passerine branch (`Taegut`), reducing the strictly convergent subset to **`5/16` (`31.2%`, 95% CI `11.0–58.7%`)** (or `4–6/16` accounting for codon 311 reversal ambiguity). Under alternative topologies T2/T3, incongruent PTCs increase to `10/16` (`62.5%`) and `11/16` (`68.8%`).
    3. **Distinct Stop Triplets (`p.Ala439Ter`, `p.Trp615Ter`, `p.Gln684Ter`):** `v1.0.0` overstated that distinct stop triplets *prove* independent origins by themselves. Re-audit shows `TAA↔TAG` (codon 684) and `TAA↔TGA` (codon 615) differ by only 1 transition (allowing stop-to-stop drift after a shared stop in principle), whereas `TAG↔TGA` (codon 439) differs by 2 substitutions. Independent origin at codons 439, 615, and 684 is supported jointly by **phylogenetic distribution + ML ASR** ($P_{\text{stop}}(\text{Neornithes}) = 0.221, 0.304, 0.366$; $P(\text{Palaeognathae TAG at 684}) = 0.974$ vs. $P(\text{Galloanserae CAA}) = 0.714$), not triplet identity alone.
    4. **Neutral Stop-Gain Monte Carlo Expectation (1,000 replicates, K80 $\kappa=4$, seed `20261122`):** Simulating neutral evolution of the 2,634-bp *ENAM* Exon 8 along the empirical IQ-TREE 3 7-bird branch lengths yields a median of **26–29 shared PTCs**, **48.1–51.8% incongruent** (95% simulation range `29.2–70.8%`), and a **median of 0 crown-Neornithes PTC synapomorphies** even when inactivation precedes the crown split (`t_stem = 0.02–0.05`). The observed `31.2–37.5%` incongruence sits squarely within neutral expectations.
    5. **Splice-Border Table `Table_S5b` Artifact Fixed:** `v1.0.0` sliced 2 bp outside raw local `BLASTn` endpoints for *AMBN* Exons 3, 4, and 9 (generating spurious non-canonical borders and calling a nonexistent 5' donor on terminal Exon 9). Re-auditing with ±80-bp intronic flanks aligned in both MAFFT and MUSCLE confirms:
       - *ENAM* Exon 6 (100% MAFFT/MUSCLE agreement): `Allmis/Galgal/Anapla = AG..GT`; `Strcam = AC..GT`; `Aptrow = AC..AT`; `Dronov = GT..AT`.
       - *AMBN* Exon 3 (100% MAFFT/MUSCLE agreement in Aves): `Anapla/Aptrow/Dronov/Strcam = AG..GT` (canonical!), while **`Galgal = AC..GT`** (chicken-specific 3' acceptor disruption!).
       - *AMBN* Exon 4 (100% agreement in Durocryptodira): `Chemyd/Chrpic = AC..GC` (shared Durocryptodira 3' acceptor disruption + minor canonical `GC` donor).
       - *AMBN* Exon 9 (terminal coding exon, no 5' splice donor) removed from `Table_S5b`.

---

## 4. Phase C — High-Permutation Statistical Re-evaluation (1,000 Permutations + 20-aa Block Null)

- **Script & Outputs:** `scripts/enamel_relics/09_permutation_sensitivity.py` → `runs/enamel_relics/v101/stats/{perm_sensitivity_per_gene_v101.tsv, retention_indices_1000perm_v101.tsv, clade_medians_equidistant_1000perm_v101.tsv, phaseC_summary.json}` (mirrored into `Table_S3b` and `Table_S3c`).
- **Scope & Execution:** `72` parallel `tblastn v2.17.0+` runs (`18` species × `2` domains × `2` SEG modes, exact `-subject` per species, `-evalue 10 -max_target_seqs 50`) evaluating `1,000` global composition-preserving shuffles per query (seed `20261119`, nested prefixes `N = 20, 200, 500, 1,000`) plus `200` local 20-aa block-preserving shuffles (seed `20261120`) for all equidistant queries and `SPARCL1_Amis`, with `2,000`-replicate bootstrap 95% CIs (seed `20261121`). Runtime: `1,184.7 s`.
- **Exact Reproduction & Precision Gain (**PASSOU**):**
  - Observed bitscores reproduced **`504/504`** identically against `v1.0.0`.
  - Median $|p_{95}(N=1000) - p_{95}(\text{orig } N=20)| = 0.86\text{ bits}$ (`95th percentile = 2.13 bits`; `max = 25.35 bits` solely in heavily `SEG`-masked `DSPP_Amis`, where `N=20` collapsed to `0.0` with a `25.87-bit` bootstrap CI width vs. `0.70 bits` at `N=1,000`).
  - Across all `504` cells, median bootstrap 95% CI width on $p_{95}$ narrowed **3.4-fold** from `2.70 bits` (`N=20`) to `0.80 bits` (`N=1,000`).
  - `12` marginal noise-floor flips occurred (`0.26–1.42 bits` excess), cleanly dropping `Pelsin` (`Psin`) and `Tingut` (`Tgua`) *Alligator* `ENAM` hits to `0.00` excess (empirical $p = 0.061$ and $0.052$) while elevating `Strhab` `ENAM_Acar` to `0.80 bits` ($p = 0.024$).
- **Effect of Query Length, Composition, `SEG`, and Clade Substitution Rates:**
  - Null $p_{95}$ (`N=1,000`, $n=14$ queries) correlates positively with query length (`-seg yes` Spearman $\rho = +0.651$, $p = 0.0117$; `-seg no` $\rho = +0.702$, $p = 0.0051$) and Shannon entropy ($\rho = +0.769$, $p = 0.0013$) and **negatively** with `SEG`-masked fraction (`-seg yes` Pearson $r = -0.795$, $p = 0.0007$; `-seg no` $r = -0.670$, $p = 0.0087$, reflecting BLAST+ composition-based score adjustment even without `SEG`).
  - Equidistant divergence (`Anolis`, ~281 Ma) does **not** imply equal substitution rates: median top-HSP *LPL* amino-acid identity is `80.0%` (Crocodylia), `73.9%` (Testudines), and `88.7%` (Aves), proving clade- and gene-specific rate heterogeneity rather than a molecular clock.
- **Core Benchmark Findings at `N = 1,000` (and 20-aa Block Null):**
  - Globular `LPL_Acar` (`0.0%` SEG): Aves median `135.6` (`133.1–143.5`) vs. Crocodylia `120.4` (`119.4–121.4`) → **`112.6%`** (`-seg yes` & `-seg no` identical; 20-aa block null `135.9`).
  - Intact bone `IBSP_Acar` (`52.8%` SEG): Aves **`74.5%`** (`40.5` vs. `54.3`, `-seg yes`; block null `42.7`) and **`48.5%`** (`48.2` vs. `99.3`, `-seg no`; `51.5%` attenuation).
  - Intact bone `SPP1_Acar` (`25.5%` SEG): Aves **`50.0%`** (`16.8`, range `7.3–23.9`, vs. Crocodylia `33.5`, `-seg yes`; block null `18.3`) and **`58.0%`** (`21.0` vs. `36.2`, `-seg no`), directly overlapping decayed avian `ENAM_Acar` (`0.0–12.7 bits`) and `AMBN_Acar` (`0.0–15.4 bits`).
  - *LPL*-normalized retention indices (equidistant `Anolis`, `-seg yes`, `N=1,000`, bootstrap 95% CIs): Crocodylia `2.27 [2.26–2.28]` to `2.80 [2.80–2.82]`; Testudines `0.00` to `0.58 [0.57–0.58]`; Aves `0.00` to `0.12 [0.11–0.12]` (all denominators `LPL_exc` `119.4–256.1 bits`, highly stable).
  - Because birds and turtles represent only **$n = 2$ independent origins of edentulism**, all clade comparisons remain strictly descriptive (zero pooled-species inferential $p$-values claimed).

