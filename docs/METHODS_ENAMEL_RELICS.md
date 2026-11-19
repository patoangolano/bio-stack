# Methods & Provenance Record — Enamel-Gene Relics in Birds and Turtles

*Auto-compiled from verified provenance logs (`runs/enamel_relics/provenance_run.jsonl`, `runs/enamel_relics/tool_versions.txt`, `runs/enamel_relics/audit_metadata.json`, and `runs/enamel_relics/ensembl_116_ggal_audit.json`). Separate from `docs/metodos.md` (which belongs to another repository project).*

## 1. Taxon Sampling & RefSeq / Ensembl Annotation Audit

Eighteen sauropsid RefSeq reference assemblies spanning all four major extant lineages — Crocodylia (*n* = 3: *Alligator mississippiensis* `GCF_030867095.1`, *Gavialis gangeticus* `GCF_001723915.1`, *Crocodylus porosus* `GCF_001723895.1`), Squamata (*n* = 3: *Anolis carolinensis* `GCF_035594765.1`, *Python bivittatus* `GCF_059951505.1`, *Podarcis muralis* `GCF_964188315.1`), Testudines (*n* = 4: *Chrysemys picta* `GCF_011386835.1`, *Chelonia mydas* `GCF_015237465.2`, *Gopherus evgoodei* `GCF_007399415.2`, *Pelodiscus sinensis* `GCF_052426305.1`), and Aves (*n* = 8: *Gallus gallus* `GCF_016699485.2`, *Anas platyrhynchos* `GCF_047663525.1`, *Strigops habroptila* `GCF_004027225.2`, *Taeniopygia guttata* `GCF_048771995.1`, *Dromaius novaehollandiae* `GCF_036370855.1`, *Struthio camelus* `GCF_040807025.1`, *Apteryx rowi* `GCF_003343035.1`, *Tinamus guttatus* `GCF_059628265.1`) — plus *Homo sapiens* `GRCh38.p14` (`GCF_000001405.40`) as mammalian outgroup reference were audited via NCBI Datasets API v2 (`2026-11-19`, full details in `runs/enamel_relics/supplementary/Table_S1_assemblies_and_annotations.tsv`).

For *Gallus gallus* `GRCg7b` (`GCF_016699485.2`), both NCBI RefSeq Annotation Release 106 (`Ggal_genomic.gff`) and Ensembl Release 116 (`Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz`, SHA-256 `407e3e270da9668eecc205b685493d157abbd8fa80bc698a65bfd3ed3a5f8009`) were audited locally to verify genome-wide absence of enamel-gene annotations and zero coordinate overlap with recovered chrZ relics.

## 2. Syntenic Interval Extraction & Controls

Two genomic intervals were extracted per species across all 18 sauropsids:
1. **Enamel-cluster locus (*NRG1–LPL*):** anchored on *NRG1* and *LPL* (`±60 kb` padding; single-anchor `±800 kb` window for *Dromaius* and *Tinamus* where one flanking symbol is unannotated in RefSeq; `Table_S2a`).
2. **Osteogenic IDP control locus (*IBSP–SPP1–SPARCL1*):** anchored on conserved SIBLING/osteogenic genes (`±100 kb` padding) to evaluate functional intrinsically disordered proteins across identical phylogenetic distances.

## 3. Compositional Profiling, SEG Sensitivity, and Null-Corrected Retention

All 16 query proteins (full RefSeq *Alligator* `ENAM` `XP_019346377.1` [1,078 aa], legacy partial `QOX06214.1` [849 aa], `AMBN` `XP_019346386.1`, `AMTN` `XP_014451634.1` and `KYO35092.1`, globular positive control `LPL` `XP_019346320.1`, cross-chromosome negative control `DSPP` `XP_059578269.1`; equidistant *Anolis* `ENAM` `XP_062825318.1`, `AMBN` `NP_001315156.1`, `AMTN` `NP_001315157.1`, `LPL` `XP_062827593.1`, `IBSP` `XP_003221343.2`, `SPP1` `XP_003221344.2`; *Alligator* `IBSP`, `SPP1`, `SPARCL1`) were profiled with `segmasker` `1.0.0` (`-window 12 -locut 2.2 -hicut 2.5`) and Shannon entropy. Twenty composition-preserving random permutations per query were generated deterministically using `numpy.random.default_rng(20261017)` (`01_audit_and_fetch.py`). Translated searches were executed with NCBI `tblastn` `2.17.0+` under both `-seg yes` and `-seg no` (`-evalue 10 -max_target_seqs 100`).

## 4. Whole-Proteome Reciprocal BLASTx & Single-Copy Verification

Every candidate locus (`bitscore >= 28` or `E <= 1e-2`) was tested by `blastx` `2.17.0+` (`-evalue 10 -seg yes -max_target_seqs 25`) against the entire *Alligator mississippiensis* RefSeq proteome (`48,040` proteins, `GCF_030867095.1`) both at the clustered locus-span level and at the top-HSP exon level (`±30 bp`) to prevent intronic hits from masking short orthologous exons (`Table_S2b`). Chicken chrZ relics were further tested for single-copy specificity via `blastn -task dc-megablast` (`-evalue 1e-3 -dust yes`).

## 5. Level 3 Splice-Aware Coding Lesions & Level 4 Dual-Aligner Shared Mutations

Splice-aware protein-to-genome alignments were performed with `miniprot` `0.18-r281` (`-j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5`) across all 18 intervals for `ENAM`, `AMBN`, `AMTN`, and `LPL`. Alignments containing introns exceeding `25,000 bp` were segmented to the primary orthologous block to filter distal low-complexity chaining artifacts (specifically *Chrysemys picta* `ENAM`'s 58.4-kb jump to `protein SON`, `e = 0.010`; `Table_S4`).

For Level 4 shared inactivating mutations, verified coding exons of *Alligator* `ENAM` (Exons 6 and 8) and `AMBN` (Exons 3, 4, and 9) were aligned with orthologous crocodylian controls (*Alligator*, *Crocodylus*, *Gavialis*; required `0` lesions in all alignments) and avian/testudine relics independently with **MAFFT `v7.525` (`--localpair --maxiterate 1000`)** and **MUSCLE `5.3` (`-align`)**. Only mutations identically called in both aligners (exact reference coordinate, mutation class, and carrier set) were retained (`Table_S5a`, `Table_S5b`).

## 6. Provenance & Verification Manifest

All SHA-256 hashes of final supplementary TSVs are recorded in `runs/enamel_relics/supplementary/MANIFEST_SHA256.tsv`. Execution provenance is logged in `runs/enamel_relics/provenance_run.jsonl` and `runs/provenance.jsonl`.
