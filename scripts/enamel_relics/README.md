# Enamel-Gene Relics in Birds and Turtles — Reproducible Analysis Pipeline

Deterministic comparative-genomics pipeline underlying the manuscript:
**"Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles"**.

## Environment & Pinned Tool Versions

All alignments and sequence manipulations were executed under WSL2 (`Ubuntu-26.04`, Linux `6.18.58.1-microsoft-standard-WSL2`) using Pixi-managed bioinformatics binaries (`~/.pixi/bin`) orchestrated from Python 3.11 (`uv` virtual environment on Windows 11). Exact versions captured at runtime in `runs/enamel_relics/tool_versions.txt`:

| Tool | Version | Role |
|---|---|---|
| NCBI BLAST+ (`tblastn`, `blastx`, `blastn`) | `2.17.0+` | Syntenic homology search (`-seg yes`/`-seg no`), whole-proteome reciprocal BLASTx, single-copy `dc-megablast` |
| `segmasker` | `1.0.0` | Low-complexity profiling (`-window 12 -locut 2.2 -hicut 2.5`) |
| `miniprot` | `0.18-r281` | Splice-aware protein-to-genome alignment (Level 3 frameshift, PTC, and `GT-AG` splice-site audit) |
| `MAFFT` | `v7.525` | `L-INS-i` multiple sequence alignment for Level 4 exon-level mutation calling |
| `MUSCLE` | `5.3` | Independent second aligner for dual-aligner consensus filtering |
| `samtools` / `seqkit` | `1.24` / `v2.14.0` | Indexed sequence extraction and FASTA inspection |
| `NCBI Datasets CLI` | `18.38.0` | Assembly metadata and annotation auditing |
| `Typst` | `0.15.1` | Manuscript PDF compilation with continuous line and page numbering |
| Python | `3.11` (`numpy`, `matplotlib`, `scipy`) | Orchestration, compositional/null statistics, tables, and vector figures |

Deterministic random seed for all 20 composition-matched query permutations: `20261017` via `numpy.random.default_rng(20261017)` (never global `np.random.seed()`).

## End-to-End Execution Order

Run from the repository root (`C:\Users\user\bio-stack`):

```powershell
# Step 1: Audit 18 study assemblies + GRCh38.p14, download NRG1-LPL enamel intervals,
#         SIBLING/osteogenic IDP intervals, Alligator proteome (48,040 proteins), GFF3s,
#         and generate 20 composition-matched shuffles per query (seed=20261017)
uv run python scripts/enamel_relics/01_audit_and_fetch.py

# Step 1b: Audit Ensembl Release 116 chicken GRCg7b GTF and sync Pelodiscus sinensis T2T (GCF_052426305.1)
uv run python scripts/enamel_relics/01b_audit_ensembl_and_psin.py

# Step 2: Run segmasker, tBLASTn (-seg yes and -seg no) across all 18 enamel & 18 IDP intervals,
#         miniprot across all 18 intervals, and build the 48,040-protein Alligator BLAST DB in WSL2
wsl -d Ubuntu-26.04 -- bash -l /mnt/c/Users/user/bio-stack/scripts/enamel_relics/02_run_alignments_wsl.sh

# Step 3: Compute compositional/SEG metrics, null-corrected retention indices,
#         whole-proteome reciprocal BLASTx (span + top-HSP), and Level 3 miniprot + splice-site audit
uv run python scripts/enamel_relics/03_reciprocal_and_compositional.py

# Step 4: Dual-aligner (MAFFT L-INS-i + MUSCLE 5) Level 4 shared inactivating mutations
#         and exon-border splice-site analysis across avian & testudine ENAM and AMBN
uv run python scripts/enamel_relics/04_level4_shared_mutations.py

# Step 5: Consolidate supplementary tables (Table S1 - Table S6) and SHA-256 manifest
uv run python scripts/enamel_relics/05_consolidate_tables.py

# Step 6: Generate publication-grade vector (SVG/PDF) and 300 DPI PNG figures (Figures 1-3)
uv run python scripts/enamel_relics/06_make_figures.py
```

## Outputs & Provenance

- **Provenance log:** `runs/enamel_relics/provenance_run.jsonl` (mirrored to `runs/provenance.jsonl`) records timestamps, script paths, parameters, and SHA-256 checksums of all generated tables.
- **Supplementary tables:** `runs/enamel_relics/supplementary/Table_S1` through `Table_S6` + `MANIFEST_SHA256.tsv`.
- **Figures:** `runs/enamel_relics/figures/fig1_synteny_map.{svg,pdf,png}`, `fig2_compositional_bias.{svg,pdf,png}`, `fig3_lesions_and_phylogeny.{svg,pdf,png}`.
- **Submission bundle:** `submission_biorxiv/` contains the compiled PDF (`manuscript_biorxiv.pdf`), Typst/Markdown sources, verified BibTeX (`references.bib`), figures, and supplementary tables.
