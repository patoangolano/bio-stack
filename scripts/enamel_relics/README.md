# Enamel-Gene Relics in Birds and Turtles — Reproducible Analysis & Audit Pipeline (v1.0.1 Candidate)

Deterministic comparative-genomics and paleogenomics pipeline underlying the manuscript:
**"Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles"**.

## Environment & Pinned Tool Versions

All alignments, phylogenetic reconstructions, and PDF compilations were executed under WSL2 (`Ubuntu-26.04`, Linux `6.18.58.1-microsoft-standard-WSL2`) using Pixi-managed binaries (`~/.pixi/bin`) orchestrated from Python 3.11 (`uv` virtual environment on Windows 11). Exact runtime versions are captured in `runs/enamel_relics/tool_versions.txt` and re-verified in `runs/enamel_relics/v101/phaseAB/phase_A_B1_B3_audit.json`:

| Tool | Version | Role |
|---|---|---|
| NCBI BLAST+ (`tblastn`, `blastx`, `blastn`, `blastp`) | `2.17.0+` | Syntenic homology search (`-seg yes`/`-seg no`, `-max_target_seqs 50`), whole-proteome reciprocal `BLASTx`/`BLASTp`, single-copy `dc-megablast` |
| `segmasker` | `1.0.0` | Low-complexity profiling (`-window 12 -locut 2.2 -hicut 2.5`) |
| `miniprot` | `0.18-r281` | Splice-aware protein-to-genome alignment (Level 3 frameshift, PTC, and `GT-AG` splice-site audit) |
| `MAFFT` | `v7.525` | `L-INS-i` multiple sequence alignment (`--localpair --maxiterate 1000`) |
| `MUSCLE` | `5.3` | Independent second aligner for dual-aligner consensus filtering |
| `IQ-TREE 3` | `3.1.4` | Maximum-likelihood ancestral state reconstruction (`-m GTR+G --ancestral`) & Palaeognathae topology comparison |
| `samtools` / `seqkit` | `1.24` / `v2.14.0` | Indexed sequence extraction and FASTA inspection |
| `NCBI Datasets CLI` | `18.38.0` | Assembly metadata and annotation auditing |
| `Typst` | `0.15.1` | Manuscript PDF compilation with continuous line and page numbering |
| Python | `3.11` (`numpy`, `scipy`, `matplotlib`, `pytest`) | Orchestration, permutation/bootstrap statistics, ASR post-processing, figures, and acceptance tests |

### Deterministic Random Seeds (all via `numpy.random.default_rng()`, never global `np.random.seed()`)
- `20261017`: Baseline 20 composition-matched query permutations (`01_audit_and_fetch.py`) and jitter layout (`06_make_figures.py`).
- `20261119`: High-permutation sensitivity suite (`N = 20, 200, 500, 1000` global composition-preserving shuffles per query; `09_permutation_sensitivity.py`).
- `20261120`: Local 20-aa block-preserving shuffles (`N = 200` per equidistant query; `09_permutation_sensitivity.py`).
- `20261121`: Bootstrap 95% confidence intervals (`B = 2,000`) on null `p95`, excess bitscores, and *LPL*-normalized retention indices (`09_permutation_sensitivity.py`).
- `20261122`: Neutral stop-gain Monte Carlo simulation (`1,000` replicates, K80 $\kappa = 4$) along empirical IQ-TREE 3 *ENAM* Exon 8 branch lengths (`10_audit_mutations_and_asr.py`).

## End-to-End Execution Order

Run from the repository root (`C:\Users\user\bio-stack`):

```powershell
# --- Core Data Acquisition & Primary Alignments (Requires NCBI/Ensembl Network Access) ---
# Note: Raw refdata/blastdb/intervals (~1.16 GB) are gitignored; tracked tables & level4 MSAs allow
# offline regeneration of supplementary tables, figures, and acceptance tests.
uv run python scripts/enamel_relics/01_audit_and_fetch.py
uv run python scripts/enamel_relics/01b_audit_ensembl_and_psin.py
wsl -d Ubuntu-26.04 -- bash -l /mnt/c/Users/user/bio-stack/scripts/enamel_relics/02_run_alignments_wsl.sh
uv run python scripts/enamel_relics/03_reciprocal_and_compositional.py
uv run python scripts/enamel_relics/04_level4_shared_mutations.py

# --- v1.0.1 Scientific & Statistical Audit Suite ---
# Phase A/B1/B3: SHA-256 verification, Ggal Psi-ENAM Exon 6 reciprocal BLASTx audit, coordinate reconciliation
uv run python scripts/enamel_relics/07_audit_phaseA_B1_B3.py

# Phase B2: Forensic native-model audit of the 4 RefSeq PROTEIN_CODING models (Table S7)
uv run python scripts/enamel_relics/08_audit_protein_coding_models.py
uv run python scripts/enamel_relics/08b_b2_native_exons_and_replication.py
uv run python scripts/enamel_relics/08c_b2_exact_ptc_and_tableS7.py

# Phase C: 1,000-permutation sensitivity + 20-aa block null + bootstrap CIs + length/composition regression
uv run python scripts/enamel_relics/09_permutation_sensitivity.py

# Phase D: Mutation re-audit, IQ-TREE 3 ML ASR (T1/T2/T3), flanked-MSA splice borders, neutral simulation
uv run python scripts/enamel_relics/10_audit_mutations_and_asr.py

# --- Consolidation, Figures, PDF Compilation & Acceptance Tests ---
uv run python scripts/enamel_relics/11_consolidate_v101_tables.py
uv run python scripts/enamel_relics/06_make_figures.py
wsl -d Ubuntu-26.04 -- bash -lc "cd /mnt/c/Users/user/bio-stack/submission_biorxiv && typst compile manuscript.typ manuscript_biorxiv.pdf"
uv run pytest -v tests/test_enamel_relics_audit.py
```

## Outputs & Provenance

- **Provenance logs:** `runs/enamel_relics/provenance_run.jsonl` (v1.0.0 baseline) and `runs/enamel_relics/v101/provenance_*.jsonl` (v1.0.1 audit stages).
- **Audit reports:** `audit/{reproducibility_audit.md, bibliography_report.md, claims_evidence_matrix.tsv, change_log.md, biorxiv_readiness_report.md}`.
- **Supplementary tables:** `submission_biorxiv/supplementary/Table_S0` through `Table_S7` + `MANIFEST_SHA256.tsv` (mirrored from `runs/enamel_relics/supplementary/`).
- **Figures:** `submission_biorxiv/figures/fig1_synteny_map.{svg,pdf,png}`, `fig2_compositional_bias.{svg,pdf,png}`, `fig3_lesions_and_phylogeny.{svg,pdf,png}`.
- **Archival DOI:** Concept DOI `https://doi.org/10.5281/zenodo.24358066` (initial `v1.0.0` snapshot: `https://doi.org/10.5281/zenodo.24358067`).
