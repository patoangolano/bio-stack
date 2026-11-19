# bioRxiv Preprint Submission Package

**Manuscript Title:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles*
**Author:** Matheus Lima (Independent Researcher, Brazil; `limatheus458@gmail.com`; ORCID: `0009-0000-3557-8337`)
**Target bioRxiv Subject Area:** **Evolutionary Biology** (secondary: **Genomics**)
**Article Type:** New Results (Empirical Comparative Genomics / Paleogenomics)

## Package Contents

| File / Directory | Purpose |
|---|---|
| `manuscript_biorxiv.pdf` | Ready-to-upload single PDF (16 pages, continuous line numbering, page numbering, 320 DPI embedded figures, 46/46 verified APA references) |
| `manuscript.typ` | Reproducible Typst `0.15.1` source (`typst compile manuscript.typ manuscript_biorxiv.pdf`) |
| `manuscript.md` | Clean GitHub-Flavored Markdown mirror synchronized with `manuscript.typ` |
| `references.bib` | 46 primary references verified against PubMed, Europe PMC, Crossref, and Retraction Watch (0 retractions, 0 expressions of concern) |
| `figures/` | Figures 1–3 exported in vector `.svg` (path-outlined), vector `.pdf` (Type 42 fonts), and 320 DPI `.png` (Okabe-Ito colorblind-safe palette) |
| `supplementary/` | Consolidated supplementary tables `Table_S1` through `Table_S6` (10 `.tsv` files) + `MANIFEST_SHA256.tsv` |

## Pre-Upload Checklist (All Automated Metadata Complete)

1. **Author ORCID:** Inserted (`0009-0000-3557-8337`, clickable link) and compiled into `manuscript_biorxiv.pdf`.
2. **Zenodo Archival DOI:** Minted (`10.5281/zenodo.24358067`, release `enamel-relics-v1.0.0`), inserted into `Data and Code Availability`, and compiled into `manuscript_biorxiv.pdf`.
3. **Upload to bioRxiv:**
   - Category: **Evolutionary Biology** (Article Type: **New Results**)
   - Main manuscript file: `manuscript_biorxiv.pdf` (16 pages, continuous line & page numbering, zero placeholders remaining)
   - Supplementary data files: zip `supplementary/` (`Table_S1`–`Table_S6` + `MANIFEST_SHA256.tsv`) or upload individual TSVs as Supplementary Material.
   - Abstract field: copy directly from `manuscript.md` (verified **236 words**, ≤ 250-word limit).
