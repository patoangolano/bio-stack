# bioRxiv Preprint Submission Package (`enamel-relics-v1.0.2` Candidate)

**Manuscript Title:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles*  
**Author:** Matheus Lima (Independent Researcher, São Paulo, Brazil; `limatheus458@gmail.com`; ORCID: `0009-0000-3557-8337`)  
**Target bioRxiv Subject Area:** **Evolutionary Biology** (secondary: **Genomics**)  
**Article Type:** New Results (Empirical Comparative Genomics / Paleogenomics)  
**Archival Identifiers:** Zenodo Concept DOI: [`10.5281/zenodo.24358066`](https://doi.org/10.5281/zenodo.24358066) (always resolves to the latest release; prior audited `enamel-relics-v1.0.1` release preserved at [`10.5281/zenodo.24360150`](https://doi.org/10.5281/zenodo.24360150); initial `enamel-relics-v1.0.0` snapshot preserved untouched at [`10.5281/zenodo.24358067`](https://doi.org/10.5281/zenodo.24358067))

## Package Contents

| File / Directory | Purpose |
|---|---|
| `manuscript_biorxiv.pdf` | Recompiled single PDF (21 pages, continuous line numbering, page numbering, 320 DPI embedded figures, 49/49 verified APA references) |
| `manuscript.typ` | Reproducible Typst `0.15.1` source (`typst compile manuscript.typ manuscript_biorxiv.pdf`) |
| `manuscript.md` | Pure GitHub-Flavored Markdown mirror synchronized verbatim with `manuscript.typ` (including Tables 1–3) |
| `references.bib` | 49 primary references verified against Crossref (`49/49 HTTP 200`), Europe PMC (`47/49`), and Retraction Watch (`0` retractions, `0` expressions of concern; see `../audit/bibliography_report.md`) |
| `figures/` | Figures 1–3 derived strictly from `supplementary/Table_S*.tsv`, exported in vector `.svg`, vector `.pdf` (Type 42 fonts), and 320 DPI `.png` (Okabe-Ito colorblind-safe palette) |
| `supplementary/` | Consolidated supplementary tables `Table_S0` through `Table_S9` (14 `.tsv` files, strict LF line endings) + `MANIFEST_SHA256.tsv` (mirrored identically to `../runs/enamel_relics/supplementary/`) |

## Pre-Submission Human Action Checklist

See [`../audit/v102/human_approval_checklist.md`](../audit/v102/human_approval_checklist.md) and [`../audit/v102/biorxiv_readiness_report.md`](../audit/v102/biorxiv_readiness_report.md) for the complete pre-submission readiness audit:

1. **Affiliation & openRxiv Screening Compliance (`PENDENTE DE DECISÃO HUMANA`):**
   - Current draft conservatively uses **Option A** (`Independent Researcher, São Paulo, Brazil` + factual Senac student *Author Note*).
   - Under openRxiv's 20 October 2026 AI & Accountability policy, independent authors without an institutional affiliation face heightened screening scrutiny. Before uploading to bioRxiv, send Draft Emails 1 (to Centro Universitário Senac) and 2 (to `bioRxiv@openrxiv.org`) provided in `../audit/v102/biorxiv_readiness_report.md` and update the affiliation line to **Option B** only if written confirmation is received.
2. **Release `enamel-relics-v1.0.2` & Zenodo DOI:**
   - After human approval and merge of branch `audit/biorxiv-final-v1.0.2` into `main`, publish GitHub Release `enamel-relics-v1.0.2` so Zenodo automatically mints the new `v1.0.2` Version DOI under Concept DOI [`10.5281/zenodo.24358066`](https://doi.org/10.5281/zenodo.24358066).
3. **Upload to bioRxiv:**
   - Category: **Evolutionary Biology** (secondary: **Genomics**); Article Type: **New Results**.
   - Main manuscript file: `manuscript_biorxiv.pdf` (21 pages, continuous line & page numbering, zero placeholders).
   - Supplementary package: compress `supplementary/` (`Table_S0`–`Table_S9` + `MANIFEST_SHA256.tsv`, 15 files total) into a zip archive or upload individual TSVs as Supplementary Material.
   - Abstract field: copy verbatim from `manuscript.md` (verified **247 words**, ≤ 250-word bioRxiv limit).
   - Author ORCID: `0009-0000-3557-8337`.
