# Pre-Submission Readiness Decision Report (`submission_readiness.md`)

**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final`  
**Manuscript:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles* (`21` pages, `submission_biorxiv/manuscript_biorxiv.pdf`)

---

## 1. Three-Dimensional Readiness Assessment

In accordance with strict scientific audit guidelines, we do **not** assign arbitrary numerical percentages of bioRxiv acceptance. Instead, readiness is evaluated separately across three distinct dimensions:

| Dimension | Overall Classification | Basis of Assessment |
|---|---|---|
| **1. Prontidão Documental (Documentary & Editorial Readiness)** | **`PASSOU` (100% Complete)** | - Single self-contained 21-page PDF (`manuscript_biorxiv.pdf`, lines `1–679`) compiled from `Typst 0.15.1` with searchable vector text, embedded 320-DPI figures (Okabe-Ito palette), and zero placeholders.<br>- Abstract verified at **248 words** ($\le 250$-word bioRxiv limit) and synchronized verbatim between `manuscript.typ` and `manuscript.md`.<br>- `49/49` primary references verified (`HTTP 200` on Crossref, `0` retractions, `1:1` citation parity).<br>- `14` supplementary TSV tables (`Table_S0`–`Table_S9`) + `MANIFEST_SHA256.tsv` verified with strict `LF` line endings and exact SHA-256 hashes. |
| **2. Validação e Consistência Científica Interna (Internal Scientific & Statistical Consistency)** | **`PASSOU / CORRIGIDO` (`0` Known Unresolved Inconsistencies)** | - **Phase A (`CORRIGIDO`):** Non-detection of crown-Neornithes synapomorphies is explicitly bounded by `0/878` 7-bird codon coverage and `16.3%` neutral uncontradicted ancestral-stop detection power (`40.8%` secondary stop-to-sense reversal corruption).<br>- **Phase B (`CORRIGIDO`):** Coverage-masked K80 simulations (`Table_S8`) are framed as statistical compatibility conditional on the K80 model, branch lengths, topology, and empirical coverage mask.<br>- **Phase C (`CORRIGIDO`):** Marginal empirical $p = 0.024$ (`Strhab` `ENAM_Acar`) is labeled unadjusted exploratory across the 504-cell grid; Clopper–Pearson CIs are qualified for single-exon site linkage; `T1` preference is anchored on external phylogenomics rather than single-exon RELL tests.<br>- **Phase D (`CORRIGIDO` / `PASSOU`):** Scope of the flanking globular *LPL* control is calibrated; graded orthology tiers (`Table_S2b`), Gnomon native-CDS audit (`Table_S7`), and reference-vs.-ancestral codon translations (`Table_S5a`) are preserved.<br>- **Out-of-Scope External Validations (`NÃO VERIFICADO`):** Second conspecific *Apteryx rowi* assembly, chromosome-level Pleurodira assemblies, branch $d_N/d_S$ modeling, and embryonic oral RNA-seq are explicitly declared as Limitations in Section 3.5. |
| **3. Prontidão Administrativa (Administrative & Portal Submission Readiness)** | **`PENDENTE` (Awaiting Human Author Portal Action)** | - The human author has already registered a bioRxiv account using their institutional Centro Universitário Senac email address, accessed the Author Area, and opened the submission form.<br>- Before clicking **Submit**, the author must make the final human choice between **Option A** (keeping `"Independent Researcher, São Paulo, Brazil"` + factual Senac student *Author Note*, currently compiled in the PDF, and pasting the explanatory note to screeners in the portal comments box) and **Option B** (`"Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil"` + explicit non-endorsement disclaimer), and complete the 7-point human verification checklist in [`affiliation_and_submission.md`](affiliation_and_submission.md). |

---

## 2. Final Readiness Decision

The manuscript, figures, supplementary tables, provenance logs, and audit documentation on branch `audit/biorxiv-pre-submission-final` are **scientifically consistent, editorially finalized, and ready for immediate human review and manual upload to the bioRxiv submission portal**.
