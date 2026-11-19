# Phase F Audit Report — Author Affiliation Options, bioRxiv Portal Instructions, and Pre-Submission Verification Checklist

**Status:** `PENDENTE` (Awaiting final human author choice between Option A and Option B on the bioRxiv submission portal)  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final`  
**Author:** Matheus Lima (ORCID: `0009-0000-3557-8337`; Correspondence Email: `limatheus458@gmail.com`)

---

## 1. F1 — Current Author Context & Official openRxiv / bioRxiv Policy Analysis

### 1.1. Factual Status of the Author's bioRxiv Account (`CONFIRMADO`)
- The author is an undergraduate student regularly enrolled in *Tecnologia em Análise e Desenvolvimento de Sistemas* (distance-learning modality, EAD) at **Centro Universitário Senac — Campus Santo Amaro, São Paulo, Brazil**.
- The research was conceived and executed independently by the author, without institutional supervision, funding, or formal endorsement from Centro Universitário Senac.
- **Portal Progress:** The author has already created an account on bioRxiv using their institutional Senac email address, accessed the **Author Area**, and initiated the manuscript submission form.
- **Important Distinction:** Successfully creating a bioRxiv account with an institutional email and starting the submission form confirms **technical portal access**, but **does not pre-approve the manuscript's affiliation line during manual post-submission screening** by openRxiv staff (`Principle 4: "Humans make the final decisions"`, `audit/evidence/f/ai.txt`).

### 1.2. Express Policy Requirements vs. Prudential Recommendations
Auditing the verbatim primary documents (`audit/evidence/f/ai.txt`, 20 October 2026; `audit/evidence/f/trust.txt`, 29 May 2026) establishes the exact boundary between **express openRxiv rules** and **prudential academic precautions**:

1. **Express openRxiv Rule (Principle 1, 20 October 2026):**
   > *"All bioRxiv and medRxiv authors must therefore be verifiable human beings and affiliated with an organization to which we can turn if there are allegations of fraud or misconduct."*
2. **No Express bioRxiv Requirement for Prior Written Institutional Authorization:**
   - Neither bioRxiv's Submission Guide nor openRxiv's October 20, 2026 policy demands a signed letter of institutional authorization to upload a preprint, provided the author's affiliation is **factual** (the author is genuinely enrolled at that institution) and does **not misrepresent institutional endorsement or funding**.
   - **Prudential Precaution vs. Mandatory Rule:** Seeking written confirmation from Centro Universitário Senac before placing the university's name on the primary affiliation line (**Option B**) is a **prudential recommendation** to avoid any internal misunderstanding with the university regarding an unendorsed extracurricular study—**not** an external rule imposed by bioRxiv.

---

## 2. F2 — Two Factual Affiliation Formulations Prepared for Author Decision

Both options below are 100% honest and explicitly disclaim institutional supervision, funding, and endorsement. **Option A** is currently compiled in `submission_biorxiv/manuscript_biorxiv.pdf`; if the author prefers **Option B** (to match the institutional Senac email used to register the bioRxiv account and satisfy Principle 1 directly in the primary affiliation field), the exact 3-line replacement in `manuscript.typ` and `manuscript.md` is provided below.

### Option A — Primary Affiliation as Independent Researcher + Factual Senac Student Author Note *(Currently Compiled in PDF)*

- **Manuscript Header (`manuscript.typ` lines 102–107):**
  - **Line 1 (Primary Affiliation):** `1 Independent Researcher, São Paulo, Brazil`
  - **Line 2 (Correspondence & ORCID):** `* Correspondence: limatheus458@gmail.com · ORCID: 0009-0000-3557-8337`
  - **Line 3 (Author Note):** *Author Note: The author is an undergraduate student enrolled in the Systems Analysis and Development program (distance learning) at Centro Universitário Senac (Santo Amaro campus, São Paulo, Brazil); this study was conducted independently with no institutional supervision, funding, or formal endorsement.*
- **Administrative Trade-off:** Completely avoids placing Senac's name in the primary affiliation index without prior consultation, while still disclosing the author's verifiable student enrollment on the title page. Because the author registered the bioRxiv account with an institutional Senac email address, screeners can verify the student bond if they read the *Author Note*, though automated metadata parsers that inspect only the portal's `Institution` field may flag `"Independent Researcher"` for manual review.

### Option B — Factual Student Affiliation at Centro Universitário Senac + Explicit Non-Endorsement Disclaimer *(Alternative Ready for Author Selection)*

- **Manuscript Header (if selected by the author):**
  - **Line 1 (Primary Affiliation):** `1 Undergraduate Student, Systems Analysis and Development Program (EAD), Centro Universitário Senac — Santo Amaro, São Paulo, Brazil`
  - **Line 2 (Correspondence & ORCID):** `* Correspondence: limatheus458@gmail.com · ORCID: 0009-0000-3557-8337`
  - **Line 3 (Author Note / Disclaimer):** *Author Note: This study was conducted independently as an extracurricular research project, without institutional supervision, funding, or formal endorsement by Centro Universitário Senac.*
- **Administrative Trade-off:** Aligns 1:1 with the author's institutional Senac email address registered on bioRxiv and satisfies Principle 1's organizational affiliation requirement directly in the primary metadata field, while explicitly stating that the work is an independent student project without institutional supervision, funding, or endorsement.

---

## 3. F3 — Field-by-Field Guide for Filling Out the bioRxiv Submission Form

When completing the active submission form in the bioRxiv Author Area, enter the following exact values:

| bioRxiv Portal Field | Exact Value to Enter | Notes |
|---|---|---|
| **Article Type** | `New Results` | Empirical comparative genomics study across 18 RefSeq assemblies. |
| **Subject Area (Primary)** | `Evolutionary Biology` | Secondary subject (if prompted): `Genomics`. |
| **Manuscript Title** | `Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles` | Match `manuscript_biorxiv.pdf` verbatim. |
| **Author Name & ORCID** | `Matheus Lima` · ORCID: `0009-0000-3557-8337` | Link/authenticate ORCID in the bioRxiv author profile. |
| **Author Institution / Affiliation** | **If keeping Option A:** `Independent Researcher (Undergraduate Student, Centro Universitário Senac — Santo Amaro), São Paulo, Brazil`<br>**If choosing Option B:** `Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil` | Do **not** include student ID numbers (*RA*), CPF, or account passwords anywhere in the manuscript or metadata. |
| **Abstract** | Copy verbatim from [`submission_biorxiv/manuscript.md`](../../submission_biorxiv/manuscript.md) (lines 14) | Verified **248 words** ($\le 250$-word bioRxiv limit). |
| **Main Manuscript File** | Upload [`submission_biorxiv/manuscript_biorxiv.pdf`](../../submission_biorxiv/manuscript_biorxiv.pdf) | 21 pages, embedded 320-DPI figures, continuous line and page numbers. |
| **Supplementary Material** | Upload `Table_S0` through `Table_S9` + `MANIFEST_SHA256.tsv` from [`submission_biorxiv/supplementary/`](../../submission_biorxiv/supplementary/) (or the zip archive `enamel_relics_supplementary_tables.zip`) | 14 `.tsv` tables + SHA-256 manifest. |
| **Data & Code Availability / External Links** | GitHub: `https://github.com/patoangolano/bio-stack`<br>Zenodo Concept DOI: `https://doi.org/10.5281/zenodo.24358066`<br>Zenodo `v1.0.2` DOI: `https://doi.org/10.5281/zenodo.24368656` | Enter in the External Data / Code Repository fields on the portal. |
| **Funding Statement** | `This research received no specific grant or external funding from any funding agency in the public, commercial, or not-for-profit sectors.` | Select `"No funding"` in the standardized funder metadata widget. |
| **Competing Interests** | `The author declares no competing financial or non-financial interests.` | Select `"No competing interests"`. |
| **Cover Letter / Comments to Screeners (Optional Box)** | *"This comparative genomics study of 18 public sauropsid RefSeq genomes was conducted independently by an undergraduate student enrolled at Centro Universitário Senac — Santo Amaro (ORCID 0009-0000-3557-8337) without institutional funding or supervision. All analysis scripts, deterministic seeds, SHA-256 manifests, primary bibliographic checks (49/49 verified references), and supplementary tables (Tables S0–S9) are openly archived on GitHub (https://github.com/patoangolano/bio-stack) and Zenodo (https://doi.org/10.5281/zenodo.24358066), and the manuscript includes a full Declaration of Generative AI and AI-Assisted Technologies in compliance with openRxiv's October 20, 2026 policy."* | Strongly recommended if submitting under Option A, as human screeners read this note directly. |

---

## 4. F4 — Human Author Scientific Accountability & Verification Checklist

The **Declaration of Generative AI and AI-Assisted Technologies** (`manuscript.typ` lines 388–390) transparently discloses that agentic AI workflows (Anthropic Claude and Google Gemini) were used under the author's direction for pipeline implementation, statistical simulations, programmatic alignment and database auditing, bibliographic verification, figure generation, and manuscript drafting, without claiming manual base-by-base alignment inspections that were executed programmatically.

Before clicking **Submit** on bioRxiv, the human author should personally review and confirm the following 7 core claims against their cited tables and figures:

- [ ] **1. Dual-database absence in chicken `GRCg7b` & 3 chrZ relics (Table 1, Figure 1B, `Table_S1`, `Table_S2a–b`):** Verify that *ΨENAM* Exon 8 (`1,985-bp` union) and *ΨAMBN* Exon 3 (`138 bp`) are reciprocal `BLASTx`-confirmed, whereas *ΨENAM* Exon 6 (`47-bp` core, truncated at codon 4) is explicitly reported as nucleotide- and synteny-confirmed (*not* reciprocal `BLASTx`-confirmed).
- [ ] **2. Four Gnomon `LOW QUALITY PROTEIN` models (`Table_S7`, Section 2.1):** Verify the distinction between the 3 chromosome-level pseudogenes (*C. picta* *ENAM*, *C. mydas* *AMBN*, *S. camelus* *AMELX*) and the 1 scaffold-level single-assembly lesion (*A. rowi* *AMELX*).
- [ ] **3. Disordered-protein homology bias (Table 2, Figure 2, `Table_S3a–c`):** Verify that equidistant *Anolis* queries (`1,000` permutations) show intact avian bone IDPs (*IBSP*, *SPP1*) losing `25.5–51.5%` of excess bitscore—overlapping decayed avian *ENAM*/*AMBN*—while globular *LPL* scores `112.6%` of crocodylian levels, and that the marginal `Strhab` `ENAM_Acar` score ($p = 0.024$) is explicitly labeled exploratory and unadjusted.
- [ ] **4. Level 3 `miniprot` coding lesions (Figure 3A, `Table_S4`):** Verify `0.00` lesions/1,000 aa across `18/18` *LPL* controls (`133/133` `GT-AG` introns), `0.11`/1,000 aa across `15` toothed enamel alignments, and `16.98`/1,000 aa across edentulous enamel pseudogenes.
- [ ] **5. Reference vs. reconstructed avian ancestral codons (Table 3, Figure 3B–C, `Table_S5a`):** Verify the explicit distinction at *ENAM* codons `601` (`TTT` [Phe] $\to$ Neornithes `TAT` [Tyr] $\to$ Palaeognathae `TAA`), `1003` (`TGG` [Trp] vs. site-independent ML `TGC` [Cys] $\to$ Psittacopasserae `TGA`), and `1022` (`Anapla` uncovered; parsimony equiparsimony vs. reversible ML Neognathae `TAG` $P_{\text{stop}} = 0.926$).
- [ ] **6. Crown-Neornithes coverage limits & `T1`–`T3` topological sensitivity (Figure 3B, Table 3, `Table_S9`):** Verify that `0/878` *ENAM* Exon 8 codons are covered in all 7 birds (`16.3%` neutral ancestral-stop detection power), so crown-Neornithes-wide synapomorphies are declared *not detected / unresolved* rather than biologically absent, and that PTC incongruence is reported across `T1` (`31.2–37.5%`), `T2` (`56.2–62.5%`), and `T3` (`62.5–68.8%`).
- [ ] **7. Model-conditional K80 simulation interpretation (`Table_S8`):** Verify that the coverage-masked K80 simulation (`16.0–19.0` median shared PTCs, `40.0–44.4%` incongruence) is presented as statistical compatibility conditioned on the K80 model and empirical coverage mask rather than causal proof.
