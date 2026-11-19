# bioRxiv Readiness, Policy Compliance, and Affiliation Report

**Audit Date:** 2026-11-19 (branch `audit/biorxiv-revision-v1.0.1`)
**Target Manuscript:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles*
**Author:** Matheus Lima (ORCID: `0009-0000-3557-8337`; Email: `limatheus458@gmail.com`)

---

## 1. Current bioRxiv / openRxiv Policy Audit (Primary Sources Retrieved 2026-11-19)

Primary policy pages were retrieved directly from `openrxiv.org` (parent non-profit operating bioRxiv and medRxiv) and saved verbatim in `audit/evidence/f/`:

1. **AI & Human Accountability Policy (20 October 2026):**
   - **Title:** *How we're approaching AI at bioRxiv and medRxiv*
   - **URL:** `https://openrxiv.org/approach-to-ai/` (saved in `audit/evidence/f/ai.html` / `ai.txt`)
   - **Key Binding Clauses:**
     - **Principle 1 ("Humans must be accountable"):** *"AI is a tool that researchers use. It is humans who are accountable for what it produces. All bioRxiv and medRxiv authors must therefore be verifiable human beings and affiliated with an organization to which we can turn if there are allegations of fraud or misconduct."*
     - **Principle 2 ("Human judgment is required"):** *"AI is now commonly used in research design, execution, analysis, and reporting. However, if an author uses AI, they must shape, review, stand behind and have intellectual input into the interpretation of the work they submit. When screening identifies evidence that authors have not critically reviewed AI-generated material, such as large numbers of fabricated references or entirely AI-generated discussion/conclusions, a submission may receive further review and we reserve the right to decline the work."*
     - **Principle 3 ("Balancing the benefits and drawbacks"):** Explicitly affirms that openRxiv does not penalize AI use (recognizing utility for drafting, especially for non-native English speakers), while reserving the right to decline submissions failing research-integrity screening.
     - **Principle 4 ("Humans make the final decisions"):** Automated tools inform screening, but human staff make final decisions; authors can appeal.
     - **Principle 5 ("Do no harm"):** Responsible screening excludes potential dual-use research of concern (DURC) or in silico biosecurity risks. (Not applicable here: basic comparative genomics of benign tooth-enamel structural genes in vertebrates.)
2. **Trust Signals & Accountability Policy (29 May 2026):**
   - **Title:** *Trust Signals and Accountability* (by Richard Sever & Tracy K. Teal)
   - **URL:** `https://openrxiv.org/trust-signals-and-accountability/` (saved in `audit/evidence/f/trust.html` / `trust.txt`)
   - **Key Binding Clause:** *"Sadly research misconduct does occur, and openRxiv frequently has to respond in disputes about authorship, data ownership or intellectual property. These must be arbitrated by the organization with which an individual is affiliated. This is one reason openRxiv requires authors to supply an affiliation: we need an organization responsible for ethical oversight to which we can refer allegations of fraud/misconduct."*
   - Also emphasizes ORCIDs as researcher disambiguation (not yet standalone identity verification), data/code audit trails (hashes, reproducible pipelines), and standardized ROR funder metadata.
3. **Funder Metadata Update (3 September 2026):** `https://openrxiv.org/updating-funder-information/` (`audit/evidence/f/fund.txt`) — standardized ROR funder metadata supported, including post-posting updates; explicit no-funding declarations remain valid.
4. **Licensing Updates (20 May 2026):** `https://openrxiv.org/preprint-licensing-updates/` (`audit/evidence/f/lic.txt`) — authors retain copyright and select from Creative Commons licenses (CC BY 4.0 recommended for open science and compatible with Zenodo CC-BY-4.0).
5. **Contact Email:** `bioRxiv@openrxiv.org` (confirmed via `biorxiv.org` submission/contact guidelines).

---

## 2. Requirement-by-Requirement Compliance Matrix

| Requirement | Status | Evidence & Notes |
|---|---|---|
| Empirical biological research (New Results, non-clinical, non-DURC) | **PASSOU** | Original 18-genome comparative genomic & paleogenomic analysis; zero pathogens/toxins/DURC. |
| Verifiable human author & ORCID | **PASSOU** | Matheus Lima, ORCID `0009-0000-3557-8337`, email `limatheus458@gmail.com`. |
| **Organizational affiliation with ethical oversight** (20 Oct 2026 & 29 May 2026 policies) | **PENDENTE** | **Critical screening risk.** Current manuscript lists `"Independent Researcher, Brazil"`. openRxiv explicitly requires affiliation with an organization capable of ethical/misconduct oversight. Author is enrolled in *Tecnologia em Análise e Desenvolvimento de Sistemas* (EAD) at Centro Universitário Senac — Santo Amaro, São Paulo, Brazil, with no documented institutional supervision, funding, or authorization to list Senac as research sponsor. Requires human decision & institutional/bioRxiv clarification (see Sections 3–4). |
| No premature "bioRxiv Preprint" posted-status label | **PASSOU** *(corrected in v1.0.1)* | v1.0.0 running header displayed `bioRxiv Preprint` despite no posted bioRxiv DOI/record. Replaced in `v1.0.1` with `Preprint Manuscript — Draft for Submission` (no claim of prior posting). |
| Abstract self-contained and ≤ 250 words | **PASSOU** | Verified word count ≤ 250 words in `manuscript.typ` / `manuscript.md` (checked by automated test suite). |
| Complete IMRaD + Limitations + Data/Code + Funding + Competing Interests | **PASSOU** | All sections present, factual, and synchronized across Typst, Markdown, and PDF. |
| Transparent Declaration of Generative AI | **PASSOU** | Explicitly discloses Anthropic Claude / Google Gemini coding & editing assistance and states human accountability, manual alignment/control inspection, and responsibility per openRxiv 20 Oct 2026 Principle 2. |
| Zero fabricated references (0 retractions) | **PASSOU** | 100% of cited references independently verified against Crossref and Europe PMC (`audit/bibliography_report.md`). |
| Continuous line numbering & page numbering | **PASSOU** | Enforced in Typst `0.15.1` (`#set par.line(...)` and `Page X of Y` footer). |
| Open data/code + permanent archive (Zenodo DOI) | **PENDENTE** *(v1.0.1 release)* | `enamel-relics-v1.0.0` archived at `10.5281/zenodo.24358067` (Concept DOI `10.5281/zenodo.24358066`), but v1.0.0 lists creator `patoangolano` and predates v1.0.1 corrections. `.zenodo.json` & `CITATION.cff` created for `v1.0.1`; manuscript cites Concept DOI `10.5281/zenodo.24358066` (resolving to latest release) while noting `v1.0.0` (`24358067`) as archival baseline. |
| Author formal submission declarations on bioRxiv portal | **PENDENTE** | Must be completed manually by the human author in the bioRxiv submission portal. |

---

## 3. Affiliation Options for Human Decision (Zero Unauthorized Institutional Claims)

Because Matheus Lima is an enrolled undergraduate student at **Centro Universitário Senac — Santo Amaro** (EAD, *Tecnologia em Análise e Desenvolvimento de Sistemas*) conducting independent research **without documented institutional supervision, funding, or formal endorsement**, the manuscript **must not** imply that Senac supervised, funded, or approved this study. Two transparent options are presented for the author's explicit choice:

### Option A — Conservative Default (Currently Applied in Draft `v1.0.1` Pending Authorization)
- **Byline affiliation:**
  `Independent Researcher, São Paulo, Brazil`
- **Mandatory Author Note / Affiliation Clarification (in Acknowledgments/Author Notes):**
  *"The author is an undergraduate student enrolled in the Systems Analysis and Development program (distance learning modality) at Centro Universitário Senac (Santo Amaro campus, São Paulo, Brazil). This study was conceived and conducted independently and received no institutional supervision, funding, or formal endorsement from Centro Universitário Senac."*
- **Pros:** 100% truthful; zero risk of misrepresenting institutional endorsement or violating Senac branding/academic policies.
- **Cons:** High risk of screening hold or decline by bioRxiv under the 20 October 2026 policy requiring an organizational affiliation for misconduct oversight unless bioRxiv confirms eligibility in advance (via Draft Email 2 below).

### Option B — Institutional Student Identification (Conditional on Written Confirmation from Senac & bioRxiv)
- **Byline affiliation:**
  `Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil`
  *(or `Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil`)*
- **Mandatory Author Note / Disclaimer:**
  *"The author is an undergraduate student at Centro Universitário Senac — Santo Amaro; this work was conducted independently and does not represent institutional research, supervision, funding, or endorsement by Centro Universitário Senac."*
- **Pros:** Directly identifies the academic institution to which the student belongs, potentially satisfying openRxiv's organizational-oversight criterion if Senac confirms its ethics/ombudsman office can receive external integrity inquiries regarding enrolled students.
- **Cons:** **Must NOT be used without written confirmation from Senac** (coordination / ethics / research office) permitting use of the institution's name for an independent extracurricular manuscript.

> **Recommendation:** Send **both** draft enquiries below *before* clicking submit on bioRxiv. Keep **Option A** in the manuscript until written responses are received.

---

## 4. Formal Draft Enquiries (DO NOT SEND AUTOMATICALLY — For Author Review & Dispatch)

### Draft 1 (Portuguese) — Consulta Formal ao Centro Universitário Senac (Coordenação do Curso / Pró-Reitoria de Pesquisa / Ouvidoria Acadêmica)

> **Assunto:** Consulta sobre identificação de vínculo estudantil em manuscrito científico independente (preprint bioRxiv)
>
> Prezados(as) responsáveis pela Coordenação do Curso de Tecnologia em Análise e Desenvolvimento de Sistemas (EAD) / Núcleo de Pesquisa e Ética Acadêmica do Centro Universitário Senac — Campus Santo Amaro,
>
> Meu nome é **Matheus Lima**, aluno regularmente matriculado no curso superior de Tecnologia em Análise e Desenvolvimento de Sistemas (modalidade EAD), vinculado ao Campus Santo Amaro (ORCID: `0009-0000-3557-8337`).
>
> Desenvolvi, de forma totalmente independente e sem financiamento externo, um estudo computacional em genômica comparativa e bioinformática intitulado *"Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles"*, com todo o código e dados públicos arquivados no GitHub e no Zenodo.
>
> Pretendo submeter o manuscrito como preprint ao servidor **bioRxiv** (mantido pela organização sem fins lucrativos openRxiv). Em razão da política editorial atual do bioRxiv (publicada em 20 de outubro de 2026), exige-se que autores indiquem uma afiliação institucional verificável para fins de responsabilidade ética e integridade científica.
>
> Venho, com total transparência, solicitar orientação institucional sobre dois pontos:
> 1. É permitido indicar meu vínculo discente de forma estritamente factual — por exemplo, *"Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil"* — acompanhado de nota explícita declarando que o trabalho é independente e não contou com orientação, financiamento ou endosso institucional do Senac?
> 2. Caso a instituição autorize essa identificação discente, qual setor/e-mail institucional deveria constar como contato formal caso o bioRxiv necessite verificar meu vínculo de matrícula ou encaminhar eventuais questões de integridade acadêmica?
>
> Caso a orientação institucional seja não utilizar o nome do Centro Universitário Senac em produções extracurriculares independentes, respeitarei integralmente essa diretriz e manterei apenas a designação *"Independent Researcher, São Paulo, Brazil"*.
>
> Agradeço antecipadamente pela atenção e coloco-me à disposição para enviar a cópia integral do manuscrito e o comprovante de matrícula.
>
> Atenciosamente,
> **Matheus Lima**
> `limatheus458@gmail.com` · ORCID: `https://orcid.org/0009-0000-3557-8337`

### Draft 2 (English) — Pre-Submission Eligibility Enquiry to bioRxiv / openRxiv (`bioRxiv@openrxiv.org`)

> **Subject:** Pre-submission inquiry on affiliation eligibility for an independently authored undergraduate bioinformatics manuscript
>
> Dear bioRxiv / openRxiv Editorial Team,
>
> I am writing to request clarification on affiliation eligibility prior to submitting an empirical comparative genomics manuscript titled *"Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles"* (target category: *Evolutionary Biology*, Article Type: *New Results*).
>
> I am a verifiable human researcher (Matheus Lima, ORCID `0009-0000-3557-8337`) currently enrolled as an undergraduate student in the Systems Analysis and Development program (distance learning) at Centro Universitário Senac — Santo Amaro, São Paulo, Brazil. This bioinformatics study was designed, executed, and audited independently using public NCBI RefSeq and Ensembl reference genomes, with full code, SHA-256 checksums, deterministic seeds, and provenance logs openly archived on GitHub and Zenodo (`https://doi.org/10.5281/zenodo.24358066`). Generative AI coding/editing tools were used transparently as documented in the manuscript's dedicated AI declaration, and all scientific design, manual sequence/alignment auditing, and verification were performed by me.
>
> In light of openRxiv's May 29, 2026 (*Trust Signals and Accountability*) and October 20, 2026 (*How we're approaching AI at bioRxiv and medRxiv*) policies regarding organizational affiliation and accountability, could you please advise:
> 1. Whether a submission from an enrolled undergraduate student conducting independent computational research is eligible for screening if institutional affiliation is listed as `"Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil"` (subject to my university's written permission) with an explicit disclaimer that the university did not supervise or fund the study; OR
> 2. Whether `"Independent Researcher, São Paulo, Brazil"` with a verified ORCID, institutional enrollment documentation available on request, and full reproducible code/data archival is acceptable for screening under current bioRxiv guidelines.
>
> Thank you for your time and guidance.
>
> Best regards,
> **Matheus Lima**
> Email: `limatheus458@gmail.com`
> ORCID: `https://orcid.org/0009-0000-3557-8337`

---

## 5. Zenodo Metadata Audit & Concept vs. Version DOI Decision

Live queries against `https://doi.org/10.5281/zenodo.24358067`, Zenodo REST API (`/api/records/24358067`), and DataCite API (`https://api.datacite.org/dois/10.5281/zenodo.24358067`) on 2026-11-19 confirmed:

| Field | Live Zenodo / DataCite Value (`24358067`) | Audit Finding |
|---|---|---|
| **Version DOI** | `10.5281/zenodo.24358067` (state: `findable`, HTTP 302 → 200) | Permanently pinned to GitHub tag `enamel-relics-v1.0.0` (commit `a425b2d`). |
| **Concept DOI** | `10.5281/zenodo.24358066` (`IsVersionOf`) | Always resolves to the **latest** published archival version of this repository on Zenodo. |
| **Archived ZIP** | `patoangolano/bio-stack-enamel-relics-v1.0.0.zip` (83,096,375 bytes, MD5 `4943fe7d17410ae4abc570679933ba46`) | Verified by direct download: archived `manuscript_biorxiv.pdf` SHA-256 = `f768fb366ce07965...` (exact match to tag `enamel-relics-v1.0.0`, which still contains the `10.5281/zenodo.XXXXXXX` placeholder). |
| **Creator Metadata** | `[{"name": "patoangolano", "affiliation": null}]`, `nameIdentifiers: []` | **Metadata defect:** GitHub handle used instead of `Lima, Matheus` + ORCID `0009-0000-3557-8337`, because `.zenodo.json` and `CITATION.cff` were absent. |
| **License** | `cc-by-4.0` (SPDX CC BY 4.0) | Appropriate and open-science compliant. |
| **Publication Date** | `2026-11-19` | Matches release publication timestamp (`2026-11-19T03:18:18Z`). |

### Resolution Applied in Branch `audit/biorxiv-revision-v1.0.1`
1. **Added `.zenodo.json` and `CITATION.cff`** at repository root with `Lima, Matheus` and ORCID `0009-0000-3557-8337` (no unauthorized institutional affiliation) so that when the human author creates release `enamel-relics-v1.0.1`, Zenodo will automatically mint a new version DOI with accurate author/ORCID metadata under the same Concept DOI `10.5281/zenodo.24358066`.
2. **DOI Citation Policy in `manuscript.typ` / `manuscript.md`:** Cite the **Concept DOI `https://doi.org/10.5281/zenodo.24358066`** as the primary persistent link (which automatically resolves to `v1.0.1` once released, avoiding the chicken-and-egg self-referential version-DOI mismatch that affected `v1.0.0`), while explicitly noting `enamel-relics-v1.0.0` (`10.5281/zenodo.24358067`) as the initial archival snapshot.
3. **Removed Premature Posted-Status Labels:** Changed running header in `manuscript.typ` from `[bioRxiv Preprint]` to `[Preprint Manuscript — Draft for Submission]` and updated `README_SUBMISSION.md` accordingly.
