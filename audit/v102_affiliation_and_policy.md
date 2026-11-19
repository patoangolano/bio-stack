# bioRxiv Readiness, Affiliation, and Policy Compliance Report (`v1.0.2`)

**Audit Date:** 2026-11-19 (branch `audit/biorxiv-final-v1.0.2`)  
**Target Manuscript:** *Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles*  
**Author:** Matheus Lima (ORCID: `0009-0000-3557-8337`; Email: `limatheus458@gmail.com`)  
**Baseline Release Audited:** `enamel-relics-v1.0.1` (Zenodo Version DOI: `10.5281/zenodo.24360150`; Concept DOI: `10.5281/zenodo.24358066`; `v1.0.0` preserved at `10.5281/zenodo.24358067`)

---

## 1. Scope Distinction: Administrative/Documentary Readiness vs. Scientific Peer Review

This report evaluates **documentary, technical, and administrative compliance** with bioRxiv / openRxiv screening criteria. Per strict scientific audit standards:
- **No numerical probability of bioRxiv acceptance is assigned**, as screening decisions depend on human editorial discretion at openRxiv (`https://openrxiv.org/approach-to-ai/`).
- **Administrative screening readiness** (format, disclosures, data availability, absence of placeholders, verifiable human accountability, organizational affiliation) is strictly distinguished from **scientific merit and formal peer review**.

---

## 2. Official bioRxiv / openRxiv Policy Audit (Primary Sources Retrieved 2026-11-19)

Primary policy statements from `openrxiv.org` are archived verbatim in `audit/evidence/f/`:

1. **AI & Human Accountability Policy (20 October 2026 — `https://openrxiv.org/approach-to-ai/`, `audit/evidence/f/ai.txt`):**
   - **Principle 1 ("Humans must be accountable"):** *"AI is a tool that researchers use. It is humans who are accountable for what it produces. All bioRxiv and medRxiv authors must therefore be verifiable human beings and affiliated with an organization to which we can turn if there are allegations of fraud or misconduct."*
   - **Principle 2 ("Human judgment is required"):** Authors using AI must direct, review, and stand behind the work they submit; submissions showing unreviewed AI fabrication (e.g., hallucinated references or unverified claims) are subject to decline.
   - **Principle 3 ("Balancing benefits and drawbacks"):** openRxiv does not penalize transparent AI use for coding, analysis, or English drafting.
   - **Principle 4 ("Humans make the final decisions"):** Human screening staff make final decisions; authors may appeal.
   - **Principle 5 ("Do no harm"):** Excludes dual-use research of concern (DURC), pathogens, or toxins (not applicable to comparative genomics of vertebrate tooth-enamel structural genes).
2. **Trust Signals & Accountability Policy (29 May 2026 — `https://openrxiv.org/trust-signals-and-accountability/`, `audit/evidence/f/trust.txt`):**
   - Requires an organizational affiliation so that disputes regarding authorship, data ownership, or research misconduct can be referred to an organization with ethical oversight.
   - Encourages ORCID identifiers, permanent code/data archival (Zenodo DOIs, SHA-256 manifests), and transparent provenance.
3. **Funder Metadata (`https://openrxiv.org/updating-funder-information/`, 3 Sep 2026) & Licensing (`https://openrxiv.org/preprint-licensing-updates/`, 20 May 2026):**
   - Explicit no-funding statements are valid; CC BY 4.0 licensing is recommended and matches the repository's `.zenodo.json` / `CITATION.cff`.

---

## 3. Affiliation Analysis: Centro Universitário Senac (EAD) vs. Independent Research

- **Factual Context (`CONFIRMADO`):** The author is an undergraduate student enrolled in *Tecnologia em Análise e Desenvolvimento de Sistemas* (distance-learning modality, EAD) at **Centro Universitário Senac — Campus Santo Amaro, São Paulo, Brazil**. The research was conducted independently without institutional supervision, funding, or formal endorsement by Senac.
- **Policy Tension (`PENDENTE DE DECISÃO HUMANA`):**
  1. Listing `"Independent Researcher, São Paulo, Brazil"` (**Option A**, currently compiled in `manuscript_biorxiv.pdf`) is 100% factual and avoids any unauthorized use of Centro Universitário Senac's name as a research sponsor, but faces an administrative screening risk under openRxiv's 20 October 2026 requirement for an organizational affiliation with misconduct oversight.
  2. Listing `"Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil"` (**Option B**, accompanied by an explicit non-endorsement/no-supervision disclaimer) identifies the educational institution to which the student is answerable under student academic-integrity rules, **but must not be adopted without prior written confirmation from Senac**.
- **Rule Enforced in `v1.0.2`:** The main affiliation line remains **Option A** (`Independent Researcher, São Paulo, Brazil`) with a factual *Author Note* disclosing the author's EAD student enrollment at Centro Universitário Senac — Santo Amaro and explicitly stating that the study received no institutional supervision, funding, or formal endorsement.

---

## 4. Revision of the Generative AI Declaration (`CORRIGIDO` in `v1.0.2`)

- **Audit Finding:** In `v1.0.1`, the AI declaration stated that *"manual inspection of alignments, control outputs, and audit reports... were performed by the human author"*. Because sequence alignment parsing, splice-border verification, permutation benchmarks, and ancestral state reconstructions were executed programmatically by AI-assisted audit scripts directed by the human author rather than by manual base-by-base human inspection of every alignment file, that phrasing overstated manual human verification.
- **Corrected Declaration in `v1.0.2` (`submission_biorxiv/manuscript.typ` & `manuscript.md`):**
  > *"During the preparation and auditing of this work, agentic AI coding and research assistants (Anthropic Claude and Google Gemini workflows) were used extensively under the author's direction to implement and execute bioinformatics pipelines, run statistical permutation and neutral-evolution simulations, programmatically audit sequence alignments, splice junctions, and database accessions, verify primary bibliographic metadata, generate figures and supplementary tables, and draft and revise manuscript text. The human author conceived the project scope, defined the analytical and audit requirements, reviewed the generated code, audit reports, and manuscript drafts, and accepts full scientific and ethical accountability for the integrity and accuracy of the final work."*

---

## 5. Readiness Classification Matrix (`PASSOU` / `FALHOU` / `PENDENTE` / `NÃO VERIFICADO`)

| Requirement | Status | Verification & Evidence |
|---|---|---|
| Scope: Empirical biological research (New Results, non-clinical, non-DURC) | **PASSOU** | Comparative genomics of 18 sauropsid reference genomes; zero pathogens, toxins, or biosecurity risks. |
| Verifiable human author & ORCID (`0009-0000-3557-8337`) | **PASSOU** | Verified in `manuscript.typ`, `manuscript.md`, `.zenodo.json`, and `CITATION.cff`. |
| Variant nomenclature & reference vs. ancestral codon separation (Phase A) | **PASSOU** | All 23 variants in `Table_S5a` audited (`XP_019346377.1` / `XP_019346386.1` reference residues separated from reconstructed avian ancestral states; codons 601 and 1003 clarified). |
| Synapomorphy vs. shared-lesion classification with missing sister-clade coverage (Phase B) | **PASSOU** | `AMBN` Exon 4 16-nt deletion (Palaeognathae) and 1-nt insertion + `AG→AC` (Chemyd+Chrpic) and `ENAM` codon 1022 (`p.Ser1022Ter`) explicitly qualified by sister/outgroup coverage in Table 3, `Table_S5a`, and text. |
| Topological sensitivity (`T1`/`T2`/`T3`) & formal IQ-TREE 3 tests (Phase C) | **PASSOU** | `Table_S9` documents 37.5% (`T1`), 62.5% (`T2`), and 68.8% (`T3`) parsimony incongruence (31.2%, 56.2%, 62.5% under reversible ML ASR) + 10,000-replicate RELL KH/SH/wSH/ELW/AU tests + single-locus linkage caveat. |
| Neutral K80 simulation coverage-mask correction & sensitivity grid (Phase D) | **PASSOU** | `Table_S8` documents why unmasked 878-codon simulation predicted 26–29 PTCs whereas empirical coverage masking predicts median 15–19 shared PTCs (matching 16 observed, two-sided *p* = 0.56–0.95) and 40–50% incongruence. |
| Bibliographic integrity (49/49 primary references) | **PASSOU** | 49/49 verified via Crossref (`HTTP 200`), 47/49 Europe PMC, 0 retractions/concerns, 1:1 citation match. |
| Factual Generative AI disclosure (no overstated manual verification claims) | **PASSOU** | Updated in `manuscript.typ` and `manuscript.md` to reflect actual human direction/review and agentic execution. |
| Format, ≤250-word Abstract, continuous line/page numbers, visual PDF check | **PASSOU** | Compiled with `typst 0.15.1`, visually inspected across figures, tables, headers, and references; 0 placeholders; no premature "posted on bioRxiv" header. |
| Versioning & Zenodo DOI policy (`v1.0.2` candidate) | **PASSOU** | Preserves tags `enamel-relics-v1.0.0` (`10.5281/zenodo.24358067`) and `enamel-relics-v1.0.1` (`10.5281/zenodo.24360150`), cites Concept DOI `10.5281/zenodo.24358066`, and does not reuse `24360150` as a fake `v1.0.2` DOI. |
| Uncorrected scientific or technical inconsistencies | **FALHOU: 0** | Zero known uncorrected technical or scientific inconsistencies remain in `v1.0.2`. |
| Organizational affiliation confirmation (Senac & openRxiv screening) | **PENDENTE** | Requires author dispatch of Draft Emails 1 & 2 (Section 6 below) and written response before choosing Option B or submitting under Option A. |
| PR merge, `enamel-relics-v1.0.2` release minting, and manual bioRxiv upload | **PENDENTE** | Requires explicit human review and action (`audit/v102/human_approval_checklist.md`). |
| Out-of-scope wet-lab / external assembly validations (declared in Section 3.5) | **NÃO VERIFICADO** | Second conspecific *Apteryx rowi* assembly (unavailable on NCBI), chromosome-level Pleurodira turtle sampling, $d_N/d_S$ branch modeling, and embryonic oral RNA-seq. Explicitly declared as Limitations. |

---

## 6. Draft Consultations for Human Dispatch (DO NOT SEND AUTOMATICALLY)

### 6.1. Consulta Formal ao Centro Universitário Senac (Português)

> **Para:** Coordenação do Curso de Tecnologia em Análise e Desenvolvimento de Sistemas (EAD) / Ouvidoria / Setor de Integridade Acadêmica — Centro Universitário Senac (Campus Santo Amaro)  
> **Assunto:** Consulta sobre forma correta de declaração de vínculo estudantil em preprint científico independente (bioRxiv)
>
> Prezados(as),
>
> Meu nome é **Matheus Lima**, estudante regularmente matriculado no curso superior de Tecnologia em Análise e Desenvolvimento de Sistemas (modalidade EAD), vinculado ao Campus Santo Amaro (ORCID: `0009-0000-3557-8337`).
>
> Conduzi, de maneira totalmente independente, extracurricular e sem financiamento ou orientação institucional, uma pesquisa computacional em bioinformática e genômica comparativa intitulada *"Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles"*, cujo código-fonte e dados auditados encontram-se publicamente arquivados no GitHub e no Zenodo (`https://doi.org/10.5281/zenodo.24358066`).
>
> Pretendo submeter o artigo como preprint ao servidor **bioRxiv** (mantido pela organização openRxiv). Como as diretrizes de governança do bioRxiv/openRxiv (atualizadas em maio e outubro de 2026) solicitam que autores informem uma organização de vínculo e o respectivo canal competente para questões de integridade acadêmica, consulto formalmente:
> 1. O Centro Universitário Senac autoriza a indicação estritamente factual do meu vínculo estudantil no cabeçalho do manuscrito — na forma *"Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil"* — acompanhada de nota expressa declarando que o estudo foi realizado de forma independente, sem orientação, financiamento, responsabilidade técnica ou endosso institucional do Senac?
> 2. Em caso afirmativo, qual setor e endereço eletrônico institucional (por exemplo, Coordenação Acadêmica, Ouvidoria ou Comitê de Ética/Integridade) pode ser indicado caso a equipe de triagem do bioRxiv deseje confirmar meu vínculo discente?
>
> Caso a instituição prefira que seu nome não figure na linha principal de afiliação de pesquisas extracurriculares independentes, manterei a afiliação principal como *"Independent Researcher, São Paulo, Brazil"*.
>
> Atenciosamente,  
> **Matheus Lima**  
> E-mail: `limatheus458@gmail.com` · ORCID: `https://orcid.org/0009-0000-3557-8337`

### 6.2. Pre-Submission Eligibility Inquiry to bioRxiv / openRxiv (`bioRxiv@openrxiv.org`, English)

> **To:** `bioRxiv@openrxiv.org`  
> **Subject:** Pre-submission inquiry: Affiliation eligibility for independent undergraduate computational biology manuscript
>
> Dear bioRxiv / openRxiv Screening and Editorial Team,
>
> I am writing to seek clarification on affiliation formatting prior to submitting a comparative genomics manuscript titled *"Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles"* (Subject Area: *Evolutionary Biology*; Article Type: *New Results*).
>
> I am an undergraduate student enrolled in the Systems Analysis and Development program (distance learning) at Centro Universitário Senac — Santo Amaro, São Paulo, Brazil (ORCID: `0009-0000-3557-8337`). This computational study was conceived and executed independently without university funding or faculty supervision. All scripts, deterministic seeds, SHA-256 manifests, primary bibliographic verifications, and raw/processed alignments are openly archived on GitHub (`https://github.com/patoangolano/bio-stack`) and Zenodo (`https://doi.org/10.5281/zenodo.24358066`), and the manuscript includes a transparent disclosure of agentic AI coding/research assistance in accordance with openRxiv's October 20, 2026 policy.
>
> To ensure full compliance with openRxiv's accountability and affiliation policies, could you please advise:
> 1. Whether listing `"Independent Researcher, São Paulo, Brazil"` as the primary affiliation—accompanied by a factual Author Note stating my undergraduate enrollment at Centro Universitário Senac — Santo Amaro (with proof of enrollment available upon request) and clarifying that the research was conducted independently without institutional supervision or endorsement—is acceptable for bioRxiv screening; OR
> 2. Whether bioRxiv requires the university name in the primary affiliation line (`"Undergraduate Student, Centro Universitário Senac — Santo Amaro, São Paulo, Brazil"`, with the same non-endorsement disclaimer) once my institution confirms the appropriate academic-integrity contact office.
>
> Thank you very much for your guidance.
>
> Sincerely,  
> **Matheus Lima**  
> Email: `limatheus458@gmail.com` · ORCID: `https://orcid.org/0009-0000-3557-8337`
