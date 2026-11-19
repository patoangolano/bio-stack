# Phase C Audit Report — Final Statistical Interpretation Review (Exploratory Multiplicity, Linked-Site Intervals, and Topological Non-Independence)

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final`  
**Primary Script & Output:** [`scripts/enamel_relics/13_audit_final_neornithes_and_stats.py`](../../scripts/enamel_relics/13_audit_final_neornithes_and_stats.py) · [`runs/enamel_relics/v102/neornithes_detectability_audit.json`](../../runs/enamel_relics/v102/neornithes_detectability_audit.json)

---

## 1. C1 — Multiplicity & Exploratory Status of Marginal Permutation *p*-Values (`p = 0.024` in Section 2.3)

### 1.1. Audit of the 504-Cell Benchmark Grid (`Table_S3c`)
- The Level 2/3 compositional sensitivity benchmark (`Table_S3c_per_gene_seg_null_18sp.tsv`) evaluates **504 cells** (`18` sauropsid species $\times$ `14` protein queries $\times$ `2` `SEG` modes), computing for each cell an empirical one-sided permutation *p*-value:
  $$p_{\text{emp}} = \frac{1 + \#\{\text{shuf} \ge \text{obs}\}}{N + 1}, \quad N = 1{,}000 \text{ (minimum achievable } p_{\text{emp}} = \frac{1}{1001} \approx 0.0010\text{)}.$$
- Inspecting all 504 cells (`neornithes_detectability_audit.json`, `table_s3c_marginal_p_cells`) shows that **35 cells** fall in the marginal range $0.005 < p_{\text{emp}} < 0.10$, including:
  - *Strigops habroptila* (`Strhab`) queried with equidistant *Anolis* `ENAM_Acar` (`-seg yes`: $\text{obs} = 30.8\text{ bits}$, $p_{95,\text{1000}} = 30.0\text{ bits}$, $\text{excess} = 0.80\text{ bits}$, **$p_{\text{emp}} = 0.024$**; `-seg no`: $p_{\text{emp}} = 0.025$).
  - Genuine cross-chromosome negative-control hits such as *Dromaius novaehollandiae* (`Dronov`) queried with `DSPP_Amis` (`-seg yes`: $\text{obs} = 27.3\text{ bits}$, $p_{95,\text{1000}} = 26.9\text{ bits}$, $\text{excess} = 0.40\text{ bits}$, **$p_{\text{emp}} = 0.043$**) and *Tinamus guttatus* (`Tingut`) queried with `DSPP_Amis` (`-seg yes`: **$p_{\text{emp}} = 0.008$**; `-seg no`: **$p_{\text{emp}} = 0.031$**), where `DSPP` does not even reside in the `NRG1–LPL` interval!

### 1.2. Statistical Classification & Revision Applied
1. **Why unadjusted $p_{\text{emp}} = 0.024$ cannot be treated as confirmatory:** Across 504 tests (or even across the 72 edentulous enamel cells), an unadjusted one-sided threshold of $\alpha = 0.05$ expects $\sim 3.6$ to $25.2$ false positives under the null hypothesis—as directly illustrated by the `DSPP_Amis` negative-control cells achieving $p_{\text{emp}} = 0.008–0.043$.
2. **Why `Strhab` *ENAM* orthology does not rely on `ENAM_Acar` ($p_{\text{emp}} = 0.024$):** In `Strhab`, orthology of *ΨENAM* is established independently by:
   - Archosaurian `ENAM_Amis` `tBLASTn` (`89.4 bits`, excess `57.1 bits`, $p_{\text{emp}} = 0.0010$, Table S3c);
   - Whole-proteome reciprocal `BLASTx` against 48,040 *Alligator* proteins (`XP_019346377.1`, `93.6 bits`, $E = 2.0\times 10^{-18}$, Table S2b); and
   - Splice-aware `miniprot` alignment spanning `866 aa` (`80.3%` of *Alligator* *ENAM*, `8` frameshifts + `9` PTCs, Table S4).
3. **Correction in Section 2.3 (`manuscript.typ` L207 / `manuscript.md` L86):** Replaced `"while confirming Strhab ENAM_Acar at 0.80 bits, empirical p = 0.024"` with:
   > `"while leaving Strhab ENAM_Acar at a marginal 0.80-bit excess, unadjusted exploratory empirical p = 0.024 across the 504-cell grid without multiplicity correction; orthology of Strhab ENAM is established independently by whole-proteome reciprocal BLASTx in Table S2b and 866-aa miniprot alignment in Table S4 rather than by this exploratory score"`

---

## 2. C2 — Non-Independence Among Linked Sites in *ENAM* Exon 8 (Clopper–Pearson 95% CIs)

- **Audit Finding:** All 16 shared avian *ENAM* Exon 8 PTCs reside within a single 2,634-bp coding exon on chromosome Z and share a common gene genealogy across the 7 aligned birds. Consequently, treating the 16 PTC sites as independent Bernoulli trials in an exact Clopper–Pearson binomial confidence interval (`6/16 = 37.5%`, `15.2–64.6%` on `T1`; `5/16 = 31.2%`, `11.0–58.7%` under reversible ML ASR) is valid only as a **conditional binomial description of finite-site sampling uncertainty within this single exon**, and cannot be extrapolated as a genome-wide or multi-locus standard error across independent pseudogenes.
- **Correction in Section 2.5 (#6), Section 3.4, and Section 3.5 (#2):** Explicitly qualified all Clopper–Pearson intervals as `"nominal linked-site Clopper–Pearson 95% CI"` and `"conditional descriptions of finite-site sampling within one locus rather than independent evolutionary loci"`.

---

## 3. C3 — Phylogenomic Grounding of Palaeognathae Topology `T1` vs. Single-Exon RELL Tests (`Table_S9`)

- **Audit Finding:** In `Table_S9`, `T1` (Novaeratitae = `Aptrow` + `Dronov`, `Strcam` basal) achieves a higher log-likelihood on `ENAM_exon8_aves` than `T2` (`Strcam` + `Dronov`; $\Delta\ln L = -22.96$ MAFFT / $-22.10$ MUSCLE, $p_{\text{SH}} \le 0.0040, p_{\text{AU}} \le 0.0015$) or `T3` (`Strcam` + `Aptrow`; $\Delta\ln L = -23.18$ / $-22.15$, $p_{\text{SH}} \le 0.0033, p_{\text{AU}} \le 0.0017$). However:
  1. Computing RELL topology tests on `ENAM_exon8_aves` does **not** constitute an independent phylogenomic validation of the species tree because `ENAM_exon8_aves` is the exact same alignment on which the 16 PTCs are being mapped.
  2. Genome-wide analyses of Palaeognathae (Sackton et al., 2019) demonstrate substantial incomplete lineage sorting (ILS) across the rapid basal ratite radiation, meaning that individual genomic loci genuinely have discordant gene trees (`T1`, `T2`, or `T3`) and rejection of `T2`/`T3` on one exon does not make `T2`/`T3` biologically impossible across the genome.
- **Correction in Section 2.5 (`manuscript.typ` L270 / `manuscript.md` L141):** Rewritten to state:
  > `"Within this single-locus alignment, ENAM Exon 8 exhibits a higher log-likelihood under T1 in both MAFFT (ln L = −15,237.86) and MUSCLE (ln L = −15,423.30) than under T2... or T3...; however, because these RELL tests are computed on the same pseudogene alignment used to score PTC congruence—and because incomplete lineage sorting is well documented at the base of Palaeognathae (Sackton et al., 2019)—we anchor T1 on external phylogenomic evidence and evaluate PTC incongruence across all three topologies."`
