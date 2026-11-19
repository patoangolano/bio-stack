# Phase B Audit Report — Methodological & Epistemological Review of the Coverage-Masked K80 Neutral Simulation

**Status:** `CORRIGIDO`  
**Audit Date:** 2026-11-19  
**Target Branch:** `audit/biorxiv-pre-submission-final`  
**Reference Table:** [`submission_biorxiv/supplementary/Table_S8_k80_simulation_sensitivity_v102.tsv`](../../submission_biorxiv/supplementary/Table_S8_k80_simulation_sensitivity_v102.tsv) (`48` parameter combinations $\times$ `1,000` replicates = `48,000` simulations, `seed 20261123`)

---

## 1. B1 — Distinguishing Model-Conditional Statistical Compatibility from Causal Proof

We audited all occurrences of causal or deterministic language surrounding the K80 Monte Carlo simulation in `manuscript.typ` and `manuscript.md`:

| Location | Prior Formulation (`v1.0.2`) | Corrected Formulation (`audit/biorxiv-pre-submission-final`) | Justification |
|---|---|---|---|
| **Abstract** | `"matching coverage-masked neutral simulations"` | `"consistent with coverage-masked K80 simulations"` | Avoids implying exact identity or unique causal proof. |
| **Introduction (#4)** | `"consistent with empirical coverage-masked K80 neutral simulations"` | `"statistically compatible with empirical coverage-masked K80 neutral simulations"` | Explicitly frames agreement as statistical compatibility. |
| **Results 2.5 (#6)** | `"Varying kappa... confirmed that incongruence rates remain high (30.8–48.0%) across all plausible neutral regimes (Table S8)."` | `"Varying kappa... shows that high PTC incongruence (30.8–48.0%) is expected across plausible neutral regimes (Table S8). Because the coverage mask is derived from the observed alignment, this agreement represents statistical compatibility conditional on the K80 model, branch lengths, topology, and empirical mask rather than an independent causal proof of purely neutral evolution."` | Replaces `"confirmed"` with conditional expectation and explicitly states the conditioning on the empirical coverage mask. |
| **Discussion 3.1 (#4)** | `"matching empirical coverage-masked neutral stop-gain expectations"` | `"statistically compatible with empirical coverage-masked K80 neutral simulations"` | Aligns summary claim with model-conditional interpretation. |
| **Discussion 3.4** | `"Combining dual-aligner consensus... demonstrates quantitatively why isolated PTCs carry a high ambiguity burden... Crucially, accounting for empirical per-species exon coverage resolves the apparent discrepancy..."` | `"Combining dual-aligner consensus... illustrates quantitatively why isolated PTCs carry a high ambiguity burden... Because the coverage mask is conditioned on the observed alignment, these simulation results establish model-conditional statistical compatibility rather than proving that memoryless K80 neutral drift is the sole evolutionary mechanism."` | Differentiates quantitative illustration and model compatibility from causal exclusion of alternative processes. |

---

## 2. B2 — Evaluation of Methodological Conditioning (Partial Circularity of Empirical Coverage Masks)

The K80 simulation in `Table_S8` conditions on three empirical inputs derived from `ENAM_exon8_aves`:
1. **The empirical per-species codon coverage mask (`cov_dual` / `cov_mafft`):** Determining which of the 878 codons are observable (`True`) vs. deleted/unaligned (`False`) in each of the 7 birds.
2. **The empirical `IQ-TREE 3` (`GTR+G`) branch lengths on topology `T1` (`T2`, `T3`):** Governing nucleotide substitution probabilities per branch.
3. **The root sequence (`Neornithes_ML` or `Allmis_ref`):** Providing the initial sense-codon composition (specifically the frequency of 1-step stop-mutable codons such as `TGG`, `CAA`, `CAG`, `CGA`, `TAT`, `TCA`, `TCG`, `GGA`).

### Why Conditioning on the Coverage Mask Is Necessary Yet Non-Independent:
- **Why it is necessary:** Our empirical pipeline can only call a shared PTC at codon $c$ if at least two birds *retain aligned nucleotides* at codon $c$ in both aligners (`cov_dual`). Comparing an empirical count observed through a `22.7–82.5%` mask against a simulation observed through a `100%` (`878/878`) mask conflates *point-mutation accumulation* with *large-scale sequence deletion/truncation*. Conditioning on `cov_dual` isolates the point-mutation process on the surviving sequence blocks.
- **Why it limits validation independence:** Because the mask `cov_dual` is taken directly from the observed alignment rather than generated de novo by an independent indel/deletion model, the simulation does not independently explain *why* `Taegut` retains only `22.7%` of Exon 8 or *why* `0/878` codons are simultaneously covered across all 7 birds. Furthermore, if local substitution rates or deletion probabilities correlate with codon composition along Exon 8, conditioning on the surviving blocks inherits that empirical pattern.
- **Effective influence on conclusions:** The 48-scenario grid in `Table_S8` already includes both the unmasked (`unmasked_full_878`) and masked (`empirical_mafft_mask`, `empirical_dual_aligner_mask`) regimes across two independent root sequences (`Allmis_ref` and `Neornithes_ML`), three $\kappa$ values (`2.0, 4.0, 6.0`), four branch-length scales (`0.65×` to `1.25×`), and three topologies (`T1`, `T2`, `T3`). Across **every** combination—including the unmasked regime—the neutral K80 model produces high parsimony incongruence (`30.8–53.6%`), showing that high PTC homoplasy/reversal is an intrinsic consequence of long avian branch lengths rather than an artifact of the coverage mask, whereas the coverage mask specifically governs the *absolute count* of observable shared PTCs (`25.0–29.0` unmasked vs. `16.0–19.0` masked, matching the `16` observed).

---

## 3. B3 — Topological & Parameter Nuance Preserved from `Table_S8`

The 48,000 simulations in `Table_S8` (verified intact by SHA-256 `4939821f6814e095454ca4de03c2a1b341391340f2a5710cc9a09f631523f08c`) were preserved without re-execution. In Section 2.5 (#6), Section 3.4, and Section 3.5 (#2), we explicitly note that:
1. On the primary topology **`T1`**, observed parsimony incongruence (`6/16 = 37.5%`) lies near the center of the coverage-masked K80 distribution (unconditional median `40.0–44.4%`, two-sided $p = 0.618–0.848$; conditional median at $N_{\text{shared}}=16$ of `37.5%`).
2. On the alternative topologies **`T2` (`10/16 = 62.5%`)** and **`T3` (`11/16 = 68.8%`)**, the observed parsimony incongruence exceeds the simulated neutral medians (`43.8–48.0%` on `T2`, two-sided $p = 0.162–0.244$; `41.7–47.6%` on `T3`, two-sided $p = 0.072–0.092$), sitting toward the upper tail of the neutral K80 distribution because forcing the five shared Novaeratitae (`Aptrow` + `Dronov`) PTCs onto non-monophyletic topologies requires five additional parallel gains or reversals.
