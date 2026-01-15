# `sumstat_814` legacy vs v2 (v2.1) comparison (2026-01-14)

This note captures the **end-to-end findings** from comparing the legacy pgscalculator pipeline outputs for `sumstat_814` against the v2.1 pipeline outputs, including the main root causes behind early “bad comparisons”, the fixes applied, and what remains unresolved.

## Runs compared

- **Legacy output (reference)**
  - Path: `/home/jesgaaopen/ibp_pipeline_pgscalculator/test-zone/out_test3_C/sumstat_814`
  - Files:
    - `main_raw_score_all.gz`
    - `variant_map.gz`
    - `augmented_sumstat.gz`
    - `extra/raw_posteriors_chrall` (combined posteriors)
    - Legacy config: `details/sbayesr.config` (notable flags include `impute_n: true`, and legacy had `rsq: 0.95` + `p_value: 0.99`)

- **v2.1 output (latest “no-impute” rerun)**
  - Path: `/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/tests/smoke/v2.1-2026-01-11/out_full814/sumstats/sumstat_814`
  - Config used: `pgscalculator/tests/smoke/v2.1-2026-01-11/config.full.noimpute.yaml`
    - Key switch: `sbayesr.impute_n: false`
    - Full chromosomes: `chromosomes: 1-22`
  - Files:
    - `scores.tsv.gz` (v2 score output)
    - `variant_map.tsv.gz` (v1-compatible artifact written by `finalize-output`)
    - `sumstat_augmented.tsv.gz`
    - `posteriors_combined.tsv`

## Findings timeline (what went wrong, then what got fixed)

### 1) Early “bad correlation” was largely due to comparing different chromosome sets

Initial comparisons were misleading because **v2 was effectively running only chr21/22** while the intent was chr1-22:

- Cause: The v2 run reused a `prep/` directory that was generated for a **restricted chromosome set**.
- Symptom: `calc-posteriors` reported issues for many chromosomes; only chr21/22 produced meaningful outputs.
- Fix: The v2 wrapper now validates that chromosome-specific prep outputs exist for *all configured chromosomes* and fails early if not.
  - Wrapper: `pgscalculator/pgscalculator-v2.sh`

This fix ensured subsequent “full-genome” comparisons were actually chr1-22 vs chr1-22.

### 2) A `posteriors_combined.tsv` column mismatch corrupted downstream augmented sumstat

After getting a full-genome run, we hit a second major correctness issue:

- Symptom: `posteriors_combined.tsv` header contained `GENO_ID` but the data rows did not, causing downstream parsers to crash and/or map the wrong columns into `sumstat_augmented.tsv.gz`.
- Fix: `finalize_output.sh` now injects `GENO_ID` from `rsid_to_genoid.tsv` when combining posteriors, and `generate_augmented_sumstat()` uses the correct column indices afterwards.

### 3) `filter-variants` could OOM due to non-streaming EAF mapping

We hit OOM (and “0 variants” failures) in the v2 `filter-variants` step:

- Root cause: `force_eaf()` could load a large LD-reference EAF table into `awk` associative arrays, which is memory-heavy.
- Fix:
  - Add a fast-path: if the input sumstat already has a complete EAF column, skip LD-ref EAF mapping entirely.
  - Reduce memory footprint by only building allele alignment maps when needed.
  - Bump driver memory to `2g` in the smoke full-chr config so driver-side parsing/fallback work doesn’t get killed.

### 4) sbayesR option mismatches materially changed the “in-analysis” set

We observed a very low overlap (and weak correlations) between legacy and v2 “in-analysis” variants and outputs. A key contributor:

- v2 config originally included sbayesR filters `--rsq` and `--p-value` that were not present (or differed) in the legacy workflow in the way we were expecting.
- After removing those from the v2 smoke config, score correlation improved (earlier improvement observed from ~0.11 to ~0.50 in an intermediate rerun).

### 5) `.ma` header consistency (GCTB sensitivity)

We also ensured the `.ma` file passed to GCTB has a consistent header (`SNP A1 A2 freq b se p N`) across chromosomes to avoid tool sensitivity to header casing/format.

### 6) “Aha”: `--impute-n` triggers GCTB’s per-SNP N outlier filter and drops many SNPs

This turned out to be an independent, high-impact behavior difference (details below).

## Important pipeline fixes made earlier (code-level summary)

These fixes were needed before comparisons were meaningful:

- **Avoided accidental reuse of incomplete `prep/`**
  - The v2 wrapper was hardened to verify *all configured chromosomes* have chromosome-specific prep outputs, so we don’t silently reuse a chr21/22-only prep in a chr1-22 run.

- **Fixed `posteriors_combined.tsv` schema**
  - `finalize_output.sh` was patched to actually populate `GENO_ID` in `posteriors_combined.tsv`, and to correctly map `POST_EFFECT` / `POST_PIP` columns into `sumstat_augmented.tsv.gz`.

- **Reduced `filter-variants` OOM risk**
  - `force_eaf()` was optimized to avoid loading LD-ref EAF into memory when EAF is already complete in the input sumstat (fast-path).
  - Smoke config driver memory was bumped to `2g` to cover cases where mapping is still required.

## Key finding: `--impute-n` triggers a large SNP drop inside GCTB

On a chromosome-level A/B test (chr1), GCTB behaves very differently depending on whether `--impute-n` is enabled:

- With **`impute_n: true`**, GCTB reported:
  - `47446 matched SNPs in the GWAS summary data`
  - `33443 SNPs with per-SNP sample size within 3 sd around the median ...`
  - `33443 SNPs on 1 chromosomes are included.`

- With **`impute_n: false`**, GCTB reported:
  - `47446 matched SNPs in the GWAS summary data`
  - `47446 SNPs on 1 chromosomes are included.`

Interpretation:

- Even if the `.ma` file provides an `N` column, enabling `--impute-n` can lead GCTB to compute an **imputed per-SNP N**, and then apply its internal **per-SNP N outlier filter** (“within 3 sd around median”).
- This can drop a substantial fraction of SNPs and changes both **posterior set size** and downstream **scores**.

This was the motivation for the full chr1-22 rerun using `impute_n: false`.

## Latest v2.1 no-impute full run summary

- **Combined posteriors**: `561,851` variants (`posteriors_combined.tsv`)
- **Scores**:
  - `scores.tsv.gz` reports `N_VARIANTS=561851` per sample (consistent with the posteriors combined size)

## Old vs new: scores correlation (sample-wise)

Comparison performed on all 2504 samples:

- Legacy: `main_raw_score_all.gz` using `SCORE1_SUM`
- v2: `scores.tsv.gz` using `SCORE_SUM`

Result:

- **Correlation**: ~`0.4163`
- **Delta range** (`new - old`): approx `[-0.274, +0.315]`

Notes:

- This is an improvement versus earlier runs where posteriors were computed from an incomplete prep or with mismatched sbayesR filtering, but it’s still far from “nearly identical.”

## Old vs new: posterior overlap and correlations

Comparison performed by RSID, allele-aware (no swaps observed in the matched set):

- Legacy posteriors (`extra/raw_posteriors_chrall`): `617,363` variants
- v2 posteriors (`posteriors_combined.tsv`): `561,851` variants

Overlap:

- **Intersection**: `309,373` RSIDs
- **Legacy-only**: `617,363 - 309,373 = 307,990`
- **v2-only**: `561,851 - 309,373 = 252,478`

Correlations (on the intersection):

- **Posterior effect correlation**: ~`0.2573`
- **PIP correlation**: ~`0.0678`

Interpretation:

- Disabling `impute_n` increases the included SNPs (relative to the “impute-n + N-outlier filter” behavior), but there remains a large mismatch in:
  - **which variants appear in posteriors**, and
  - **PIP estimates** even when RSID matches.

## Current open questions / next steps

The main remaining gap is explaining the low posterior/PIP agreement and the modest score correlation. High priority next checks:

- **Augmented sumstat diffs (old vs new)**
  - For matched variants: allele harmonization checks, distribution of `B/SE/Z/P/EAF`, and `IN_ANALYSIS` overlap.
  - Identify top outliers contributing most to score differences.

- **GCTB configuration parity**
  - Confirm whether legacy and v2 runs use the *same* inputs beyond flags (LD reference, genotype panel, inclusion list).
  - Validate whether legacy `impute_n: true` was actually used in the legacy run’s GCTB invocation (config says yes), and decide whether the correct “legacy-like” target is:
    - keep `impute_n: true` (accept the N-outlier filtering), or
    - treat it as a legacy artifact and standardize on `impute_n: false` for determinism.

- **Variant-map and mapping differences**
  - Compare `variant_map.gz` (legacy) vs `variant_map.tsv.gz` (v2) coverage and NA rates.
  - Verify v2’s `GENO_ID` mapping consistency for downstream scoring.

## Repro notes

- The latest v2 run was executed via the v2 wrapper (`pgscalculator/pgscalculator-v2.sh`) with SLURM driver submission and:
  - `pgscalculator/tests/smoke/v2.1-2026-01-11/config.full.noimpute.yaml`
  - Input sumstat dir: `/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_814`
  - Output dir: `/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/tests/smoke/v2.1-2026-01-11/out_full814`

