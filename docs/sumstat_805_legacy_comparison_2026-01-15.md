# `sumstat_805` legacy vs v2 (v2.1) comparison (2026-01-15)

This note captures the comparison between the legacy pgscalculator outputs from the parameter-search runs and the v2.1 pipeline outputs for `sumstat_805`.

## Runs compared

- **Legacy output (reference)**
  - Path: `/home/jesgaaopen/ibp_pipeline_pgscalculator/test-zone/out_parameter_search/sbayesr_default.config/sumstat_805`
  - Files used:
    - `calc_posteriors/allchr.posteriors` (combined posteriors)
  - Notes:
    - No legacy `scores` or `augmented sumstat` artifacts were present in this run output.

- **v2.1 output (best-parity config)**
  - Path: `/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/tests/smoke/v2.1-2026-01-11/out_paramsearch_805_default/sumstats/sumstat_805`
  - Config used: `pgscalculator/tests/smoke/v2.1-2026-01-11/config.sumstat_805.paramsearch_default.yaml`
  - Files used:
    - `posteriors_combined.tsv`
    - `scores.tsv.gz`
    - `sumstat_augmented.tsv.gz`
    - `variant_map.tsv.gz`

## Posterior overlap and correlations

Comparison performed by **RSID + A1 + A2** (allele-aware, no swaps required in the matched set):

- Legacy posteriors: `610,107` variants (`allchr.posteriors`)
- v2 posteriors: `367,330` variants (`posteriors_combined.tsv`)
- Overlap (RSID + A1 + A2): `307,284` variants

Correlations on the overlap:

- **Posterior effect correlation**: `0.9500`
- **PIP correlation**: `0.4758`

## Notes / interpretation

- The overlap set is smaller than either full posterior set, indicating substantial filtering or input differences between the legacy parameter-search run and the v2.1 run.
- Despite the smaller overlap, **effect sizes agree strongly** (r ~ 0.95), while **PIP agreement remains moderate** (r ~ 0.48).
- Since the legacy parameter-search run does not include scores or augmented sumstats, **score-level comparisons are not available** for this sumstat.

## Repro notes

- Comparison computed from:
  - Legacy: `/home/jesgaaopen/ibp_pipeline_pgscalculator/test-zone/out_parameter_search/sbayesr_default.config/sumstat_805/calc_posteriors/allchr.posteriors`
  - v2: `/faststorage/project/ibp_pipeline_pgscalculator/pgscalculator/tests/smoke/v2.1-2026-01-11/out_paramsearch_805_default/sumstats/sumstat_805/posteriors_combined.tsv`
