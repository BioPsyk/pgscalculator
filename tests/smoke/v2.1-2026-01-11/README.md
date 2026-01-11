## pgscalculator v2.1 local smoke test (2026-01-11)

This folder is a lightweight, reproducible smoke test for the current `pgscalculator-v2.sh` config format.

It reuses the same **genotypes** and **LD reference** paths as the existing `test-v2/test_scenarios.sh`.

### Paths you may want to override

- **Sumstat input folder**: defaults to `sumstat_5759`, override by setting `SUMSTAT_DIR`.
  - Default: `/faststorage/project/ibp_data_open/cleansumstats/version_1.12.0/sumstat_5759`

### Run (interactive, no SLURM)

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_data_open/cleansumstats/version_1.12.0/sumstat_5759"

bash test-v2-2026-01-11/run_local.sh
```

### Run (SLURM)

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_data_open/cleansumstats/version_1.12.0/sumstat_5759"

bash test-v2-2026-01-11/submit_slurm.sh
```

### Outputs

- **Outdir**: `test-v2-2026-01-11/out/`
- **Per-sumstat**: `out/sumstats/<sumstat_name>/...`
- **Driver log**: printed on submission, and written under `out/sumstats/<sumstat_name>/logs/slurm/`

