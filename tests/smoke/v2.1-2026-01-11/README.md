## pgscalculator v2.1 local smoke test (2026-01-11)

This folder is a lightweight, reproducible smoke test for the current `pgscalculator-v2.sh` config format.

It reuses the same **genotypes** and **LD reference** paths as the existing `test-v2/test_scenarios.sh`.

### Paths you may want to override

- **Sumstat input folder**: defaults to `sumstat_5759`, override by setting `SUMSTAT_DIR`.
  - Default: `/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759`

### Run (interactive, no SLURM)

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"

bash tests/smoke/v2.1-2026-01-11/run_local.sh
```

### Run (SLURM)

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"

bash tests/smoke/v2.1-2026-01-11/submit_slurm.sh
```

### Run (SLURM, fresh outdir, all chromosomes)

Use a **new output directory** and run **all steps via --sbatch** (prep driver + prep array, then per-sumstat driver + sumstat/weights_sbayesr/weights_benchmark/score arrays + finalize job):

```bash
cd /faststorage/project/ibp_pipeline_pgscalculator/pgscalculator

export SUMSTAT_DIR="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0/sumstat_5759"

bash tests/smoke/v2.1-2026-01-11/run_sbatch_full.sh
```

- **Outdir:** `tests/smoke/v2.1-2026-01-11/out_sbatch`. The script removes this directory at start for a **from-scratch run** (prep + sumstat/weights/score/finalize all run from scratch).
- Submits prep, then the per-sumstat driver with **`--dependency=afterok:PREP_JOBID`** so the sumstat driver runs only after prep completes. Driver runs sumstat, weights, score arrays, then submits a **finalize** job (combine-scores + finalize-output) and waits for it.

### Outputs

- **Config:** `config.yaml` (single smoke config).
- **Outdir (local):** from config `outdir` (e.g. `tests/smoke/v2.1-2026-01-11/out_new`).
- **Outdir (full sbatch):** `tests/smoke/v2.1-2026-01-11/out_sbatch`.
- **Per-sumstat:** `<outdir>/sumstats/<sumstat_name>/...`
- **Driver / finalize logs:** under `<outdir>/sumstats/<sumstat_name>/logs/slurm/`

### Manual checks for finalize (sort/join)

After a successful run, you can confirm the finalize step (all joins via `LC_ALL=C sort` + `join`) by inspecting:

| File | What to check |
|------|----------------|
| `<outdir>/sumstats/<sumstat_name>/posteriors_combined.tsv` | Header has `GENO_ID`; first data column is RSID, second is GENO_ID (or NA). No in-memory join. |
| `<outdir>/sumstats/<sumstat_name>/sumstat_augmented.tsv.gz` | Columns include original sumstat + `GENO_ID`, `POST_EFFECT`, `POST_PIP`, `IN_ANALYSIS`. Row count matches formatted sumstat (join keeps all rows with NA where no match). |
| `<outdir>/sumstats/<sumstat_name>/augmented_sumstat.gz` | Header: `RSID`, `EffectAllele`, `OtherAllele`, `B`, `SE`, `Z`, `P`, `MAF`, `postEffect`, `benchEffect`. Same number of rows as variant_map; MAF column is NA when `references.maf_file` is false (as in this smoke config). |

### Work directory + cleanup

- **Per-sumstat working files** (intermediate step outputs): `out/sumstats/<sumstat_name>/work/`
- **Per-sumstat tmp**: `out/sumstats/<sumstat_name>/tmp/`
- **Prep tmp**: `out/prep/tmp/`

By default (for development), work/tmp are **kept**. If you want to remove them after a successful run, pass:
- `--cleanup` to `pgscalculator-v2.sh` (works for both local runs and SLURM driver jobs)

### Warnings you may see

If you run against an **existing outdir** where format-sumstat or filter-variants was already completed in a previous (possibly partial) run, you may see:

1. **`chrN: missing filtered sumstat: .../work/filtered/chrN_filtered.tsv (writing empty posteriors and continuing)`**  
   calc-posteriors looks for per-chromosome filtered files. If filter-variants was marked completed from an older run that only produced some chromosomes (e.g. chr21–22), or the step dir was migrated from a legacy layout, those files may be missing for other chromosomes. The step writes placeholder posteriors and continues.

2. **`Total variants mapped: 0`**  
   format-posteriors had no non-placeholder posteriors to map (follows from missing filtered inputs).

3. **`chrN: No variants in posteriors file, creating empty score`**  
   calc-score creates an empty score file when the posteriors file for that chromosome is missing or placeholder-only (again, downstream of missing filtered inputs).

**Cause:** Reusing `out/` from a partial or legacy run. Steps that were "already completed" are skipped, so the on-disk state (e.g. only chr21–22 in `work/filtered/`) is used as-is.

**Clean run (no warnings):** To regenerate all chromosomes and avoid these warnings, either:

- Use a **fresh outdir** (e.g. set a different `outdir` in the config or override with `-o`), or  
- **Remove the sumstat work dir** and re-run sumstat,weights,score:
  ```bash
  rm -rf tests/smoke/v2.1-2026-01-11/out_new/sumstats/sumstat_5759/work
  bash tests/smoke/v2.1-2026-01-11/run_local.sh
  ```
- Or run the per-sumstat steps with **`--force`** so format-sumstat and filter-variants re-run and produce all chr1–22 files (slower).
sts/smoke/v2.1-2026-01-11/run_local.sh
  ```
- Or run the per-sumstat steps with **`--force`** so format-sumstat and filter-variants re-run and produce all chr1–22 files (slower).
