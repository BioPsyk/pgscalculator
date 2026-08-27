# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Per-method augmented sumstats:** `augmented_sbayesr.gz` and `augmented_ldpred2.gz`, each restricted to its own LD-reference variant set, keyed by rsid, with its own `benchEffect` column kept inside the file.
- **Per-method benchmarks:** `calc-benchmark` runs once per active method (`work/benchmark_<method>/`), and finalize emits `bench_score_sbayesr.gz` / `bench_score_ldpred2.gz`.
- **LDpred2 diagnostics:** `details/ldpred2/summary.tsv` (mode, seed, h2/p/alpha estimates, LDSC intercept/h2, match/QC/chain counts) and `details/ldpred2/chains.png` (overlaid p/h2 sampling paths of all kept LDpred2-auto chains).

### Removed

- **Legacy output symlinks:** no longer create `scores.gz`, `bench_score.gz`, or `augmented_sumstat.gz` aliases pointing at the sBayesR per-method files. Use `scores_sbayesr.gz`, `bench_score_sbayesr.gz`, and `augmented_sbayesr.gz` directly.
- **Legacy v1 Nextflow pipeline.** Removed `main.nf`, `nextflow.config`, `modules/`, `conf/prscs.config`, `conf/sbayesr.config`, `lib/`, the v1 `pgscalculator.sh` wrapper, `concatenate_plink_maf/` (+ `tests/unit/test_concatenate_plink_maf.sh`), and `scripts/kill-nextflow.sh`. The full v1 pipeline (1.3.2) is archived at git tag **`v1.3.2`** — run `git checkout v1.3.2` to recover it. `conf/init-docker-config.sh` and `assets/` are retained (still used by v2).
- `docker/Dockerfile` no longer installs Nextflow (removed the `java_builder` stage and `.nextflow`/`NXF_OFFLINE` setup); the runtime user was renamed `nextflow` → `pgsuser`. Base image bumped to **`0.8.0`** (`docker/VERSION`) and the `README.md` pull/build references updated to match; `0.8.0` must be built + pushed to Docker Hub. The previous `0.7.0` image still runs v2 unchanged (same toolchain).

### Changed

- **v2 is now the default pipeline.** `README-v2.md` → `README.md`, `pgscalculator-v2.sh` → `pgscalculator.sh`, and `VERSION.v2` → `VERSION` (single root version file, `2.2.0`). All in-repo references (docs, tests, `bin/lib/steps/*`, the SLURM driver's self-invocation, `Dockerfile.deploy`, `scripts/init-containerization.sh`) were updated to the new names.
- `finalize-output` produces one self-contained augmented file per discovered method instead of a single combined `augmented_sumstat.gz` with `postEffect_<method>` columns. This avoids result files containing variants outside a method's LD reference.
- LDpred2 chain diagnostics moved from `logs/ldpred2_chains.png` to `details/ldpred2/chains.png` (now overlaying all kept chains rather than the first only).
- **Interactive multi-method runs:** the `sumstat`/`weights`/`score` step groups are now method-aware, so a plain `run --all` (or `run --steps sumstat,weights,score`) runs the full sBayesR *and* LDpred2 paths locally (one pass), matching the SLURM driver. Previously interactive `run --all` was sBayesR-only and the LDpred2 path required the driver.
- **Parallel LDpred2 inclusion list:** `prep-inclusion-list-ldpred2` now builds its per-chromosome variant-map partials independently (a new `prep-inclusion-list-ldpred2-combine` step concatenates them), mirroring `prep-inclusion-list`. The SLURM driver builds the LDpred2 partial in the same per-chromosome prep array task as the sBayesR one, so genome-wide prep no longer runs the LDpred2 map join as a single ~18-min sequential pass.
- **LDpred2 sample size:** `run_ldpred2.R` now prefers the per-variant `N` column (the LDpred2-recommended input) and uses the metadata-derived scalar (`--effective-sample-size` / case-control) only to fill rows with missing `N`. Previously the metadata scalar overrode the entire per-variant column. The N source split is recorded in `details/ldpred2/summary.tsv` (`n_per_variant`, `n_metadata_fallback`).

### Fixed

- `calc-ldpred2` no longer aborts under `set -u` when `input:` (cleansumstats metadata) is not configured.
- **Incremental finalize:** `finalize-output` is now method-aware. Previously a single `details/.completed` marker meant a Day-2 run (adding a new method) was skipped entirely, so the new method's `augmented_<method>.gz` was never written. Finalize now re-runs whenever a discovered method is missing its augmented file.
- **Single-chromosome local prep:** `prep-inclusion-list` / `prep-inclusion-list-ldpred2` no longer fall into per-chromosome "array task" (partial-only) mode just because the config lists a single chromosome (e.g. `chromosomes: 22`). The wrapper now sets `single_chr_task` only when `--_chr` is passed (driver array task), so a normal local `--steps prep` on a single-chromosome config runs the full prep including the combine.
- **Silent prerequisite failures:** the `pgscalculator.sh` prep/sumstat/posteriors prerequisite checks captured the helper's status with `var=$(...); if [[ $? -ne 0 ]]`, which under an inherited `errexit` aborted at the assignment before printing guidance. They now use `if ! var=$(...)`, so the "Run … first" instructions are actually shown.

### Docs

- `README-v2.md`: "What's New in v2.2", pipeline-step/step-group tables, and the output tree updated to the per-method augmented/benchmark model (`augmented_<method>.gz`, `bench_score_<method>.gz`, `work/benchmark_<method>/`, `details/ldpred2/`). The **Incremental runs** section now shows the correct Day-2 command (`--steps sumstat,weights,score,finalize`, since `filter-variants` must run for the new method), documents the byte-identity guarantee for Day-1 sBayesR outputs, and references the chr22 acceptance test.
- `docs/plans/general-design.md`: "Posterior methods (v2.2)", incremental/sequential-runs, per-method failure isolation, the output directory tree, and the augmented/benchmark output schemas updated from the single-combined-file model to per-method files. The historical finalize-join recipe now carries a v2.2 supersession note. `config.template.yaml` confirmed already current (`methods`, `whichn`, `sbayesr.ld_build`, full `ldpred2` block, per-method SLURM profiles).

### Tests

- Unit: `test_variant_map_for_ldpred2.sh` (LDpred2 variant-map join: direct/swap/complement/miss, af→a2freq reconciliation, array-task partial vs. full-run+combine equivalence, opt-in skip); method-aware finalize completion guard added to `test_incremental_finalize.sh`; prep-step classifier regression guard added to `test_driver_dispatch.sh`.
- Fixed missing executable bits on `test_driver_dispatch.sh`, `test_format_posteriors_ldpred2.sh`, `test_phase9_discovery.sh` (the runner stopped before reaching them, so they never ran).
- Integration: `tests/smoke/v2.2-2026-05-26/integration_incremental_chr22.sh` — real-data chr22 incremental acceptance (Day-1 sBayesR → Day-2 add LDpred2; asserts `scores_sbayesr.gz`/`augmented_sbayesr.gz` byte-identical and both methods' outputs + diagnostics present).

## [2.2.0] - 2026-05-26

pgscalculator **v2** release (wrapper `pgscalculator.sh`, CLI `bin/pgscalculator`). The v1
Nextflow pipeline version remains in root `VERSION` (unchanged).

### Added

- **LDpred2** posterior method (`calc-ldpred2`, `run_ldpred2.R`) via `bigsnpr` / LDpred2-auto
- **`methods:`** config key and **`--methods`** CLI (`sbayesr`, `ldpred2`, or both)
- **Prep (opt-in):** `prep-ldref-ldpred2`, `prep-inclusion-list-ldpred2` when `ldpred2.ld_dir` is set
- **Per-method work dirs:** `filtered_<method>/`, `posteriors_mapped_<method>/`, `scores_<method>/`
- **Per-method outputs:** `scores_sbayesr.gz`, `scores_ldpred2.gz`; `scores.gz` → `scores_sbayesr.gz` when sBayesR-only
- **Discovery-mode finalize:** `augmented_sumstat.gz` with `postEffect_<method>` columns only when mapped posteriors exist on disk
- **Incremental runs:** run sBayesR and LDpred2 on different days in the same `outdir` without re-running sBayesR
- **SLURM:** `weights_ldpred2` (single genome-wide job), `score_sbayesr` / `score_ldpred2` profiles
- **Container:** multi-stage `r_builder` with `pak` + LDpred2 R stack (requires image **0.7.0+**, see `docker/VERSION`; independent of pipeline `VERSION`)
- **Tests:** unit tests for methods config, driver dispatch, format-posteriors LDpred2, incremental finalize; smoke `tests/smoke/v2.2-2026-05-26/`

### Changed

- `filter-variants` and `format-posteriors` are method-parameterised (`--method sbayesr|ldpred2`)
- `combine-scores` emits one gzipped score file per discovered method
- Legacy `work/filtered/` and `work/posteriors_mapped/` layouts are migrated automatically

### Container (0.7.0 — only if you need LDpred2 or have not rebuilt since Phase 1)

The **pipeline** version is `2.2.0` (`VERSION`). The **image** version is `0.7.0` (`docker/VERSION`). It was bumped once when the `r_builder` stage and `bigsnpr` stack landed (commit `2d6299c`); phases 2–11 did not change the Dockerfile.

Rebuild only when `docker/Dockerfile` / `docker/VERSION` change, not on every pipeline release:

```bash
./scripts/docker-build.sh
singularity build sif/ibp-pgscalculator-base_version-0.7.0.sif docker-daemon://ibp-pgscalculator-base:0.7.0
# or: singularity pull sif/ibp-pgscalculator-base_version-0.7.0.sif docker://biopsyk/ibp-pgscalculator:0.7.0-amd64
./scripts/test-r-packages.sh --singularity sif/ibp-pgscalculator-base_version-0.7.0.sif
```

sBayesR-only runs can keep using an older image (e.g. `0.6.0`) until you enable LDpred2.

## [1.3.2] - 2025-09-22
### Fixed
- **Critical sorting issue in variant_map_for_sbayesr.sh that caused incomplete variant mappings**
- Missing sort operation for snp2 file before first join operation, which led to join failures
- Incomplete or incorrect results when processing variant mapping for SBayesR workflow

### Added
- Comprehensive sorting validation tests to detect and prevent sorting issues
- Dedicated test suite for variant mapping sorting validation (`test_variant_map_sorting_validation.sh`)
- Helper functions for sorting validation in unit tests (`_check_file_sorted`, `_validate_join_prerequisites`)
- Detailed technical documentation of sorting issues and fixes
- Performance testing for large datasets (100+ variants) in sorting validation

### Changed
- Enhanced existing variant mapping tests with sorting validation capabilities
- Improved test coverage for edge cases (empty files, single records, different genome builds)

## [1.3.1] - 2025-09-13
### Added
- Robust column-name-based MAF extraction script for better maintainability
- PSAM file standardization to IID-only format to ensure consistent score calculations
- Comprehensive unit tests for MAF extraction functionality

### Fixed
- Zero scores in combined output due to FID/IID format inconsistencies between genotype files
- Incorrect MAF values and missing NCHROBS in raw_maf_chrall output
- Pipeline sensitivity to varying plink2 .afreq output column formats
- Variant ID mangling when using genotype files with missing rsIDs

### Changed
- MAF extraction now uses column names instead of positional indices for robustness
- Score calculation simplified after implementing consistent PSAM standardization

## [1.3.0] - 2025-06-29
### Added
- **PLINK2 dosage format support as default genotype input format**
- **Automatic PLINK1 to PLINK2 conversion with format detection**
- Multi-architecture support for Docker images (amd64 and arm64)
- Docker manifest support for seamless cross-platform deployment
- Enhanced genotype file processing with improved variant ID handling
- Comprehensive unit tests for genotype format conversion and processing

### Changed
- **Default genotype format changed from PLINK1 to PLINK2 dosage format**
- Updated Docker build and push scripts to handle multi-arch builds
- Improved Docker image distribution with platform-specific tags
- Enhanced variant mapping and benchmark scoring for PLINK2 format
- Improved debugging and channel tracing capabilities

## [1.2.7] - 2025-05-19
### Fixed
- Intermediate scores to be in method specific folders

## [1.2.6] - 2024-11-29
### Added
- Enhanced FAQ documentation with detailed explanation of benchmark calculations

## [1.2.5] - 2024-11-29
### Fixed

- renaming a duplicate header in augmentet output

## [1.2.4] - 2024-09-02
### Changed

- remove memory upper limits in nextflow.config. Replaced by setting plink and sort memory variables for each specific process in nextflow.config

### Fixed
- Remove the tmp in mount init script

## [1.2.3] - 2024-08-23
### Changed

- Removed memory restriction set in nextflow.config, and removed all labels from all processes

## [1.2.2] - 2024-08-16
### Changed

- README.md to keep only the absolutely most important to run using singularity
- Extra documentation moved to their own doc files in docs/

### Fixed
- Missing mount dir for supplied snplists in pgscalculator.sh

## [1.2.1] - 2024-04-30
### Added

- Memory constraints on plink2, so that the overall memory footprint will be lower.
- Checker for missing or duplicate variant ids in the genotype data, which will fill in with chr:pos:a1:a2 if missing, and add an extra number if still a duplicate.

## [1.2.0] - 2024-04-24
### Added

- b38 support and improved allele flip management. Also added an early filter on b37 NA coordinates.

## [1.1.0] - 2024-04-10
### Added

- The sbayesr workflow has been improved in all aspects: a companion variant map file, an augmentated sumstat file, and in general much better channeling and output structure.

## [1.0.0] - 2023-11-10
### Added

- Everything to run an sbayesr workflow.

## [0.1.0] - 2023-10-02
### Added

- Basic containerized setup, and a general reformatter for the different prs softwares

