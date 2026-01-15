% Variant Mapping Plan (chr/pos + alleles)

## Purpose
Define the required behavior for variant ID mapping so that sumstat, LD reference, and genotype
identifiers are linked unambiguously and conversions never rely on RSIDs alone.

## Canonical mapfile schema
Single shared `chr` and `pos`, with per-source SNP IDs and allele columns:
- `chr`, `pos`
- `sumstat_snpid`, `sumstat_effect`, `sumstat_other`
- `geno_snpid`, `geno_a1`, `geno_a2`
- `ldref_snpid`, `ldref_a1`, `ldref_a2`

Notes:
- No compound `chrpos` field is required.
- No sorting requirement is imposed for runtime use.
- The mapfile may include non-key frequency columns used for EAF filling/auditing.

## Filtering order (planned)
1) **prep** builds the base `variant_map.tsv` from the **intersection** of genotype + LD reference variants
   (by `chr/pos + alleles`), including LD reference EAF (`ldref_a2freq`).
   - Also write chromosome-specific mapfiles (e.g., `prep/variant_map/chrN.tsv`) for parallel use.
2) **sumstat** attaches `sumstat_*` columns via `chr/pos + alleles` into **chromosome-specific**
   sumstat mapfiles (using the prep `prep/variant_map/chrN.tsv` files).
3) **sumstat** concatenates the per-chromosome sumstat mapfiles into a single sumstat-specific
   mapfile for output/auditing.
3) **sumstat** reduces to the **sumstat intersection of the mapfile** to produce the posterior-input
   sumstat (same row count as the matched subset, not the full union), and fills missing `EAF`
   from `ldref_a2freq` during this reduction (allele-aware, no `EAF_1KG`).
4) **sumstat** applies user-provided inclusion lists (if any), sequentially in this order:
   genotype list -> sumstat list -> ldref list.
5) **filter-variants** runs the following sub-steps:
   - Ensure N is present (from sumstat or metadata fallback).
   - Filter bad values (pass1) before derivations:
     - Missing/NA or zero `B` (beta).
     - Missing/NA or zero `SE`.
     - Missing/NA `EAF` or boundary `EAF` (0 or 1).
   - Derive B/SE as needed.
   - Filter bad values (pass2) after derivations:
     - Missing/NA or zero `B` (beta).
     - Missing/NA or zero `SE`.

## Build workflow
### Prep step (sumstat-agnostic)
1) Build `variant_map.tsv` from the **intersection** of genotype + LD reference variants.
2) Populate geno/ldref columns and add LD reference EAF (`ldref_a2freq`) to the mapfile.
3) No sumstat columns are added at prep.
4) Always derive SNP inclusion lists from the mapfile:
   - Use the prep mapfile when the prep step needs an inclusion list.
   - Use the sumstat-annotated mapfile when the sumstat step needs an inclusion list.
   - It is acceptable to derive a union inclusion list that contains all three sources, which
     can be useful for building the posterior-calculation input.

### Format-sumstat step (GRCh37 coordinate mapping)

1. Paste `cleaned_GRCh37.gz` and `cleaned_GRCh38.gz` side-by-side.
2. In a single awk pass:
   - Filter out rows where CHR or POS (b37) are empty/NA (liftover failed).
   - Split output by chromosome into per-chromosome files (`formatted/chrN.tsv`).
3. Output: Per-chromosome sumstat files ready for parallel processing.

### Sumstat step (sumstat-specific)
All steps after GRCh37 coordinate mapping should support **per-chromosome parallelism**.
Each chromosome process should only load its corresponding `prep/variant_map/chrN.tsv`.

1) Create a sumstat-specific copy of the mapfile.
2) Attach `sumstat_*` columns using `chr/pos + alleles` to match against the map.
3) Fill missing sumstat `EAF` from mapfile `ldref_a2freq` when needed (allele-aware).
4) The mapfile number of rows remains; same as in prep mapfile ; no reduction to the sumstat intersection.
5) Apply a **user-provided inclusion list** (replacing INFO/MAF filtering):
   - Users may provide **three separate inclusion lists**, one for each ID space:
     `gt` (genotype), `ss` (sumstat), and `ld` (ldref).
   - If multiple lists are provided, apply them sequentially in this order: genotype -> sumstat -> ldref.
   - Filter the sumstat **after** the mapfile reduction (i.e., after reducing to the sumstat intersection of the mapfile).
   - The filter should be applied to the **posterior input sumstat file**.

## Runtime usage (ID conversions)
All downstream conversions use the mapfile directly (no chr/pos matching at this stage):
- Before posterior calc: `sumstat_snpid` -> `ldref_snpid`
- After posterior calc: `ldref_snpid` -> `geno_snpid`

## EAF handling (current behavior)
- LD reference EAF is extracted in `prep-ldref` into `prep/references/ldref_eaf.tsv`
  (`RSID`, `A1`, `A2`, `A2Freq`) from the `.info` files.
- In `filter-variants`, `force_eaf` fills missing `EAF` in the sumstat:
  - Priority: sumstat `EAF` (if present) -> `ldref_eaf` (preferred). `EAF_1KG` is
    explicitly **not** used.
  - Mapping to `ldref_eaf` is by **RSID**, with allele alignment:
    - If sumstat effect allele matches LD `A2`, use `A2Freq`.
    - If it matches LD `A1`, use `1 - A2Freq`.
    - Otherwise, use `A2Freq` as-is.
- Filtering uses `EAF` in `filter_bad_values` to drop missing values and boundaries (0/1),
  and `EAF` is also used for deriving `B`/`SE` when needed.

## EAF handling (mapfile fit)
Source ideas consistent with the mapfile plan:
- Keep the authoritative EAF values in the sumstat file, but use the mapfile to route
  `sumstat_snpid` -> `ldref_snpid` for any fallback fill from LD reference.
- Store LD reference EAF (`ldref_a2freq`) in the **prep mapfile** so the sumstat step can
  fill `EAF` while it reduces to mapfile variants. This avoids any dependence on `EAF_1KG`.
- Keep LD_EAF in output mapfile.

## Final output
Always write the sumstat-annotated mapfile as `variant_map.tsv.gz` for auditing and back-tracing.

## SLURM submission system (`--sbatch`)

### Overview
The `--sbatch` flag provides integrated SLURM job submission. It handles two distinct
workflows: **prep jobs** and **per-sumstat driver jobs**.

### Submission modes

1. **Prep job** (`--sbatch --steps prep`):
   - Submits a single SLURM job for the prep step.
   - Resources configured via `slurm.prep: { mem, cpus, time }`.
   - Logs written to `<outdir>/prep/logs/slurm/`.

2. **Per-sumstat driver job** (`--sbatch --steps sumstat,posteriors,score -i <sumstat>`):
   - Submits a lightweight **driver job** that orchestrates the full per-sumstat workflow.
   - The driver job runs `sumstat` directly (not chromosome-parallel).
   - For `posteriors` and `score` steps, the driver submits **SLURM array jobs**
     (one task per chromosome) and waits for completion before proceeding.
   - Resources for the driver configured via `slurm.driver: { mem, cpus, time }`.
   - Logs written to `<outdir>/sumstats/<sumstat_name>/logs/slurm/`.

### Driver job architecture

```
User submission                     SLURM cluster
     │
     ▼
./pgscalculator-v2.sh --sbatch --steps sumstat,posteriors,score -i sumstat_814
     │
     └──► Driver job (pgs_sumstat_814_driver) ──────────────────────────────────►
              │
              ├── [1] Runs format-sumstat directly (GRCh37 mapping + chr split)
              │
              ├── [2] Submits sumstat array (1-22%max_parallel)
              │       └── Per-chr: mapfile join, EAF fill, filter, reduce
              │       └── Waits for all tasks to complete
              │
              ├── [3] Submits posteriors array (1-22%max_parallel)
              │       └── Waits for all tasks to complete
              │
              ├── [4] Submits score array (1-N%max_parallel)
              │       └── Waits for all tasks to complete
              │
              └── [5] Runs combine-scores, finalize-output directly
```

### Array job details
- **One task per chromosome**: Each array task runs a single chromosome.
- **max_parallel**: Configurable via `slurm.<step>.max_parallel` (default: 22).
  Set to 1 to serialize chromosome processing.
- **Job array syntax**: Uses SLURM `--array=1-N%max_parallel`.
- **Chromosome file**: The driver writes a temporary `<job_name>.chromosomes.txt` file
  and each task reads its chromosome from line `$SLURM_ARRAY_TASK_ID`.
- **Score step special handling**: The array runs only `calc-score`; the driver runs
  `combine-scores` and `finalize-output` after the array completes.

### Config structure

```yaml
slurm:
  account: ibp_pipeline_pgscalculator
  partition: normal

  # Lightweight driver job (submits arrays for sumstat/posteriors/score)
  driver: { mem: 1g, cpus: 1, time: '02:00:00' }

  # Prep step (sumstat-agnostic, single job)
  prep: { mem: 16g, cpus: 4, time: '02:00:00' }

  # Sumstat array (one task per chromosome, after format-sumstat splits)
  sumstat: { mem: 1g, cpus: 1, time: '00:30:00', max_parallel: 22 }

  # Posteriors array (one task per chromosome)
  posteriors: { mem: 20g, cpus: 6, time: '04:00:00', max_parallel: 22 }

  # Score array (one task per chromosome with mapped posteriors)
  score: { mem: 10g, cpus: 4, time: '01:00:00', max_parallel: 22 }
```

### Log locations
- **Prep logs**: `<outdir>/prep/logs/slurm/pgs_prep_<jobid>.out`
- **Driver logs**: `<outdir>/sumstats/<name>/logs/slurm/pgs_<name>_driver_<jobid>.out`
- **Array logs**: `<outdir>/sumstats/<name>/logs/slurm/pgs_<name>_<step>_<arrayid>_<taskid>.out`

### Constraints
- **Prep must run alone**: `--sbatch --steps prep` cannot be combined with other steps.
- **Sumstat requires -i**: Per-sumstat steps require `-i <sumstat_dir>`.
- **No interactive arrays**: SLURM arrays are submitted and the driver waits for
  completion using `squeue` polling.

### Failure handling (chromosome-level)

The system uses a **continue-on-failure** strategy: individual chromosome failures
do not abort the entire run. Instead, placeholder files are written and downstream
steps can proceed (producing partial results).

#### Current behavior by failure scenario

| Scenario | Driver behavior | Step behavior | Final output |
|----------|-----------------|---------------|--------------|
| **One chromosome fails** | Logs warning, continues | Writes placeholder `.snpRes` (header-only) + `FAILED` marker | Partial score (missing that chr) |
| **Several chromosomes fail** | Logs warning per failure, continues | Placeholders for each failed chr | Partial score (missing failed chrs) |
| **All chromosomes fail** | Logs warning, continues to next step | All placeholders | Empty `scores.tsv.gz` (header-only) |

#### Placeholder mechanism

When a chromosome fails at any step, the step writes:
1. **Placeholder output file**: Header-only file (e.g., `chr5.snpRes` with just the column headers).
2. **FAILED marker**: File indicating failure reason (e.g., `work_chr5/FAILED` containing `sbayesr_failed`).

This allows downstream steps to:
- Detect the failure (by checking for `FAILED` markers or header-only files).
- Skip processing for that chromosome without crashing.
- Produce partial results from successful chromosomes.

#### Failure markers by step

| Step | Placeholder file | Failure marker |
|------|------------------|----------------|
| filter-variants | `filtered/chr<N>_filtered.tsv` | `filtered/FAILED_chr<N>` |
| calc-posteriors | `posteriors/chr<N>.snpRes` | `posteriors/work_chr<N>/FAILED` |
| format-posteriors | `posteriors_mapped/chr<N>.snpRes` | `posteriors_mapped/FAILED_chr<N>` |
| calc-score | `scores/chr<N>.sscore` | `scores/FAILED_chr<N>` |
| combine-scores | (no placeholder; writes header-only `scores.tsv.gz` if all empty) | — |

#### Driver-level monitoring

The driver job:
1. **Submits array** and polls `squeue` until all tasks finish.
2. **Counts failures**: Uses `sacct` to detect FAILED/TIMEOUT/OOM tasks.
3. **Logs warning**: `"Warning: SLURM array job X for step 'posteriors' had failed/unknown task(s)."`.
4. **Continues**: Does not abort; downstream arrays are submitted.
5. **Sanity-checks outputs**: After each array, verifies expected files exist; logs warning if missing.

#### Failure reasons tracked in FAILED markers

- `filtered_sumstat_missing`: No filtered sumstat input for this chromosome.
- `format_for_sbayesr_failed`: Could not generate sbayesR input file.
- `ma_missing`: sbayesR input file was not created.
- `ldref_missing`: LD reference files not found for chromosome.
- `sbayesr_failed`: gctb/sbayesR exited non-zero.
- `sbayesr_no_output`: gctb ran but produced no `.snpRes` output.
- `posteriors_missing`: Posteriors file missing for format step.
- `mapped_posteriors_missing`: Mapped posteriors missing for score step.

#### Implications

- **Partial results are valid**: A run with 20/22 successful chromosomes produces a
  usable (if incomplete) PGS score.
- **All-fail produces empty output**: If every chromosome fails, `scores.tsv.gz` contains
  only headers; downstream analysis should detect this.
- **Logs are essential**: Check `FAILED` markers and SLURM logs to diagnose failures.
- **No automatic retry**: Failed chromosomes are not retried; manual re-run is required.

### Intended workflow

```bash
# Step 1: Prep (shared across all sumstats)
./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch

# Step 2: Per-sumstat (can submit many in parallel; each gets its own driver job)
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_814 --sbatch
./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_815 --sbatch
# ... etc
```

---

## Criteria checklist
- Mapfile uses a single `chr` and `pos` for all sources.
- Mapfile contains all three source SNP IDs and their allele columns.
- Prep step is sumstat-agnostic and builds a union map for geno + ldref only.
- SNP inclusion lists are always derived from the mapfile (prep or sumstat mapfile).
- Sumstat step attaches `sumstat_*` via `chr/pos + alleles`.
- Sumstat `EAF` is filled only from LD reference EAF (`ldref_a2freq`) when missing;
  `EAF_1KG` is not used.
- INFO/MAF reference files are replaced by a user-provided inclusion list, with an
  explicit ID-space specifier (`ss`, `ld`, or `gt`), applied after mapfile reduction.
- Conversions to/from LD reference and genotype IDs only use the mapfile.
- Final output includes the full mapfile (`variant_map.tsv.gz`).

