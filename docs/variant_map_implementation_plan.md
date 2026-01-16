% Variant Map Implementation Plan (Remaining Work)

## Scope
Track only the **remaining** implementation updates needed to align code with
`docs/variant_map_mapping_plan.md`.

## Already implemented
- User inclusion lists (`filters.inclusion_list.{gt,ss,ld}`) parsed and mounted.
- Sumstat reduction produces `sumstat_for_posteriors.tsv.gz` and fills missing EAF from LDref.
- LDref SNP IDs flow into sbayesR input via `LDREF_SNPID`, and posteriors map LDref -> genotype.
- Output `variant_map.tsv.gz` written during finalize.
- **Format-sumstat single-pass chromosome split** (2026-01-15):
  - `format_sumstat.sh` now outputs per-chromosome files (`formatted/chrN.tsv`).
  - NA filtering and chromosome splitting in a single awk pass.
- **Per-chromosome filter-variants** (2026-01-15):
  - `filter_variants.sh` supports `specific_chr` parameter.
  - `run_filter_variants_chr()` processes a single chromosome.
  - `build_sumstat_map_and_reduce_chr()` filters mapfile to target chromosome.
  - Concatenation helpers for per-chr mapfiles and filtered sumstats.
- **Sumstat SLURM array in driver** (2026-01-15):
  - Driver runs format-sumstat directly, then submits filter-variants as array.
  - `slurm.sumstat` config section with defaults: `{ mem: 1g, cpus: 1, time: '00:30:00', max_parallel: 22 }`.
  - Output sanity checking for sumstat array.

## Remaining updates

### 1) Prep mapfile should be **intersection**, not union
Files to update:
- `bin/lib/steps/prep_inclusion_list.sh`

Changes:
- Switch base `prep/variant_map.tsv` to **intersection** of genotype + LD reference
  (chr/pos + alleles).
- Stop emitting unmatched geno-only or ldref-only rows in `prep/variant_map.tsv`.

### 2) Emit per-chromosome mapfiles in prep
Files to update:
- `bin/lib/steps/prep_inclusion_list.sh`

Changes:
- Write `prep/variant_map/chrN.tsv` (one file per chromosome), same schema as `prep/variant_map.tsv`.
- Ensure inclusion list derivation uses the base mapfile (intersection).

### 3) Update tests / validation checklist
- Validate prep mapfile size ~= intersection (geno ∩ ldref), not union.
- Validate per-chromosome mapfiles exist and sumstat mapfile is concatenated from them.
- Confirm EAF fill (LDref) still works under per-chromosome processing.


