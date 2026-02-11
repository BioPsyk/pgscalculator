% General Design

## Purpose
Define the required behavior for variant ID mapping so that sumstat, LD reference, and genotype
identifiers are linked unambiguously and conversions never rely on RSIDs alone.

## Canonical mapfile schema
Single shared `chr` with **dual position columns** for build support, plus per-source SNP IDs and allele columns:
- `chr`, `pos_b37`, `pos_b38`
- `sumstat_snpid`, `sumstat_effectallele`, `sumstat_otherallele`
- `geno_snpid`, `geno_a1`, `geno_a2`
- `ldref_snpid`, `ldref_a1`, `ldref_a2`

Notes:
- No compound `chrpos` field is required.
- No sorting requirement is imposed for runtime use.
- The mapfile may include non-key frequency columns used for EAF filling/auditing.
- Both `pos_b37` and `pos_b38` are always populated from liftover reference files during prep.
- The mapfile contains **all** LD reference variants that could be lifted over; genotype columns are `NA` where no genotype match, and sumstat columns are `NA` where no sumstat match (per-sumstat step).

## Genome build handling

### Config options

```yaml
# Genome build of the genotype files (default: GRCh37)
genotype_build: GRCh37  # or GRCh38

# Liftover reference file (required if genotype_build: GRCh38)
# This file provides chr:pos mapping between GRCh38 and GRCh37, derived from dbSNP/cleansumstats
liftover_reference: /path/to/references/liftover/dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt.gz
```

### Liftover reference file

The liftover reference is derived from the cleansumstats pipeline (dbSNP-based) and stored in
`references/liftover/`. This file provides chr:pos mapping between builds for all known variants.

**File locations:**
```
references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz  # sorted on col1 (b37)
references/liftover/dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt.gz  # sorted on col1 (b38)
```

**Format (space-delimited, sorted on column 1 with `LC_ALL=C`):**

For `dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz` (b37 → b38):
```
10:1045940 10:1000000 rs1831596373 A C
...
```
- Column 1: `chr:pos_b37` (sort key)
- Column 2: `chr:pos_b38`

For `dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt.gz` (b38 → b37):
```
10:1000000 10:1045940 rs1831596373 A C
...
```
- Column 1: `chr:pos_b38` (sort key)
- Column 2: `chr:pos_b37`

**Common columns (both files):**
1. `chr:pos` of source build - **sort key**
2. `chr:pos` of target build
3. `rsid` - dbSNP RSID (matches sumstat RSIDs)
4. `a1` - Allele 1
5. `a2` - Allele 2 (may contain multiple alleles, e.g., "C,G")

**Note:** Each file is pre-sorted on column 1 using `LC_ALL=C sort`, enabling
efficient `join` operations with genotype positions in the corresponding build.

### Build assumptions
- **LD reference**: Always GRCh37 (sbayesR ukb_10k_hm3 is GRCh37).
- **Sumstats**: Both builds available via cleansumstats (`cleaned_GRCh37.gz`, `cleaned_GRCh38.gz`).
- **Genotypes**: Configurable via `genotype_build` (default: GRCh37).
- **Liftover reference**: Used to pre-compute augmented LD reference with both positions.

### Pre-computed augmented LD reference (optimization)

To avoid redundant liftover operations, we pre-compute an augmented LD reference that
contains both `pos_b37` and `pos_b38` for all LD reference variants. This augmentation
is performed as part of `--steps prep` (specifically within `prep-ldref`) and is cached
for reuse across subsequent runs.

**Augmentation within `prep-ldref`** (runs once, cached)

**IMPORTANT: Efficiency considerations**
- The liftover reference files are already sorted on column 1 using `LC_ALL=C sort`.
- Do NOT re-sort the liftover file - use it directly with `join`.
- Do NOT process chromosomes separately - this would require 22 scans of the huge liftover file.
- Instead: concatenate all LD ref variants across chromosomes, join once, then split results.

```bash
# Step 1: Concatenate all LD reference files and sort on chr:pos_b37
for chr in $(seq 1 22); do
  cat prep/ldref/chr${chr}_ld_rsids
done | LC_ALL=C sort -k1,1 > prep/ldref/all_ld_rsids_sorted.tsv

# Step 2: Join with pre-sorted liftover file (single pass!)
# Input: all LD ref positions (sorted on chr:pos_b37)
# Liftover: chr:pos_b37, chr:pos_b38, rsid, a1, a2 (pre-sorted on col1, space-separated)
# Output: Augmented LD reference with both positions
# Note: join uses whitespace as default separator - no need for tr

LC_ALL=C join \
  prep/ldref/all_ld_rsids_sorted.tsv \
  <(zcat dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz) \
  > prep/ldref_augmented/all_ld_augmented.tsv

# Step 3: Split by chromosome
awk -F'\t' '{
  chr = $1; sub(/:.*/, "", chr)
  print >> "prep/ldref_augmented/chr" chr "_ld_augmented.tsv"
}' prep/ldref_augmented/all_ld_augmented.tsv

# Output columns: chr:pos_b37, ldref_a1, ldref_a2, ldref_rsid, chr:pos_b38, liftover_rsid, liftover_a1, liftover_a2
```

**Key insight**: By concatenating all LD ref variants first, we only scan the liftover file ONCE
instead of 22 times. This is critical for performance since the liftover file is ~900M lines.

**Output location:**
```
prep/ldref_augmented/chr{N}_ld_augmented.tsv
```

**Augmented LD reference schema:**
- `chr:pos_b37` - Original LD reference position (GRCh37)
- `chr:pos_b38` - Lifted position (GRCh38) from liftover reference
- `rsid` - dbSNP RSID
- `a1`, `a2` - Alleles

**Benefits:**
- One-time computation per LD reference (not per genotype dataset)
- Prep step just joins genotypes against pre-augmented LD reference
- No liftover join needed during per-genotype prep
- Sumstat matching uses `pos_b38` directly from pre-computed map

### Matching strategy by genotype build

| Genotype build | Prep matching key | pos_b37 source | pos_b38 source |
|----------------|-------------------|----------------|----------------|
| **GRCh37** (default) | chr:pos_b37 + alleles | augmented LD ref | augmented LD ref |
| **GRCh38** | chr:pos_b38 + alleles | augmented LD ref | augmented LD ref |

**Note:** Both position columns always come from the pre-augmented LD reference.
The liftover join is done once during `prep-ldref` (part of `--steps prep`), cached for reuse.

**Key insight**: By pre-augmenting the LD reference with both positions, the per-genotype prep
step only needs a simple join - no liftover processing required during the main workflow.

### Implementation details

**When running `--steps prep`:**

1. **prep-ldref** (includes augmentation):
   - Extract positions and RSIDs from LD reference (always GRCh37).
   - Join with liftover file (`dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz`) to add `pos_b38`.
   - Output: `prep/ldref_augmented/chr{N}_ld_augmented.tsv` with both positions.
   - **Cached**: If augmented files already exist, skip the liftover join.

2. **prep-genotypes**: Extract positions from genotypes as-is (in their native build).
   - Output: `chr:pos`, `a1`, `a2`, `variant_id` per chromosome.

3. **prep-inclusion-list** (now simplified):
   - If `genotype_build == GRCh37`:
     - Match genotypes ↔ augmented LD ref by **chr:pos_b37 + alleles**.
   - If `genotype_build == GRCh38`:
     - Match genotypes ↔ augmented LD ref by **chr:pos_b38 + alleles**.
   - Both `pos_b37` and `pos_b38` come from the pre-augmented LD ref.
   - **All** augmented LD ref variants are kept (left join from LD ref); genotype columns are `NA` where no genotype match.
   - No liftover join needed during this step.

4. **format-sumstat**: Simplified - just split `cleaned_GRCh38.gz` by chromosome.
   - No coordinate mapping needed since variant_map has both positions.

5. **filter-variants**: Match sumstat to variant_map using `pos_b38` (sumstat native coordinates).
   - Both `pos_b37` and `pos_b38` are already populated from prep step.
   - No coordinate conversion needed - direct matching on GRCh38.

6. **calc-score**: Uses `geno_snpid` which is build-independent.

### Augmented LD reference join strategy (within prep-ldref)

The augmentation (performed once during `--steps prep`) uses efficient unix `join` with
`LC_ALL=C` for fast matching. The liftover reference `dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz`
is pre-sorted on column 1 (`chr:pos_b37`).

**Conceptual workflow:**
```bash
# For each chromosome, join LD ref (b37) with liftover to add b38 positions
for chr in {1..22}; do
  # LD ref format: chr:pos_b37, a1, a2, rsid (already sorted or sort here)
  # Liftover format: chr:pos_b37, chr:pos_b38, rsid, a1, a2 (pre-sorted on col1)
  
  LC_ALL=C join -1 1 -2 1 \
    <(LC_ALL=C sort -k1,1 ldref/chr${chr}_ld_rsids.tsv) \
    <(zcat dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz | grep "^${chr}:" | tr ' ' '\t') \
    > ldref_augmented/chr${chr}_ld_augmented.tsv
done

# Output: chr:pos_b37, ldref_a1, ldref_a2, rsid, chr:pos_b38, liftover_rsid, liftover_a1, liftover_a2
```

### Genotype ↔ Augmented LD ref matching (prep-inclusion-list)

The mapfile must retain **all** augmented LD ref rows (left join from LD ref). Use a join that keeps every augmented LD ref variant and sets genotype columns to `NA` where no genotype match.

**For GRCh37 genotypes:**
```bash
# Match on chr:pos_b37 (column 1 in both files)
LC_ALL=C join -1 1 -2 1 \
  <(LC_ALL=C sort -k1,1 genotypes_b37.tsv) \
  <(LC_ALL=C sort -k1,1 ldref_augmented/chr${chr}_ld_augmented.tsv) \
  > variant_map_chr${chr}.tsv
```

**For GRCh38 genotypes:**
```bash
# Match on chr:pos_b38 (column 1 in genotypes, column 5 in augmented LD ref)
# Need to re-key augmented LD ref on pos_b38 for join
LC_ALL=C join -1 1 -2 1 \
  <(LC_ALL=C sort -k1,1 genotypes_b38.tsv) \
  <(awk -F'\t' '{print $5, $0}' ldref_augmented/chr${chr}_ld_augmented.tsv | LC_ALL=C sort -k1,1) \
  > variant_map_chr${chr}.tsv
```

**Allele matching considerations:**
- After position join, verify alleles match (with strand flip support).
- The liftover reference `a2` column may contain multiple alleles (e.g., "C,G");
  match if genotype allele is any of the listed alleles.

### Fallback behavior

If `genotype_build: GRCh38` is set but liftover_reference is not provided or file is missing:
1. Log an error: "liftover_reference required when genotype_build is GRCh38".
2. Exit with non-zero status.

If liftover produces few matches (e.g., <10% of genotype variants):
1. Log a warning about low liftover rate.
2. Continue processing with matched variants.

## Filtering order (planned)
1) **prep** (`--steps prep`) builds base mapfiles from **all** LD reference variants that could be lifted over (augmented LD ref),
   including LD reference EAF (`ldref_a2freq`) and both position columns (`pos_b37`, `pos_b38`).
   - `prep-ldref` creates augmented LD reference with both positions (via liftover, cached).
   - `prep-genotypes` extracts genotype positions in their native build.
   - `prep-inclusion-list` joins genotypes onto augmented LD ref (left join from LD ref); genotype columns are `NA` where no match.
   - **All** augmented LD ref variants are kept in the mapfile.
   - Write **chromosome-specific mapfiles** only (e.g., `prep/variant_map/chrN.tsv`) for parallel use.
   - Do **not** create a combined prep mapfile during processing.
2) **sumstat** attaches `sumstat_*` columns via `chr/pos_b38 + alleles` into **chromosome-specific**
   sumstat mapfiles (using the prep `prep/variant_map/chrN.tsv` files).
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
### Prep step (`--steps prep`, sumstat-agnostic)
1) `prep-ldref`: Create augmented LD reference with both `pos_b37` and `pos_b38`.
   - Joins LD ref (b37) with liftover file to add `pos_b38`.
   - Cached: if augmented files exist, skip liftover join.
2) `prep-genotypes`: Extract positions from genotypes in their native build.
3) `prep-inclusion-list`: Build per-chromosome mapfiles from **all** augmented LD ref variants (left join from LD ref).
   - Match by `pos_b37` (if genotype_build=GRCh37) or `pos_b38` (if genotype_build=GRCh38); genotype columns `NA` where no match.
   - Both `pos_b37` and `pos_b38` come from augmented LD ref.
4) Populate geno/ldref columns, dual positions, and LD reference EAF (`ldref_a2freq`).
3) No sumstat columns are added at prep, and no combined prep mapfile is required.
4) The prep step should support **per-chromosome parallelism** when creating the
   base mapfiles, using **SLURM arrays** (config: `slurm.prep.max_parallel`).
5) Always derive SNP inclusion lists from the mapfile:
   - Use the prep mapfile when the prep step needs an inclusion list.
   - Use the sumstat-annotated mapfile when the sumstat step needs an inclusion list.
   - It is acceptable to derive a union inclusion list that contains all three sources, which
     can be useful for building the posterior-calculation input.
6) Timing/logging: each step logs start + completion with elapsed time, and the overall
   pipeline logs a total elapsed time.

### Format-sumstat step (simplified with dual-position variant_map)

Since the variant_map now contains both `pos_b37` and `pos_b38` (populated via liftover references
during prep), the sumstat can be matched directly using its native GRCh38 coordinates.

1. Read `cleaned_GRCh38.gz` directly (no need to paste with `cleaned_GRCh37.gz`).
2. In a single awk pass:
   - Split output by chromosome into per-chromosome files (`formatted/chrN.tsv`).
3. Output: Per-chromosome sumstat files ready for parallel processing.

**Note:** The `cleaned_GRCh37.gz` file is no longer needed for coordinate mapping since the
variant_map provides both positions. This simplifies processing and saves time.

### Sumstat step (sumstat-specific)
All sumstat processing steps should support **per-chromosome parallelism**.
Each chromosome process should only load its corresponding `prep/variant_map/chrN.tsv`.

1) Create sumstat-specific per-chromosome mapfiles.
2) Attach `sumstat_*` columns using `chr/pos_b38 + alleles` to match against the variant_map.
   - Sumstat uses native GRCh38 coordinates; variant_map has `pos_b38` from prep.
3) Fill missing sumstat `EAF` from mapfile `ldref_a2freq` when needed (allele-aware).
4) The mapfile row set is unchanged (all LD ref liftover variants); sumstat_* columns are attached, `NA` where no sumstat match.
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

## Scoring (no separate inclusion list)
There is **no separate inclusion list for scoring**. Only variants that were sent to the posterior step (filtered sumstat matched to the mapfile) are in the posteriors file; those variants are already restricted to genotype-matched mapfile rows. The score step uses the **posteriors file** to obtain variant IDs for plink `--extract` and `--score`. No prep-derived inclusion list is required for scoring.

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

## Output files (v2)

### sumstat_augmented.tsv.gz vs augmented_sumstat.gz

| | **sumstat_augmented.tsv.gz** | **augmented_sumstat.gz** |
|--|------------------------------|---------------------------|
| **Row set** | Variant map intersection (one row per variant present in both the mapfile and the formatted sumstat — **not** all sumstat variants; inner join on sumstat ⋈ variant_map). | Variant map (one row per variant in the LD reference that could be lifted over). |
| **Columns** | All columns from the formatted sumstat plus GENO_ID, POST_EFFECT, POST_PIP, IN_ANALYSIS. | Reduced: RSID, EffectAllele, OtherAllele, B, SE, Z, P, MAF, postEffect, benchEffect only. |
| **Purpose** | Full augmented sumstat for internal/debug use; row set restricted to mapfile variants for performance and relevance. | User-facing, compact file for auditing; same row set as variant_map so users can join on RSID. |
| **Built from** | Per-chr matched files from filter-variants (already mapfile-restricted, with LDREF_SNPID + GENO_ID) ⋈ posteriors (left join). ⋈ variant_map (inner join) ⋈ posteriors. | variant_map ⋈ sumstat (B,SE,Z,P) ⋈ posteriors ⋈ MAF; uses sumstat_augmented as source for B,SE,Z,P. |

### augmented_sumstat.gz

Per-sumstat file for auditing and back-tracing. **Same row set as the variant map** (one row per variant in `variant_map.tsv.gz`), so users can join on `RSID`. Contains only the sumstat columns needed for interpretation plus calculated MAF, posterior effect, and (when available) benchmark effect. CHR, POS, and genoID are omitted because **RSID is the key** for lookups in the variant map.

**Schema:**
```
RSID  EffectAllele  OtherAllele  B  SE  Z  P  MAF  postEffect  [benchEffect]
```

- **From sumstat:** `RSID`, `EffectAllele`, `OtherAllele`, `B`, `SE`, `Z`, `P` (NA where sumstat did not match).
- **Added:** `MAF` (calculated from genotypes; NA where no genotype match).
- **Added:** `postEffect` (posterior effect).
- **Added when feature exists:** `benchEffect` (benchmark effect).

Users who need chr/pos or genotype IDs join on `RSID` against the variant map.

**Memory use (finalize-output):** Building `augmented_sumstat.gz` used to load the full `sumstat_augmented.tsv.gz` (B, SE, Z, P by RSID) into awk arrays keyed by RSID, then stream the variant map and look up. On large sumstats (e.g. hundreds of MB uncompressed) that caused finalize to OOM even with 16g. **Fix:** use a Unix **sort + join** pipeline: extract (key, B, SE, Z, P) from the sumstat in a streaming way and sort by key; sort the variant map by sumstat_snpid; `join` on the key so the sumstat is never loaded into memory. Only posteriors and MAF (small) are kept in memory for the final column assembly. This keeps finalize memory low and independent of sumstat size.

### main_raw_score_all.gz

Primary PGS score output. Name and column set are unchanged from v1.

**Schema:**
```
IID  ALLELE_CT  NAMED_ALLELE_DOSAGE_SUM  SCORE1_AVG  SCORE1_SUM  FILE_SUM
```

- `IID`: Sample identifier (FID is not included).
- `ALLELE_CT`, `NAMED_ALLELE_DOSAGE_SUM`, `SCORE1_AVG`, `SCORE1_SUM`, `FILE_SUM`: as in current v1.

### variant_map.tsv.gz

Sumstat-specific variant map for joining output files. **Column 1 is RSID** (from the liftover reference), so users can use it as a single key to map between augmented sumstat, score inputs, and any other outputs.

**Row set:** **All** variants in the LD reference that could be lifted over (full augmented LD ref). Variants not found in genotype have genotype columns set to `NA`; variants not matched when joining the sumstat have sumstat columns set to `NA`.

**Schema (RSID as col1, then canonical mapfile columns):**
```
rsid  chr  pos_b37  pos_b38  sumstat_snpid  sumstat_effectallele  sumstat_otherallele  geno_snpid  geno_a1  geno_a2  ldref_snpid  ldref_a1  ldref_a2  ldref_a2freq
```

- `rsid`: From liftover reference; **key for user-facing joins** between output files.
- Remaining columns: same as current v2 sumstat-specific variant map (chr, dual positions, per-source SNP IDs and alleles, LD ref EAF).

## Benchmark weights and scores (v2)

The benchmark provides a comparison PGS using **observed GWAS effects** (no Bayesian shrinkage) on a **pruned variant set**, following the same ideas as v1.

### Purpose
- Compare the main (sBayesR) score against a simple sum-of-effects score on LD-pruned variants.
- Use the same genotype data and the same filtered sumstat (effect estimates), but no posterior step.

### Benchmark weights
- **Source:** Filtered sumstat (same as posterior input): variant ID, allele, and effect (B or BETA).
- **Mapping:** RSID (or sumstat SNP ID) → genotype ID via the variant map (or prep-derived RSID→geno list for genotypes-matched variants only).
- **No shrinkage:** Weights are the raw GWAS effects; no sBayesR posteriors.

### Benchmark score computation (v1-style, carried into v2)
1. **Input:** Filtered sumstat (per chromosome) and variant map (or RSID→geno list) so that each sumstat row has a genotype ID.
2. **Restrict to genotype-matched variants:** Only variants present in genotypes are scored.
3. **MAF filter:** Apply a MAF threshold (e.g. 0.05) so that very rare variants are excluded from the benchmark.
4. **LD pruning:** Run plink `--indep-pairwise` (e.g. 250 50 0.25) on the genotype-matched variant set to obtain a pruned list.
5. **Benchmark sumstat:** Build a 3-column file (genotype ID, allele, effect) for **pruned variants only**.
6. **Score:** Run plink `--score` with that file (no dosage denominator in v1-style; or cols=scoresums) to obtain per-sample benchmark scores.
7. **Combine:** Aggregate per-chromosome benchmark scores into a single benchmark score file (e.g. `benchmark.sscore` or combined table).

### Outputs
- **Benchmark weights:** Effectively the filtered sumstat restricted to pruned, genotype-matched variants, with genotype ID for scoring.
- **Benchmark scores:** Per-sample scores (e.g. `benchmark.sscore` or equivalent), comparable to the main score for correlation/QC.

### Integration with augmented sumstat
- When the benchmark step has been run, **benchEffect** can be filled in the augmented sumstat (effect used in the benchmark score for that variant, or NA if not in the pruned set). This allows users to compare posterior vs benchmark effect per variant.

### Benchmark configuration
- **maf_threshold:** MAF filter for benchmark variant set (e.g. 0.05).
- **indep_pairwise:** LD pruning parameters for plink `--indep-pairwise` as `[window_kb, step, r2]` (e.g. 250, 50, 0.25). Configurable in `config.yaml` under `benchmark:`.

## Weights step (renamed from posteriors)

The step that produces scoring weights is named **weights** (not "posteriors") so that it covers both (1) **sBayesR posterior weights** and (2) **benchmark weights** (observed effects, LD-pruned). Both are run in the same phase and can be executed in parallel.

### Two array jobs, independent resource config
- **sBayesR weights array:** Runs `calc-posteriors` and `format-posteriors` per chromosome (posterior effects from sBayesR, mapped to genotype IDs). Resource configuration: **`slurm.weights_sbayesr`** (mem, cpus, time, max_parallel).
- **Benchmark weights array:** Runs `calc-benchmark` per chromosome (observed effects, MAF filter, LD pruning, benchmark score). Resource configuration: **`slurm.weights_benchmark`** (mem, cpus, time, max_parallel).

The driver submits **two separate SLURM array jobs** when `--steps weights` is requested, so sBayesR and benchmark can be tuned independently (e.g. more memory/cpus for sBayesR, lighter resources for benchmark). Both arrays run in the same step phase; the score step uses the sBayesR weights (posteriors_mapped); benchmark scores are written alongside.

### Step group and CLI
- **Step group:** `weights` (replaces `posteriors`). Running `--steps weights` runs: calc-posteriors, format-posteriors, calc-benchmark (all per-chromosome where applicable).
- **Concrete steps** remain: `calc-posteriors`, `format-posteriors`, `calc-benchmark` for direct invocation.

## Final output
Only at final output time, combine per-chromosome files into consolidated outputs
(`variant_map.tsv.gz`, augmented sumstat, main score file, etc.) for auditing and back-tracing.

## Separate finalize step (final output creation) — implemented

**Problem:** Today `score` includes calc-score (array) plus combine-scores and finalize-output run **in the driver process** after the score array. Finalize-output is memory-heavy (large awk over variant_map + full augmented sumstat); with a small driver (e.g. 5g) it can OOM on large sumstats.

**Idea:** Introduce a dedicated **finalize** step group so final output creation is a separate phase and can be given its own resources (including higher memory when submitted as a SLURM job).

### 1. Step group and order

- **score** = `calc-score` only (per-chromosome scoring; array when using --sbatch).
- **finalize** = `combine-scores` + `finalize-output` (single, non-array: merge chr scores, then build all final deliverables).
- **STEP_GROUP_ORDER:** `prep` → `sumstat` → `weights` → `score` → **`finalize`**.

So:

- `--steps score` runs only the calc-score step (and, in current implementation, nothing merges scores yet; see below).
- `--steps score,finalize` (or `--steps sumstat,weights,score,finalize`) runs calc-score then combine-scores then finalize-output.
- `--all` would run prep, sumstat, weights, score, **finalize**.

Concrete steps inside **finalize**:

- **combine-scores:** Merge per-chr `work/scores/chr*.sscore` → `scores.tsv.gz`, `main_raw_score_all.gz`.
- **finalize-output:** Combine posteriors, build `sumstat_augmented.tsv.gz`, `augmented_sumstat.gz`, `variant_map.tsv.gz`, copy config, run summary, stepwise details.

#### Joins in the finalize step

All joins in finalize are done so that **no large dataset is fully loaded into memory**; large inputs are either streamed or joined via Unix `sort` + `join`.

| Sub-step | Inputs | Join strategy | Memory |
|----------|--------|----------------|--------|
| **combine-scores** | Per-chr `chr*.sscore` | No join; concatenate/aggregate score files. | One stream at a time. |
| **Combine posteriors** | Per-chr `chr*.snpRes`, `ldref_to_genoid.tsv` | No join; concatenate chr files. Optional GENO_ID via small in-memory lookup (RSID → geno_snpid). | Small (lookup only). |
| **sumstat_augmented.tsv.gz** | Per-chr matched files (`filtered/chr*_matched.tsv`), posteriors_combined | **Approach B (reuse filter-variants output):** Concatenate per-chr matched files (already mapfile-restricted, ~1M rows total, with LDREF_SNPID + GENO_ID). Rearrange so LDREF_SNPID is col 1, sort, left-join with posteriors. No re-sorting of the full formatted sumstat. | Concatenation + sort of ~1M rows; no large file processed. |
| **augmented_sumstat.gz** | sumstat_augmented.tsv.gz, variant_map, posteriors, MAF | **Sort + Unix join:** (1) Extract (key, B, SE, Z, P) from sumstat_augmented in a streaming awk; sort by key. (2) Sort variant_map by col4 (sumstat_snpid). (3) Unix `join` (tab-separated, -1 4 -2 1) on key → one row per variant_map row with B,SE,Z,P. (4) One awk adds MAF and postEffect via small in-memory lookups (posteriors by ldref_snpid, MAF by geno_snpid). Row set = variant_map; the large sumstat is never loaded. | Small (posteriors + MAF only); sumstat only in sort/join temp files on disk. |
| **variant_map.tsv.gz** | variant_map.tsv (from prep or sumstat) | No join; reorder columns so RSID (ldref_snpid) is col1, then gzip. | One line at a time. |

Since `sumstat_augmented.tsv.gz` is built from the per-chr matched files (already mapfile-restricted during filter-variants), its row count matches the mapfile (~1M variants) rather than the full sumstat (~17M). The expensive chr:pos+alleles matching was already done per-chromosome during filter-variants; finalize just concatenates and joins with posteriors. This keeps both `sumstat_augmented.tsv.gz` itself and the downstream `augmented_sumstat.gz` (which reads it) fast.

#### Plan: All joins via Unix join

**Goal:** Every join in finalize uses Unix `sort` + `join` only; no in-memory hash lookups for join keys. Memory stays flat regardless of file size.

**Conventions:** Tab-separated. **All sort and join in finalize must use `LC_ALL=C`** for locale-independent, byte-order-consistent ordering and matching (same as prep): run `LC_ALL=C sort ...` and `LC_ALL=C join ...` for every sort and join. Example: `LC_ALL=C sort -t $'\\t' -k<keycol>,<keycol>`; `LC_ALL=C join -t $'\\t' ...`. Strip headers before sort/join; prepend header to final output. One temp dir for intermediates.

**1. Combine posteriors (add GENO_ID)** — One join. Concatenate chr posterior files (skip headers), output RSID, A1, A2, FREQ, EFFECT, SE, PIP; `LC_ALL=C sort` by RSID. `LC_ALL=C sort` ldref_to_genoid by RSID. `LC_ALL=C join -1 1 -2 1` on RSID → RSID, GENO_ID, A1, A2, FREQ, EFFECT, SE, PIP. Prepend header. No in-memory lookup.

**2. sumstat_augmented.tsv.gz** — **Approach B: reuse filter-variants output.** Filter-variants saves per-chr matched files (`filtered/chr*_matched.tsv`) containing all sumstat columns + LDREF_SNPID + GENO_ID, already restricted to mapfile variants (~1M rows). At finalize: concatenate matched files, rearrange so LDREF_SNPID is col 1, `LC_ALL=C sort`, left-join with posteriors (`-a 1 -e NA`), emit sumstat cols + GENO_ID + POST_EFFECT + POST_PIP + IN_ANALYSIS; gzip. No re-sorting of the full formatted sumstat (~17M). **Key design rule: `sumstat_augmented.tsv.gz` must contain only mapfile variants, not the full sumstat.** This avoids sorting 17M+ rows and keeps downstream `augmented_sumstat.gz` (which reads `sumstat_augmented.tsv.gz`) fast.

**3. augmented_sumstat.gz** — Three joins. (a) variant_map `LC_ALL=C sort` by sumstat_snpid; sumstat_augmented: extract key, B, SE, Z, P (streaming), `LC_ALL=C sort` by key. `LC_ALL=C join -1 4 -2 1` → variant_map + B, SE, Z, P; output with ldref_snpid in known col, `LC_ALL=C sort` by ldref_snpid. (b) posteriors: RSID, EFFECT; `LC_ALL=C sort` by RSID. `LC_ALL=C join` on ldref_snpid/RSID → + postEffect; `LC_ALL=C sort` by geno_snpid. (c) MAF `LC_ALL=C sort` by geno_snpid. `LC_ALL=C join` on geno_snpid → + MAF. One awk to emit final 10 columns; prepend header; gzip. No in-memory lookups.

**4. variant_map.tsv.gz** — No join; reorder cols so RSID is col 1, gzip.

**Efficiency:** (1) Sort each input file once; reuse the same sorted file if it is joined multiple times (e.g. posteriors sorted by RSID can feed both sumstat_augmented and augmented_sumstat pipelines). (2) When writing join output, put the *next* join key in column 1 so the next step is a single `LC_ALL=C sort -k1,1` with no column reordering. (3) Use one temp dir for all intermediates; delete at end of finalize-output. (4) For sumstat/formatted inputs, use a single streaming awk to extract or normalize (key, rest); avoid reading the full file into memory. (5) combine-scores stays concatenate-only (no join). (6) Use `join -a 1 -e NA` only where the driver row set must be preserved (e.g. variant_map as driver in augmented_sumstat.gz); inner join elsewhere is simpler and sufficient where missing rows are acceptable. **(7) Critical: `sumstat_augmented.tsv.gz` must contain only mapfile variants, not the full sumstat. The primary strategy (Approach B) achieves this by reusing the per-chr matched files saved during filter-variants — no re-join with the formatted sumstat needed.  Including all sumstat variants (~17M) would make sorting and all downstream joins needlessly slow; the mapfile intersection (~1M) is the relevant row set.** **(8) Filter-variants saves `filtered/chr*_matched.tsv` (pre-QC matched sumstat with LDREF_SNPID + GENO_ID) alongside the filtered outputs. This small I/O cost during filter-variants eliminates the most expensive finalize operation (sort+join of the full formatted sumstat).**

### 2. SLURM behaviour (--sbatch)

- **Driver** runs: format-sumstat → sumstat array → weights_sbayesr array → weights_benchmark array → score array. It does **not** run combine-scores or finalize-output in-process.
- After the **score array** completes, the driver **submits one SLURM job** for the **finalize** step (single job, no array), with its own resources, e.g.:
  - **slurm.finalize:** `{ mem: 16g, cpus: 1, time: 1:00:00 }` (configurable; higher mem for large sumstats).
- The driver **waits** for this finalize job to complete (same pattern as for arrays), then exits. Optionally it could run cleanup if `--cleanup` is set.

So:

- Driver stays small and only orchestrates; it never runs the heavy finalize-output.
- Final output creation always runs in a dedicated job with `slurm.finalize`, avoiding driver OOM.

### 3. Config

```yaml
slurm:
  # ... existing keys ...
  score:    { mem: 10g, cpus: 4, time: '01:00:00', max_parallel: 22 }
  finalize: { mem: 16g, cpus: 1, time: '1:00:00' }   # single job, no array
```

### 4. CLI (implemented: explicit finalize)

- **Default for per-sumstat pipeline:** Users today run `--steps sumstat,weights,score`. With the new split, to get the same end-to-end result they would run `--steps sumstat,weights,score,finalize` (or we keep “score” in the CLI to mean “score + finalize” for backward compatibility and only split internally for SLURM; see below).
- **Option A (explicit finalize):** Require `finalize` in `--steps` when final outputs are desired. Docs and examples use `sumstat,weights,score,finalize`.
- **Option B (score implies finalize):** When user requests `--steps score`, the pipeline runs both the score group (calc-score) and the finalize group (combine-scores, finalize-output). For **--sbatch** only, we still submit finalize as a **separate job** with slurm.finalize, so the driver never runs finalize in-process. So CLI stays “score” but SLURM gains a dedicated finalize job.

Recommendation: **Option B** — keep `--steps score` meaning “score + finalize” for CLI/docs; implement the split only in the SLURM path (submit finalize as a separate job with slurm.finalize after the score array). That way we fix the OOM and resource story without changing user-facing step names.

### 5. Summary

| Aspect | Current | Proposed |
|--------|--------|----------|
| score group | calc-score, combine-scores, finalize-output | calc-score only (finalize = separate group) |
| finalize group | (none) | combine-scores, finalize-output |
| --sbatch after score array | Driver runs combine-scores + finalize-output (can OOM) | Driver submits one **finalize job** (slurm.finalize); driver waits. |
| Config | slurm.score only | slurm.score + **slurm.finalize** (e.g. mem: 16g) |

**Status:** Finalize step group, slurm.finalize, and driver finalize job are implemented. CLI uses explicit `finalize` in `--steps`.

### What needs to be done (finalize and outputs)

1. **Plan: All joins via Unix join** — **Implemented.** Every join in finalize uses `LC_ALL=C sort` + `LC_ALL=C join` only; no in-memory hash lookups for join keys.
   - **Combine posteriors:** One `LC_ALL=C sort` + `LC_ALL=C join` on RSID (ldref_to_genoid).
   - **sumstat_augmented.tsv.gz:** **Approach B** — concatenate per-chr matched files from filter-variants (already mapfile-restricted, with LDREF_SNPID + GENO_ID), one `LC_ALL=C sort` + `LC_ALL=C join` with posteriors; one awk reorders columns and adds IN_ANALYSIS. 
   - **augmented_sumstat.gz:** Three Unix joins with LC_ALL=C (variant_map ⋈ sumstat extract, then ⋈ posteriors, then ⋈ MAF); one awk emits final 10 columns. No in-memory MAF/postEffect lookups.

2. **benchEffect in augmented_sumstat.gz** — When the benchmark step has been run, fill the benchEffect column from benchmark outputs (effect used in benchmark score per variant); currently output as `NA`.

## SLURM submission system (`--sbatch`)

### Overview
The `--sbatch` flag provides integrated SLURM job submission. It handles two distinct
workflows: **prep jobs** and **per-sumstat driver jobs**.

### Submission modes

1. **Prep job** (`--sbatch --steps prep`):
   - Submits a lightweight **driver job** for prep.
   - The driver runs `prep-ldref` (including liftover augmentation if not cached),
     then `prep-genotypes`, then submits a **prep array** (one task per chromosome)
     for `prep-inclusion-list`.
   - After the array completes, the driver combines per-chromosome maps and
     creates the final inclusion list.
   - Resources configured via `slurm.prep: { mem, cpus, time, max_parallel }`.
   - Logs written to `<outdir>/prep/logs/slurm/`.

2. **Per-sumstat driver job** (`--sbatch --steps sumstat,weights,score -i <sumstat>`):
   - Submits a lightweight **driver job** that orchestrates the full per-sumstat workflow.
   - The driver job runs `sumstat` directly (not chromosome-parallel).
   - For `weights`, the driver submits **two** SLURM array jobs (sBayesR weights and benchmark weights),
     each independently configurable. For `score`, it submits one array. The driver waits for each array
     to complete before proceeding.
   - Resources for the driver configured via `slurm.driver: { mem, cpus, time }`.
   - Logs written to `<outdir>/sumstats/<sumstat_name>/logs/slurm/`.

### Driver job architecture

```
User submission                     SLURM cluster
     │
     ▼
./pgscalculator-v2.sh --sbatch --steps sumstat,weights,score -i sumstat_814
     │
     └──► Driver job (pgs_sumstat_814_driver) ──────────────────────────────────►
              │
              ├── [1] Runs format-sumstat directly (chr split only, no liftover needed)
              │
              ├── [2] Submits sumstat array (1-22%max_parallel)
              │       └── Per-chr: mapfile join, EAF fill, filter, reduce
              │       └── Waits for all tasks to complete
              │
              ├── [3] Submits weights_sbayesr array (1-22%max_parallel)
              │       └── Per-chr: calc-posteriors, format-posteriors
              │       └── Waits for all tasks to complete
              │
              ├── [4] Submits weights_benchmark array (1-22%max_parallel)
              │       └── Per-chr: calc-benchmark
              │       └── Waits for all tasks to complete
              │
              ├── [5] Submits score array (1-N%max_parallel)
              │       └── Waits for all tasks to complete
              │
              └── [6] Runs combine-scores, finalize-output directly
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

  # Lightweight driver job (submits arrays for sumstat/weights/score)
  driver: { mem: 1g, cpus: 1, time: '02:00:00' }

  # Prep driver + array (sumstat-agnostic)
  prep: { mem: 2g, cpus: 1, time: '02:00:00', max_parallel: 22 }

  # Sumstat array (one task per chromosome, after format-sumstat splits)
  sumstat: { mem: 1g, cpus: 1, time: '00:30:00', max_parallel: 22 }

  # Weights: two array jobs, independently configurable
  weights_sbayesr:  { mem: 20g, cpus: 6, time: '01:00:00', max_parallel: 22 }
  weights_benchmark: { mem: 2g, cpus: 2, time: '00:30:00', max_parallel: 22 }

  # Score array (one task per chromosome with mapped posteriors)
  score: { mem: 10g, cpus: 4, time: '01:00:00', max_parallel: 22 }

# Benchmark calculation (MAF and LD pruning)
benchmark:
  maf_threshold: 0.05
  indep_pairwise: [250, 50, 0.25]   # window_kb, step, r2 for plink --indep-pairwise
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
3. **Logs warning**: `"Warning: SLURM array job X for step 'weights_sbayesr' (or 'weights_benchmark') had failed/unknown task(s)."`.
4. **Continues**: Does not abort; downstream arrays are submitted.
5. **Sanity-checks outputs**: After each array, verifies expected files exist; logs warning if missing.

#### Failure reasons tracked in FAILED markers

- `filtered_sumstat_missing`: No filtered sumstat input for this chromosome.
- `format_for_sbayesr_failed`: Could not generate sbayesR input file.
- `ma_missing`: sbayesR input file was not created.
- `ldref_missing`: LD reference files not found for chromosome.
- `sbayesr_failed`: gctb/sbayesR exited non-zero.
- `sbayesr_no_output`: gctb ran but produced no `.snpRes` output.
- `posteriors_missing`: Posteriors file missing for format-posteriors step.
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
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score -i /path/to/sumstat_814 --sbatch
./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score -i /path/to/sumstat_815 --sbatch
# ... etc
```

---

## Criteria checklist
- Mapfile uses a single `chr` with **dual position columns** (`pos_b37`, `pos_b38`) for build support.
- Mapfile contains all three source SNP IDs and their allele columns.
- Prep step is sumstat-agnostic; mapfile contains all LD ref liftover variants (genotype columns NA where no match).
- Prep step supports **GRCh38 genotypes** via pre-augmented LD reference (liftover done once in `prep-ldref`).
- Matching is always by **chr:pos + alleles** (using appropriate position column based on genotype build).
- SNP inclusion lists are always derived from the mapfile (prep or sumstat mapfile).
- Sumstat step attaches `sumstat_*` via `chr/pos_b38 + alleles` (native GRCh38 coordinates).
- Sumstat `EAF` is filled only from LD reference EAF (`ldref_a2freq`) when missing;
  `EAF_1KG` is not used.
- INFO/MAF reference files are replaced by a user-provided inclusion list, with an
  explicit ID-space specifier (`ss`, `ld`, or `gt`), applied after mapfile reduction.
- Conversions to/from LD reference and genotype IDs only use the mapfile.
- Final output includes the full mapfile (`variant_map.tsv.gz`) with both position columns.
- Output files (v2): augmented_sumstat.gz has **same row set as variant map**, RSID + sumstat effect columns + MAF + postEffect (+ benchEffect when added); main_raw_score_all.gz has IID and score columns (no FID); variant_map.tsv.gz has RSID (liftover-derived) as col1 and **all LD ref liftover variants**, with genotype/sumstat columns NA where no match.
- Scoring uses posteriors file only (no separate inclusion list for scoring).
- Benchmark: LD pruning + MAF filter + plink --score with observed effects; benchmark weights/scores as in v1-style.