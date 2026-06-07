# Reference data for pgscalculator v2

Where each reference data set comes from and how to install it.

All commands assume you are running from the repo root unless noted.
Reference paths in `config.yaml` can point anywhere on your filesystem; the
layout shown here (`references/<set>/`) matches the examples in
`config.template.yaml` but is purely a convention. The `references/`
directory is `.gitignore`d so binary data is never committed.

## Overview

| Reference set                       | Required for                              | Source                | Size on disk |
|-------------------------------------|-------------------------------------------|-----------------------|--------------|
| Liftover (dbSNP)                    | `prep-ldref` (dual-position variant map)  | cleansumstats         | ~21 GB       |
| cleansumstat RSID map               | `prep-inclusion-list`                     | cleansumstats         | ~14 GB       |
| sBayesR LD (UKB 10k HM3)            | `calc-posteriors` (sBayesR)               | cnsgenomics           | ~few GB      |
| LDpred2 LD — HapMap3+ (default)     | `calc-ldpred2` (LDpred2)                  | figshare + Google Drive | ~15 GB     |
| LDpred2 LD — HapMap3 (alternative)  | `calc-ldpred2` (LDpred2)                  | figshare              | ~7 GB        |

Sizes are uncompressed, on-disk. The LDpred2 sets are independent of each
other — you only need the one you configured.

## 1. Liftover reference (required)

Used by `prep-ldref` to attach both GRCh37 and GRCh38 positions to every LD
reference variant, so the rest of the pipeline is genome-build-agnostic at
variant-map time.

These files come from the upstream cleansumstats pipeline (a BioPsyk tool)
and are not public downloads. Copy them in (~10 GB each):

```bash
mkdir -p references/liftover
gzip -c <cleansumstats>/All_20180418_GRCh37_GRCh38.sorted.bed \
    > references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz
gzip -c <cleansumstats>/All_20180418_GRCh38_GRCh37.sorted.bed \
    > references/liftover/dbsnp_cleansumstat_reference_GRCh38_GRCh37.txt.gz
```

Point your `config.yaml` at the file matching the build of your genotypes:

```yaml
genotype_build: GRCh37
liftover_reference: /path/to/references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz
```

## 2. cleansumstat RSID map (required for variant inclusion)

Used by `prep-inclusion-list` to map between dbSNP rsids in the LD reference
and sumstat IDs. Built from the same source files as the liftover reference
plus the sBayesR LD reference variant list. The full build procedure is
documented in `docs/outdated/setup_references.md` (the awk-based extraction
loop) and is run once per project.

## 3. sBayesR LD reference (UKB 10k, HapMap3, sparse)

Source: cnsgenomics (the GCTB authors).

```bash
mkdir -p references/ld-sbayesr/ukb
cd references/ld-sbayesr/ukb
wget https://cnsgenomics.com/data/GCTB/band_ukb_10k_hm3.zip
unzip band_ukb_10k_hm3.zip
rm band_ukb_10k_hm3.zip
```

After unzip you should see
`references/ld-sbayesr/ukb/band_ukb_10k_hm3/band_chr{1..22}.ldm.sparse{,.info}`.

In `config.yaml`:

```yaml
lddir: /path/to/references/ld-sbayesr/ukb/band_ukb_10k_hm3
```

## 4. LDpred2 LD reference

LDpred2 supports two precomputed European LD reference sets from the
bigsnpr authors (Privé et al.). Choose explicitly via the
`ldpred2.ld_variant_set` config key:

- `hm3_plus` (default, recommended) — 1,444,196 HapMap3+ variants
- `hm3` (alternative) — 1,054,330 HapMap3 variants

The pipeline does **not** silently substitute one for the other. If the
configured set is missing or unloadable the run errors out, and the
operator must explicitly opt into the alternative.

### Recommended: HapMap3+ (1,444,196 variants, ~15 GB)

Used for new analyses per the 2023 LDpred2-auto paper (Privé et al.).
Larger variant coverage → better PGS power for well-powered GWAS.

```bash
mkdir -p references/ld-ldpred2/hm3_plus
cd references/ld-ldpred2/hm3_plus

# 1. Variant metadata (44 MB, figshare)
curl -L -o map_hm3_plus.rds \
    https://ndownloader.figshare.com/files/37802721

# 2. 22 chromosome LD blocks (~15 GB, Google Drive).
#    The figshare article for HM3+ only carries the map; the LD matrices
#    live on G-drive. gdown handles the interstitial confirmation flow.
#    Install once:  pip install --user gdown  (or via conda)
gdown --id 17dyKGA2PZjMsivlYb_AjDmuZjM1RIGvs
unzip ldref_hm3_plus.zip
rm ldref_hm3_plus.zip
```

After unzip you should see 22 × `LD_with_blocks_chr{N}.rds` plus
`map_hm3_plus.rds` directly under `references/ld-ldpred2/hm3_plus/`.

`config.yaml`:

```yaml
ldpred2:
  ld_variant_set: hm3_plus
  ld_dir: /path/to/references/ld-ldpred2/hm3_plus
```

### Alternative: HapMap3 (1,054,330 variants, ~7 GB)

Use for cross-validation against legacy LDpred2 studies, when HM3+ is
unavailable, or when the variant set you care about has poor HM3+
coverage. Single direct figshare URL — no Google Drive flow.

```bash
mkdir -p references/ld-ldpred2/hm3
cd references/ld-ldpred2/hm3

# 1. Variant metadata (35 MB, figshare)
curl -L -o map_hm3.rds \
    https://ndownloader.figshare.com/files/36360900

# 2. 22 chromosome LD blocks bundled in one zip (7.7 GB, figshare)
curl -L -o ldref_with_blocks.zip \
    https://ndownloader.figshare.com/files/36363087
unzip ldref_with_blocks.zip
# The zip extracts into an ldref/ subdirectory; flatten so all files
# sit directly under references/ld-ldpred2/hm3/ (matches the hm3_plus
# layout the pipeline expects):
mv ldref/LD_with_blocks_chr*.rds .
rm -f ldref/map.rds   # duplicate of map_hm3.rds we just downloaded
rmdir ldref
rm ldref_with_blocks.zip
```

`config.yaml`:

```yaml
ldpred2:
  ld_variant_set: hm3
  ld_dir: /path/to/references/ld-ldpred2/hm3
```

### Verifying the LDpred2 download

A quick check in the container — change the path for the set you
downloaded:

```bash
singularity exec sif/ibp-pgscalculator-base_version-X.Y.Z.sif Rscript -e '
suppressPackageStartupMessages(library(bigsnpr))
dir <- "references/ld-ldpred2/hm3_plus"  # or .../hm3
map_name <- if (grepl("hm3_plus$", dir)) "map_hm3_plus.rds" else "map_hm3.rds"
m  <- readRDS(file.path(dir, map_name))
LD <- readRDS(file.path(dir, "LD_with_blocks_chr22.rds"))
stopifnot(inherits(LD, "dsCMatrix"))
cat("OK:", nrow(m), "variants,", nrow(LD), "x", ncol(LD), "chr22 LD\n")
'
```

Expected `nrow(m)` is 1,444,196 (HM3+) or 1,054,330 (HM3); expected
chr22 LD is 21233 × 21233 (HM3+) or 15449 × 15449 (HM3).

#### Schema note (relevant only if you read the map outside the pipeline)

The two map files share most columns but have minor differences:

| Column        | HM3+ (`map_hm3_plus.rds`) | HM3 (`map_hm3.rds`) |
|---------------|---------------------------|----------------------|
| LD-block id   | `block_id`                | `group_id`           |
| hg17 position | not present               | `pos_hg17`           |
| hg18 position | `pos_hg18`                | `pos_hg18`           |
| hg38 position | `pos_hg38`                | `pos_hg38`           |

The pipeline's LDpred2 R script aliases `group_id → block_id` on read so
both sets work without further config.

## 5. Test data

For the smoke tests and example runs:

- `tests/example_data/sumstats/` — committed, small test sumstats.
- `references/genotypes_test/` — example genotypes, *not* in git. Sourced
  from the BioPsyk project storage; see the team's onboarding notes.

## 6. Legacy / not currently used by v2

The repo retains setup notes and a few downloaded directories for tools
that are not invoked by the v2 driver but were used by the v1 Nextflow
pipeline and may be revived later:

- **PRS-CS** LD references (`references/ld-prscs/`) — sourced from the
  PRS-CS Dropbox links documented in
  [`docs/outdated/setup_references.md`](outdated/setup_references.md).
- **SBayesRC** LD matrices (`references/sbayesrc/`) — downloaded via
  `references/sbayesrc/download_sbayesrc_data.sh` (HapMap3 from
  `gctbhub.cloud.edu.au`, ~3 GB; plus a 54 GB imputed-variant set).

These references can stay on disk but are not consumed by the v2
codepath today. When/if PRS-CS or SBayesRC is integrated into v2, the
relevant section will be promoted into this document.
