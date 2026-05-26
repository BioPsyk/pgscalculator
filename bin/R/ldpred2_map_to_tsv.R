#!/usr/bin/env Rscript
# Convert an LDpred2 LD-reference map (map_hm3.rds or map_hm3_plus.rds) to a
# TSV consumable by the rest of the pgscalculator pipeline, and validate the
# 22 per-chromosome LD-block files alongside it.
#
# Invoked by bin/lib/steps/prep_ldref_ldpred2.sh.
#
# Input:
#   --ld-dir   <dir>   containing LD_with_blocks_chr{1..22}.rds
#   --map-file <path>  e.g. .../map_hm3_plus.rds or .../map_hm3.rds
#   --out-tsv  <path>  e.g. <outdir>/prep/ldref_ldpred2/map.tsv
#
# Output schema (tab-separated, single header line):
#   chr  pos_b37  pos_b38  a0  a1  rsid  af_UKBB  ld  block_id
#
# The bigsnpr ecosystem ships HM3 maps with a `group_id` column and HM3+ maps
# with a `block_id` column; this script aliases the former into the latter so
# downstream code only needs to know one name.

suppressPackageStartupMessages({
    library(argparser)
    library(data.table)
})

p <- arg_parser(
    "Validate an LDpred2 LD reference directory and dump its map to TSV."
)
p <- add_argument(p, "--ld-dir",
                  help = "Directory containing LD_with_blocks_chr{N}.rds",
                  type = "character")
p <- add_argument(p, "--map-file",
                  help = "Path to map_hm3.rds or map_hm3_plus.rds",
                  type = "character")
p <- add_argument(p, "--out-tsv",
                  help = "Output TSV path",
                  type = "character")
argv <- parse_args(p)

if (!dir.exists(argv$ld_dir)) {
    stop(sprintf("ld-dir does not exist: %s", argv$ld_dir))
}
if (!file.exists(argv$map_file)) {
    stop(sprintf("map-file does not exist: %s", argv$map_file))
}

# 1) 22 LD-block files
missing <- character(0)
for (chr in 1:22) {
    f <- file.path(argv$ld_dir, sprintf("LD_with_blocks_chr%d.rds", chr))
    if (!file.exists(f)) missing <- c(missing, basename(f))
}
if (length(missing) > 0) {
    stop(sprintf("Missing %d LD-block file(s) in %s: %s",
                 length(missing),
                 argv$ld_dir,
                 paste(missing, collapse = ", ")))
}
cat(sprintf("[OK] 22 LD_with_blocks_chr*.rds present in %s\n", argv$ld_dir))

# 2) Map
m <- readRDS(argv$map_file)
cat(sprintf("[OK] map loaded: %d variants, %d cols (%s)\n",
            nrow(m),
            ncol(m),
            paste(colnames(m), collapse = ", ")))

# 3) HM3 ships `group_id`; HM3+ ships `block_id`. Normalise on the latter.
if (!("block_id" %in% colnames(m)) && ("group_id" %in% colnames(m))) {
    m$block_id <- m$group_id
    cat("[OK] aliased group_id -> block_id (HM3 schema)\n")
}

required <- c("chr", "pos", "a0", "a1", "rsid",
              "af_UKBB", "ld", "block_id", "pos_hg38")
absent <- setdiff(required, colnames(m))
if (length(absent) > 0) {
    stop(sprintf("Map is missing required column(s): %s",
                 paste(absent, collapse = ", ")))
}

# 4) Emit TSV. bigsnpr ships HM3 and HM3+ both with `pos` = GRCh37.
out <- data.table(
    chr      = m$chr,
    pos_b37  = m$pos,
    pos_b38  = m$pos_hg38,
    a0       = m$a0,
    a1       = m$a1,
    rsid     = m$rsid,
    af_UKBB  = m$af_UKBB,
    ld       = m$ld,
    block_id = m$block_id
)

out_dir <- dirname(argv$out_tsv)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
fwrite(out, argv$out_tsv, sep = "\t", quote = FALSE)
cat(sprintf("[OK] wrote %s (%d rows)\n", argv$out_tsv, nrow(out)))
