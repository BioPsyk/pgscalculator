#!/usr/bin/env Rscript
# Compare sample-level PGS from scores_sbayesr.gz and scores_ldpred2.gz (§11 sanity check).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: compare_methods.R <scores_sbayesr.gz> <scores_ldpred2.gz> [min_r=0.6]")
}

read_scores <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  d <- read.delim(con, check.names = FALSE)
  iid_col <- grep("^IID$|^#IID$", names(d), value = TRUE)[1]
  score_col <- grep("^SCORE_SUM$|^SCORE1_SUM$", names(d), value = TRUE)[1]
  if (is.na(iid_col) || is.na(score_col)) {
    stop("Could not find IID and SCORE_SUM in ", path)
  }
  setNames(d[[score_col]], d[[iid_col]])
}

min_r <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.6
s1 <- read_scores(args[[1]])
s2 <- read_scores(args[[2]])
common <- intersect(names(s1), names(s2))
if (length(common) < 2) {
  stop("Fewer than 2 overlapping samples between score files")
}

r <- cor(s1[common], s2[common], use = "complete.obs")
cat(sprintf("Pearson r (n=%d samples): %.4f\n", length(common), r))
if (is.na(r) || r < min_r) {
  stop(sprintf("Correlation %.4f is below threshold %.2f", r, min_r))
}
cat("OK: correlation above threshold\n")
