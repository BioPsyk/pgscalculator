#!/usr/bin/env Rscript
# Quick correlation plot of per-individual PGS across methods + benchmarks.
suppressWarnings(suppressMessages({}))

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
png_out <- args[2]

read_score <- function(path, col) {
  df <- read.table(gzfile(path), header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  data.frame(IID = df$IID, score = df[[col]], stringsAsFactors = FALSE)
}

s1 <- read_score(file.path(outdir, "scores_sbayesr.gz"),       "SCORE_SUM")
s2 <- read_score(file.path(outdir, "scores_ldpred2.gz"),       "SCORE_SUM")
b1 <- read_score(file.path(outdir, "bench_score_sbayesr.gz"),  "SCORE1_SUM")
b2 <- read_score(file.path(outdir, "bench_score_ldpred2.gz"),  "SCORE1_SUM")

m <- Reduce(function(a, b) merge(a, b, by = "IID"),
            list(setNames(s1, c("IID", "sBayesR")),
                 setNames(s2, c("IID", "LDpred2")),
                 setNames(b1, c("IID", "bench_sBayesR")),
                 setNames(b2, c("IID", "bench_LDpred2"))))

mat <- as.matrix(m[, c("sBayesR", "LDpred2", "bench_sBayesR", "bench_LDpred2")])
cat(sprintf("N individuals: %d\n", nrow(mat)))
cr <- cor(mat, method = "pearson")
cat("\nPearson correlation matrix:\n")
print(round(cr, 4))
crs <- cor(mat, method = "spearman")
cat("\nSpearman correlation matrix:\n")
print(round(crs, 4))

# Scatter-matrix with Pearson r in upper panels, histograms on the diagonal.
panel_cor <- function(x, y, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- sprintf("r = %.3f", r)
  cex <- 1.2 + abs(r) * 1.6
  text(0.5, 0.5, txt, cex = cex, col = ifelse(r >= 0, "#1f6feb", "#d1242f"))
}
panel_hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  br <- h$breaks; nB <- length(br)
  y <- h$counts; y <- y / max(y)
  rect(br[-nB], 0, br[-1], y, col = "#c8d6f0", border = "white")
}
panel_pts <- function(x, y, ...) {
  points(x, y, pch = 16, cex = 0.35, col = adjustcolor("#333333", 0.35))
  abline(lm(y ~ x), col = "#d1242f", lwd = 1.5)
}

png(png_out, width = 1500, height = 1500, res = 150)
pairs(mat,
      labels = c("sBayesR", "LDpred2", "bench\nsBayesR", "bench\nLDpred2"),
      lower.panel = panel_pts, upper.panel = panel_cor, diag.panel = panel_hist,
      main = "Per-individual PGS correlation across methods (genome-wide, n=2504)")
invisible(dev.off())
cat(sprintf("\nWrote plot: %s\n", png_out))
