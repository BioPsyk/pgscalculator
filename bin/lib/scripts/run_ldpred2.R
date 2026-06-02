#!/usr/bin/env Rscript
# Genome-wide LDpred2-auto / LDpred2-inf posterior estimation.
#
# Reads per-chr filtered sumstats from work/filtered_ldpred2/, matches against
# the precomputed HM3/HM3+ LD reference, runs LDSC + LDpred2, and writes
# per-chr .snpRes files under work/posteriors_ldpred2/.
#
# Invoked by bin/lib/steps/calc_ldpred2.sh (see docs/plans/ldpred2-integration.md §5.5).

suppressPackageStartupMessages({
    library(argparser)
    library(bigsnpr)
    library(bigreadr)
    library(bigsparser)
    library(data.table)
    library(ggplot2)
})

options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

fail_with <- function(marker_file, msg) {
    cat("[ERROR]", msg, "\n", sep = "")
    if (nzchar(marker_file)) {
        writeLines(msg, marker_file)
    }
    quit(status = 1, save = "no")
}

write_empty_snpres_dir <- function(out_dir) {
    hdr <- paste(
        "Id", "Name", "Chrom", "Position", "A1", "A2",
        "A1Frq", "A1Effect", "SE", "PIP", "LastSampleEff"
    )
    for (chr in 1:22) {
        f <- file.path(out_dir, sprintf("chr%d.snpRes", chr))
        writeLines(c(hdr), f)
    }
}

write_chr_snpres <- function(out_dir, df_chr, beta_chr, postp_chr) {
    chr <- unique(df_chr$chr)[1]
    out_file <- file.path(out_dir, sprintf("chr%d.snpRes", chr))
    hdr <- paste(
        "Id", "Name", "Chrom", "Position", "A1", "A2",
        "A1Frq", "A1Effect", "SE", "PIP", "LastSampleEff"
    )
    lines <- vapply(seq_len(nrow(df_chr)), function(i) {
        se_val <- "NA"
        pip_val <- if (is.na(postp_chr[i])) "NA" else sprintf("%.8f", postp_chr[i])
        sprintf(
            "%d %s %d %d %s %s %.6f %.6f %s %s %.6f",
            i,
            df_chr$rsid[i],
            df_chr$chr[i],
            df_chr$pos[i],
            df_chr$a1[i],
            df_chr$a0[i],
            df_chr$af_UKBB[i],
            beta_chr[i],
            se_val,
            pip_val,
            beta_chr[i]
        )
    }, character(1))
    writeLines(c(hdr, lines), out_file)
}

get_betas_auto <- function(multi_auto) {
    keep <- which(vapply(multi_auto, function(auto) {
        all(auto$path_h2_est > 0, na.rm = TRUE) &&
            all(auto$path_p_est < 0.5, na.rm = TRUE) &&
            all(auto$path_p_est > 1e-4, na.rm = TRUE)
    }, logical(1)))
    if (length(keep) == 0) {
        return(list(beta = NULL, postp = NULL, keep = keep))
    }
    beta_mat <- sapply(multi_auto[keep], function(auto) auto$beta_est)
    postp_mat <- sapply(multi_auto[keep], function(auto) auto$postp_est)
    list(
        beta = rowMeans(beta_mat),
        postp = rowMeans(postp_mat),
        keep = keep
    )
}

p <- arg_parser("Run genome-wide LDpred2 on filtered summary statistics")
p <- add_argument(p, "--sumstat-dir", type = "character",
                  help = "Directory with chr{N}_filtered.tsv files")
p <- add_argument(p, "--ld-dir", type = "character",
                  help = "Directory with LD_with_blocks_chr{N}.rds")
p <- add_argument(p, "--ld-meta", type = "character",
                  help = "Path to map_hm3.rds or map_hm3_plus.rds")
p <- add_argument(p, "--out-dir", type = "character",
                  help = "Output directory for chr{N}.snpRes files")
p <- add_argument(p, "--mode", type = "character", default = "auto",
                  help = "LDpred2 mode: auto | inf")
p <- add_argument(p, "--shrink-corr", type = "numeric", default = 0.95,
                  help = "Correlation shrinkage before LDpred2")
p <- add_argument(p, "--allow-jump-sign", type = "logical", default = FALSE,
                  help = "Allow sign flips in LDpred2-auto")
p <- add_argument(p, "--hyper-p-max", type = "numeric", default = 0.2,
                  help = "Max grid value for p (auto mode)")
p <- add_argument(p, "--hyper-p-length", type = "integer", default = 30L,
                  help = "Grid length for p (auto mode)")
p <- add_argument(p, "--seed", type = "integer", default = 1L,
                  help = "Random seed")
p <- add_argument(p, "--ncores", type = "integer", default = 1L,
                  help = "Number of CPU cores")
p <- add_argument(p, "--genotype-build", type = "character", default = "GRCh37",
                  help = "Deprecated: genome build of the genotype files (no longer used for LD matching)")
p <- add_argument(p, "--ld-build", type = "character", default = "GRCh37",
                  help = "Genome build of the LD reference's native positions (GRCh37 or GRCh38)")
p <- add_argument(p, "--merge-by-rsid", type = "logical", default = FALSE,
                  help = "Match sumstat to LD map by RSID instead of chr:pos")
p <- add_argument(p, "--effective-sample-size", type = "numeric", default = NA,
                  help = "Fallback effective sample size")
p <- add_argument(p, "--n-cases", type = "numeric", default = NA,
                  help = "Case count for h2 scaling")
p <- add_argument(p, "--n-controls", type = "numeric", default = NA,
                  help = "Control count for h2 scaling")
p <- add_argument(p, "--plot-file", type = "character", default = "",
                  help = "Optional path for LDpred2-auto diagnostic plot (chains.png)")
p <- add_argument(p, "--summary-file", type = "character", default = "",
                  help = "Optional path for the LDpred2 run diagnostics table (summary.tsv)")
argv <- parse_args(p)

# Diagnostics collected through the run and flushed to --summary-file at the end
# (long format: metric<TAB>value). Written even on a successful partial run so
# the user can sanity-check h2 / p / variant counts.
diag_env <- new.env(parent = emptyenv())
diag_env$rows <- list()
add_diag <- function(metric, value) {
    if (is.null(value) || length(value) == 0) value <- NA
    diag_env$rows[[length(diag_env$rows) + 1L]] <- c(metric, as.character(value))
}
write_summary <- function() {
    path <- argv$summary_file
    if (is.null(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
        return(invisible(FALSE))
    }
    if (length(diag_env$rows) == 0) return(invisible(FALSE))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    con <- file(path, "w")
    on.exit(close(con), add = TRUE)
    writeLines("metric\tvalue", con)
    for (r in diag_env$rows) writeLines(paste(r[1], r[2], sep = "\t"), con)
    cat(sprintf("[OK] Wrote run diagnostics: %s\n", path))
    invisible(TRUE)
}
fmt_num <- function(x, digits = 6) {
    if (is.null(x) || length(x) != 1 || is.na(x) || !is.finite(x)) return("NA")
    trimws(formatC(as.numeric(x), format = "g", digits = digits))
}

required_args <- c("sumstat_dir", "ld_dir", "ld_meta", "out_dir")
for (arg_name in required_args) {
    if (is.na(argv[[arg_name]]) || !nzchar(argv[[arg_name]])) {
        stop(sprintf("--%s is required", gsub("_", "-", arg_name)))
    }
}

if (!argv$mode %in% c("auto", "inf")) {
    stop("--mode must be 'auto' or 'inf'")
}

# Clamp requested cores to what is actually available to the process. bigsnpr's
# assert_cores() aborts (rather than degrading) if ncores exceeds the cgroup/OS
# core count, so a threads/cpus mismatch would otherwise fail a long-running job.
ncores <- argv$ncores
avail_cores <- tryCatch(bigparallelr::nb_cores(), error = function(e) 1L)
if (is.na(ncores) || ncores < 1L) ncores <- 1L
if (ncores > avail_cores) {
    cat(sprintf("[OK] Requested %d cores but only %d available; using %d\n",
                ncores, avail_cores, avail_cores))
    ncores <- avail_cores
}

sumstat_dir <- argv$sumstat_dir
ld_dir <- argv$ld_dir
ld_meta <- argv$ld_meta
out_dir <- argv$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

add_diag("mode", argv$mode)
add_diag("seed", argv$seed)
add_diag("ncores_used", ncores)
add_diag("shrink_corr", fmt_num(argv$shrink_corr))
add_diag("ld_meta", basename(ld_meta))

failed_match <- file.path(out_dir, "FAILED_match")
failed_sfbm <- file.path(out_dir, "FAILED_sfbm")
failed_ldsc <- file.path(out_dir, "FAILED_ldsc")
failed_ldpred2 <- file.path(out_dir, "FAILED_ldpred2")
for (f in c(failed_match, failed_sfbm, failed_ldsc, failed_ldpred2)) {
    if (file.exists(f)) unlink(f)
}

# --- Load and stack per-chr filtered files ------------------------------------
parts <- list()
for (chr in 1:22) {
    f <- file.path(sumstat_dir, sprintf("chr%d_filtered.tsv", chr))
    if (!file.exists(f)) next
    dt <- fread(f, showProgress = FALSE)
    if (nrow(dt) == 0) next
    parts[[length(parts) + 1L]] <- dt
}

if (length(parts) == 0) {
    fail_with(failed_match, "No chr*_filtered.tsv files found in sumstat-dir")
}

raw <- rbindlist(parts, use.names = TRUE, fill = TRUE)
need_cols <- c("CHR", "POS", "EffectAllele", "OtherAllele", "B", "SE", "N")
missing_cols <- setdiff(need_cols, names(raw))
if (length(missing_cols) > 0) {
    fail_with(
        failed_match,
        sprintf("Filtered sumstat missing column(s): %s", paste(missing_cols, collapse = ", "))
    )
}

rsid_col <- if ("LDREF_SNPID" %in% names(raw)) "LDREF_SNPID" else if ("RSID" %in% names(raw)) "RSID" else NA
if (is.na(rsid_col)) {
    fail_with(failed_match, "Filtered sumstat needs LDREF_SNPID or RSID column")
}

sumstats <- data.table(
    chr = as.integer(raw$CHR),
    pos = as.integer(raw$POS),
    a1 = toupper(as.character(raw$EffectAllele)),
    a0 = toupper(as.character(raw$OtherAllele)),
    rsid = as.character(raw[[rsid_col]]),
    beta = as.numeric(raw$B),
    beta_se = as.numeric(raw$SE),
    n_eff = as.numeric(raw$N)
)

sumstats <- sumstats[
    !is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(beta_se) & beta_se > 0
]

eff_n <- argv$effective_sample_size
n_cases <- argv$n_cases
n_controls <- argv$n_controls
if (length(eff_n) == 1 && !is.na(eff_n)) {
    sumstats[, n_eff := eff_n]
} else if (length(n_cases) == 1 && length(n_controls) == 1 &&
           !is.na(n_cases) && !is.na(n_controls)) {
    n_eff_derived <- 4 / (1 / n_cases + 1 / n_controls)
    sumstats[, n_eff := n_eff_derived]
}

sumstats <- sumstats[!is.na(n_eff) & n_eff > 0]

if (nrow(sumstats) < 100) {
    fail_with(failed_match, sprintf("Too few valid sumstat rows after QC: %d", nrow(sumstats)))
}

# --- LD reference map ---------------------------------------------------------
# The filtered sumstat positions are GRCh38 (cleansumstats standard; filter-variants
# matches the formatted sumstat on pos_b38). So we match against whichever LD-map
# column is in GRCh38. `ld_build` declares the build of the LD map's native `pos`:
#   - ld_build GRCh38 -> native `pos` is already GRCh38, match on it directly
#   - ld_build GRCh37 -> native `pos` is GRCh37, so match on the `pos_hg38` column
# This keeps matching positional (via the dbSNP-reconciled variant map) and never
# relies on the LD reference's own rsids, which may diverge from the dbSNP backbone.
map_ldref <- readRDS(ld_meta)
ld_build <- toupper(if (length(argv$ld_build) == 1 && !is.na(argv$ld_build)) {
    argv$ld_build
} else {
    "GRCH37"
})
if (ld_build %in% c("GRCH38", "HG38")) {
    if (!("pos" %in% colnames(map_ldref))) {
        fail_with(failed_match, "LD map missing pos column")
    }
} else {
    if (!("pos_hg38" %in% colnames(map_ldref))) {
        fail_with(failed_match,
                  "LD map missing pos_hg38 (needed to match a GRCh38 sumstat against a GRCh37-native LD reference)")
    }
    map_ldref$pos <- map_ldref$pos_hg38
}

if (!("block_id" %in% colnames(map_ldref)) && ("group_id" %in% colnames(map_ldref))) {
    map_ldref$block_id <- map_ldref$group_id
}

join_by_pos <- !isTRUE(argv$merge_by_rsid)
add_diag("ld_build", ld_build)
add_diag("match_by", if (join_by_pos) "chr:pos" else "rsid")
add_diag("n_sumstat_input", nrow(sumstats))
cat(sprintf("[OK] Matching %d sumstat rows to LD reference (join_by_pos=%s)\n",
            nrow(sumstats), join_by_pos))

df_beta <- tryCatch(
    snp_match(sumstats, map_ldref, join_by_pos = join_by_pos, match.min.prop = 0),
    error = function(e) {
        fail_with(failed_match, conditionMessage(e))
        NULL
    }
)

# Keep `_NUM_ID_` (row index into the LD map) — it is required below to subset the
# per-chromosome LD matrices. Only drop the duplicated sumstat-side columns.
drops <- c("_NUM_ID_.ss", "rsid.ss")
df_beta <- df_beta[, !(names(df_beta) %in% drops), drop = FALSE]

if (nrow(df_beta) < 100) {
    fail_with(
        failed_match,
        sprintf("snp_match retained too few variants: %d", nrow(df_beta))
    )
}
cat(sprintf("[OK] snp_match retained %d variants\n", nrow(df_beta)))
add_diag("n_matched", nrow(df_beta))

# --- Allele-frequency QC (bigsnpr vignette) -----------------------------------
sd_ldref <- sqrt(2 * df_beta$af_UKBB * (1 - df_beta$af_UKBB))
sd_ss <- 2 / sqrt(df_beta$n_eff * df_beta$beta_se^2 + df_beta$beta^2)
is_bad <- sd_ss < 0.5 * sd_ldref |
    sd_ss > sd_ldref + 0.1 |
    sd_ss < 0.05 |
    sd_ldref < 0.05
n_bad <- sum(is_bad, na.rm = TRUE)
if (n_bad > 0) {
    cat(sprintf("[OK] QC: removing %d variants with inconsistent allelic SD\n", n_bad))
    df_beta <- df_beta[!is_bad, , drop = FALSE]
}

add_diag("n_qc_removed", n_bad)
if (nrow(df_beta) < 100) {
    fail_with(failed_match, "Too few variants remain after allelic SD QC")
}
add_diag("n_final", nrow(df_beta))

# --- Build genome-wide SFBM ---------------------------------------------------
# The SFBM backing file is memory-mapped by bigsnpr. Many shared/network
# filesystems (e.g. the /faststorage Lustre/GDK mounts) reject mmap with
# "Error when mapping file: Invalid argument", so place the temporary backing on
# node-local, mmap-capable storage (TMPDIR if set, else /tmp). It is removed at
# exit; only the final chr*.snpRes outputs are written to out_dir.
# NB: TMPDIR may be set-but-empty (e.g. exported as ""), in which case
# Sys.getenv() returns "" rather than the `unset` fallback. Guard against that and
# fall back to R's resolved session tempdir(), which is always valid and writable.
sfbm_tmpdir <- Sys.getenv("TMPDIR")
if (!nzchar(sfbm_tmpdir) || !dir.exists(sfbm_tmpdir)) {
    sfbm_tmpdir <- tempdir()
}
tmp_file <- tempfile(tmpdir = sfbm_tmpdir, pattern = "ldpred2_corr_")
on.exit(unlink(paste0(tmp_file, c("", ".sbk", ".rds")), force = TRUE), add = TRUE)
ld_size <- 0
corr <- NULL

for (chr in sort(unique(df_beta$chr))) {
    ind.chr <- which(df_beta$chr == chr)
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

    ld_file <- file.path(ld_dir, sprintf("LD_with_blocks_chr%d.rds", chr))
    if (!file.exists(ld_file)) {
        fail_with(failed_sfbm, sprintf("Missing LD file: %s", ld_file))
    }

    num_ldref_snps <- sum(map_ldref$chr == chr)
    ld_size <- ld_size + num_ldref_snps
    cat(sprintf("[OK] chr%d: loading LD for %d / %d reference SNPs\n",
                chr, length(ind.chr), num_ldref_snps))

    corr_chr <- tryCatch(
        readRDS(ld_file)[ind.chr3, ind.chr3],
        error = function(e) {
            fail_with(failed_sfbm, conditionMessage(e))
            NULL
        }
    )

    if (is.null(corr)) {
        corr <- tryCatch(
            as_SFBM(corr_chr, tmp_file, compact = TRUE),
            error = function(e) {
                fail_with(failed_sfbm, conditionMessage(e))
                NULL
            }
        )
    } else {
        tryCatch(
            corr$add_columns(corr_chr, nrow(corr)),
            error = function(e) {
                fail_with(failed_sfbm, conditionMessage(e))
            }
        )
    }
}

# --- LD score regression ------------------------------------------------------
cat("[OK] Running LD score regression\n")
ldsc <- tryCatch(
    with(df_beta, snp_ldsc(
        ld, ld_size,
        chi2 = (beta / beta_se)^2,
        sample_size = n_eff,
        blocks = NULL,
        ncores = ncores
    )),
    error = function(e) {
        fail_with(failed_ldsc, conditionMessage(e))
        NULL
    }
)

h2_est <- ldsc[["h2"]]
add_diag("ld_ref_size", ld_size)
add_diag("ldsc_intercept", fmt_num(ldsc[["int"]], 4))
add_diag("ldsc_h2", fmt_num(h2_est, 6))
if (!is.finite(h2_est) || is.na(h2_est)) {
    fail_with(failed_ldsc, sprintf("LDSC returned non-finite h2: %s", h2_est))
}
cat(sprintf("[OK] LDSC: intercept=%.4f h2=%.4f\n", ldsc[["int"]], h2_est))

# --- LDpred2 ------------------------------------------------------------------
postp_est <- rep(NA_real_, nrow(df_beta))
beta_est <- NULL

if (argv$mode == "inf") {
    cat("[OK] Running LDpred2-inf\n")
    beta_est <- tryCatch(
        snp_ldpred2_inf(corr, df_beta, h2 = h2_est),
        error = function(e) {
            fail_with(failed_ldpred2, conditionMessage(e))
            NULL
        }
    )
    add_diag("h2_est", fmt_num(h2_est, 6))
} else {
    cat("[OK] Running LDpred2-auto\n")
    set.seed(argv$seed)
    multi_auto <- tryCatch(
        snp_ldpred2_auto(
            corr, df_beta,
            h2_init = h2_est,
            vec_p_init = seq_log(1e-4, argv$hyper_p_max, length.out = argv$hyper_p_length),
            allow_jump_sign = argv$allow_jump_sign,
            shrink_corr = argv$shrink_corr,
            ncores = ncores
        ),
        error = function(e) {
            fail_with(failed_ldpred2, conditionMessage(e))
            NULL
        }
    )

    auto_res <- get_betas_auto(multi_auto)
    n_chains_total <- length(multi_auto)
    n_chains_kept <- length(auto_res$keep)
    add_diag("n_chains_total", n_chains_total)
    add_diag("n_chains_kept", n_chains_kept)
    if (n_chains_kept == 0) {
        fail_with(failed_ldpred2, "All LDpred2-auto chains diverged")
    }
    beta_est <- auto_res$beta
    postp_est <- auto_res$postp

    # Summarise hyper-parameters across the kept (stable) chains.
    kept <- multi_auto[auto_res$keep]
    chain_scalar <- function(field) {
        v <- vapply(kept, function(a) {
            x <- tryCatch(a[[field]], error = function(e) NA_real_)
            if (is.null(x) || length(x) != 1) NA_real_ else as.numeric(x)
        }, numeric(1))
        v[is.finite(v)]
    }
    p_vals <- chain_scalar("p_est")
    h2_vals <- chain_scalar("h2_est")
    alpha_vals <- chain_scalar("alpha_est")
    p_est_final <- if (length(p_vals)) median(p_vals) else NA_real_
    h2_est_final <- if (length(h2_vals)) median(h2_vals) else NA_real_
    alpha_est_final <- if (length(alpha_vals)) median(alpha_vals) else NA_real_
    add_diag("p_est", fmt_num(p_est_final, 6))
    add_diag("h2_est", fmt_num(h2_est_final, 6))
    add_diag("alpha_est", fmt_num(alpha_est_final, 6))
    cat(sprintf("[OK] LDpred2-auto kept %d/%d chains; p=%s h2=%s\n",
                n_chains_kept, n_chains_total, fmt_num(p_est_final), fmt_num(h2_est_final)))

    if (nzchar(argv$plot_file)) {
        # Overlay the post-burn-in p / h2 sampling paths of every kept chain so
        # the user can eyeball convergence/agreement across chains.
        dta <- do.call(rbind, lapply(seq_along(auto_res$keep), function(j) {
            a <- kept[[j]]
            data.frame(
                x = seq_along(a$path_p_est),
                path_p_est = a$path_p_est,
                path_h2_est = a$path_h2_est,
                chain = factor(auto_res$keep[j])
            )
        }))
        plt <- tryCatch(
            plot_grid(
                ggplot(dta, aes(y = path_p_est, x = x, color = chain)) +
                    geom_point(size = 0.6, show.legend = FALSE) +
                    theme_bigstatsr() +
                    geom_hline(yintercept = p_est_final, col = "blue", linetype = "dashed") +
                    scale_y_log10() +
                    labs(y = "p", x = "iteration"),
                ggplot(dta, aes(y = path_h2_est, x = x, color = chain)) +
                    geom_point(size = 0.6, show.legend = FALSE) +
                    theme_bigstatsr() +
                    geom_hline(yintercept = h2_est_final, col = "blue", linetype = "dashed") +
                    labs(y = "h2", x = "iteration"),
                ncol = 1,
                align = "hv"
            ),
            error = function(e) {
                cat(sprintf("[WARN] chain plot failed: %s\n", conditionMessage(e)))
                NULL
            }
        )
        if (!is.null(plt)) {
            ggsave(argv$plot_file, plt, width = 6, height = 8)
            cat(sprintf("[OK] Wrote chain diagnostics: %s\n", argv$plot_file))
        }
    }
}

if (is.null(beta_est) || length(beta_est) != nrow(df_beta)) {
    fail_with(failed_ldpred2, "LDpred2 returned an invalid beta vector")
}

# --- Write per-chr .snpRes (sbayesR-compatible layout for format-posteriors) --
for (chr in 1:22) {
    idx <- which(df_beta$chr == chr)
    if (length(idx) == 0) {
        writeLines(
            "Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP LastSampleEff",
            file.path(out_dir, sprintf("chr%d.snpRes", chr))
        )
        next
    }
    write_chr_snpres(
        out_dir,
        df_beta[idx, , drop = FALSE],
        beta_est[idx],
        postp_est[idx]
    )
}

add_diag("n_posteriors_written", nrow(df_beta))
write_summary()

cat(sprintf("[OK] Wrote posteriors to %s\n", out_dir))
invisible(TRUE)
