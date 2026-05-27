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
p <- add_argument(p, "--shrink-corr", type = "numeric", default = 0.95)
p <- add_argument(p, "--allow-jump-sign", type = "logical", default = FALSE)
p <- add_argument(p, "--hyper-p-max", type = "numeric", default = 0.2)
p <- add_argument(p, "--hyper-p-length", type = "integer", default = 30L)
p <- add_argument(p, "--seed", type = "integer", default = 1L)
p <- add_argument(p, "--ncores", type = "integer", default = 1L)
p <- add_argument(p, "--genotype-build", type = "character", default = "GRCh37",
                  help = "GRCh37 or GRCh38 (selects pos column for matching)")
p <- add_argument(p, "--merge-by-rsid", type = "logical", default = FALSE)
p <- add_argument(p, "--effective-sample-size", type = "numeric", default = NA)
p <- add_argument(p, "--n-cases", type = "numeric", default = NA)
p <- add_argument(p, "--n-controls", type = "numeric", default = NA)
p <- add_argument(p, "--plot-file", type = "character", default = "")
argv <- parse_args(p)

required_args <- c("sumstat_dir", "ld_dir", "ld_meta", "out_dir")
for (arg_name in required_args) {
    if (is.na(argv[[arg_name]]) || !nzchar(argv[[arg_name]])) {
        stop(sprintf("--%s is required", gsub("_", "-", arg_name)))
    }
}

if (!argv$mode %in% c("auto", "inf")) {
    stop("--mode must be 'auto' or 'inf'")
}

sumstat_dir <- argv$sumstat_dir
ld_dir <- argv$ld_dir
ld_meta <- argv$ld_meta
out_dir <- argv$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

if (!is.na(argv$`effective-sample-size`)) {
    sumstats[, n_eff := argv$`effective-sample-size`]
} else if (!is.na(argv$`n-cases`) && !is.na(argv$`n-controls`)) {
    n_eff_derived <- 4 / (1 / argv$`n-cases` + 1 / argv$`n-controls`)
    sumstats[, n_eff := n_eff_derived]
}

sumstats <- sumstats[!is.na(n_eff) & n_eff > 0]

if (nrow(sumstats) < 100) {
    fail_with(failed_match, sprintf("Too few valid sumstat rows after QC: %d", nrow(sumstats)))
}

# --- LD reference map ---------------------------------------------------------
map_ldref <- readRDS(ld_meta)
build <- toupper(argv$`genotype-build`)
if (build %in% c("GRCH38", "HG38")) {
    if (!("pos_hg38" %in% colnames(map_ldref))) {
        fail_with(failed_match, "LD map missing pos_hg38 for GRCh38 matching")
    }
    map_ldref$pos <- map_ldref$pos_hg38
} else if (!("pos" %in% colnames(map_ldref))) {
    fail_with(failed_match, "LD map missing pos column")
}

if (!("block_id" %in% colnames(map_ldref)) && ("group_id" %in% colnames(map_ldref))) {
    map_ldref$block_id <- map_ldref$group_id
}

join_by_pos <- !isTRUE(argv$`merge-by-rsid`)
cat(sprintf("[OK] Matching %d sumstat rows to LD reference (join_by_pos=%s)\n",
            nrow(sumstats), join_by_pos))

df_beta <- tryCatch(
    snp_match(sumstats, map_ldref, join_by_pos = join_by_pos, match.min.prop = 0),
    error = function(e) {
        fail_with(failed_match, conditionMessage(e))
        NULL
    }
)

drops <- c("_NUM_ID_.ss", "_NUM_ID_", "rsid.ss")
df_beta <- df_beta[, !(names(df_beta) %in% drops), drop = FALSE]

if (nrow(df_beta) < 100) {
    fail_with(
        failed_match,
        sprintf("snp_match retained too few variants: %d", nrow(df_beta))
    )
}
cat(sprintf("[OK] snp_match retained %d variants\n", nrow(df_beta)))

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

if (nrow(df_beta) < 100) {
    fail_with(failed_match, "Too few variants remain after allelic SD QC")
}

# --- Build genome-wide SFBM ---------------------------------------------------
tmp_file <- tempfile(tmpdir = out_dir, pattern = "ldpred2_corr_")
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
        ncores = argv$ncores
    )),
    error = function(e) {
        fail_with(failed_ldsc, conditionMessage(e))
        NULL
    }
)

h2_est <- ldsc[["h2"]]
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
} else {
    cat("[OK] Running LDpred2-auto\n")
    set.seed(argv$seed)
    multi_auto <- tryCatch(
        snp_ldpred2_auto(
            corr, df_beta,
            h2_init = h2_est,
            vec_p_init = seq_log(1e-4, argv$`hyper-p-max`, length.out = argv$`hyper-p-length`),
            allow_jump_sign = argv$`allow-jump-sign`,
            shrink_corr = argv$`shrink-corr`,
            ncores = argv$ncores
        ),
        error = function(e) {
            fail_with(failed_ldpred2, conditionMessage(e))
            NULL
        }
    )

    auto_res <- get_betas_auto(multi_auto)
    if (length(auto_res$keep) == 0) {
        fail_with(failed_ldpred2, "All LDpred2-auto chains diverged")
    }
    beta_est <- auto_res$beta
    postp_est <- auto_res$postp

    if (nzchar(argv$`plot-file`)) {
        auto <- multi_auto[[1]]
        dta <- data.frame(
            path_p_est = auto$path_p_est,
            path_h2_est = auto$path_h2_est,
            x = seq_along(auto$path_p_est)
        )
        plt <- plot_grid(
            ggplot(dta, aes(y = path_p_est, x = x)) +
                geom_point() +
                theme_bigstatsr() +
                geom_hline(aes(yintercept = auto$p_est), col = "blue") +
                scale_y_log10() +
                labs(y = "p"),
            ggplot(dta, aes(y = path_h2_est, x = x)) +
                geom_point() +
                theme_bigstatsr() +
                geom_hline(aes(yintercept = auto$h2_est), col = "blue") +
                labs(y = "h2"),
            ncol = 1,
            align = "hv"
        )
        ggsave(argv$`plot-file`, plt, width = 6, height = 8)
        cat(sprintf("[OK] Wrote chain diagnostics: %s\n", argv$`plot-file`))
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

cat(sprintf("[OK] Wrote posteriors to %s\n", out_dir))
invisible(TRUE)
