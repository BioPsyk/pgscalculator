#!/usr/bin/env bash
# Smoke-test that the R packages added in the multi-stage Dockerfile (`r_builder`
# stage) load cleanly inside the built container. Run this AFTER docker-build.sh
# (or after singularity-build.sh — see usage below).
#
# Usage (docker):
#   ./scripts/test-r-packages.sh
#   ./scripts/test-r-packages.sh ibp-pgscalculator-base:0.7.0    # explicit image
#
# Usage (singularity / apptainer):
#   ./scripts/test-r-packages.sh --singularity tmp/ibp-pgscalculator-base_version-0.7.0.sif

set -euo pipefail

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${script_dir}/init-containerization.sh"

# Resolve target image / runtime
mode="docker"
target=""
if [[ "${1:-}" == "--singularity" ]]; then
    mode="singularity"
    target="${2:-}"
    [[ -z "$target" ]] && { echo "Error: --singularity requires a .sif path"; exit 2; }
    [[ ! -f "$target" ]] && { echo "Error: .sif not found: $target"; exit 2; }
elif [[ -n "${1:-}" ]]; then
    target="$1"
else
    target="${image_tag}"
fi

# R packages we expect to be present (from r_builder).
# Add new packages here when extending the Dockerfile.
read -r -d '' R_SCRIPT << 'EOF' || true
pkgs <- c("bigsnpr", "bigsparser", "bigstatsr", "bigreadr", "runonce",
          "argparser", "stringr", "ggplot2", "cowplot", "data.table",
          "tibble", "tidyr", "dplyr")
status <- 0
for (p in pkgs) {
  ok <- suppressPackageStartupMessages(
    requireNamespace(p, quietly = TRUE)
  )
  cat(sprintf("%-15s %s  (%s)\n",
              p,
              if (ok) "OK" else "MISSING",
              if (ok) as.character(packageVersion(p)) else "-"))
  if (!ok) status <- 1
}

# Hard pins from the Dockerfile.
pinned <- list(
  bigsnpr    = "1.12.0",
  bigsparser = "0.6.0"
)
for (p in names(pinned)) {
  if (requireNamespace(p, quietly = TRUE)) {
    v <- packageVersion(p)
    if (v < pinned[[p]]) {
      cat(sprintf("FAIL: %s version %s is older than required %s\n",
                  p, v, pinned[[p]]))
      status <- 1
    }
  }
}

# Tiny functional check: build an SFBM from a 3x3 sparse correlation matrix
# and call snp_ldpred2_inf on a toy df_beta. Catches ABI mismatches and
# missing shared libs at the function level (not just `library()` load).
if (status == 0) {
  cat("\n--- functional check ---\n")
  ok <- tryCatch({
    suppressPackageStartupMessages(library(bigsnpr))
    suppressPackageStartupMessages(library(Matrix))
    M <- Matrix::Matrix(diag(3), sparse = TRUE)
    tmp <- tempfile()
    corr <- bigsparser::as_SFBM(M, tmp, compact = TRUE)
    df_beta <- data.frame(
      beta    = c(0.1, -0.05, 0.02),
      beta_se = c(0.01,  0.01, 0.01),
      n_eff   = c(1e4,   1e4,  1e4)
    )
    beta <- snp_ldpred2_inf(corr, df_beta, h2 = 0.1)
    file.remove(paste0(tmp, ".sbk"))
    cat("snp_ldpred2_inf returned", length(beta), "effects: OK\n")
    TRUE
  }, error = function(e) {
    cat("FAIL: functional check errored:", conditionMessage(e), "\n")
    FALSE
  })
  if (!ok) status <- 1
}

quit(status = status, save = "no")
EOF

echo ">> Smoke-testing R packages in ${target} (mode=${mode})"
echo

case "$mode" in
    docker)
        docker run --rm -i "$target" R --slave --no-init-file <<< "$R_SCRIPT"
        ;;
    singularity)
        singularity exec "$target" R --slave --no-init-file <<< "$R_SCRIPT"
        ;;
esac

echo
echo ">> R packages OK"
