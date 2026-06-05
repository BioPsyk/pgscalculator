# =============================================================================
# b3sum builder (not packaged for Debian)
# =============================================================================
FROM rust:1.84-slim-bookworm AS rust_builder

RUN cargo install --version 0.3.8 b3sum

# =============================================================================
# R packages builder (isolated; rebuilds only when this stage changes)
# =============================================================================
# This stage exists so that R-package installation — which is slow because of
# bigsnpr's compiled dependencies — is cached independently from the rest of
# the final image. Changes to apt-get lines, PRScs/gctb/plink versions, or
# bin/ in the final stage no longer bust the R package cache.
#
# The final stage uses the SAME base image (`eclipse-temurin:11`) and apt-installs
# the SAME `r-base` package, so the R version in r_builder matches the runtime
# R version (avoids ABI mismatch).
FROM eclipse-temurin:11 AS r_builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libbz2-dev \
        liblzma-dev \
        libgomp1 \
        build-essential \
        pkg-config \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install pak: parallel installer that pulls prebuilt binaries from Posit
# Public Package Manager where available. Cuts bigsnpr cold install from
# ~25 min to ~3-5 min.
RUN R -e "install.packages('pak', repos = sprintf( \
    'https://r-lib.github.io/p/pak/stable/%s/%s/%s', \
    .Platform\$pkgType, R.Version()\$os, R.Version()\$arch))"

# All R packages needed for LDpred2.
# Pin: bigsnpr >= 1.12.0 (validated LDpred2-auto, paper-infer)
#      bigsparser >= 0.6 (compact SFBM format)
RUN R -e "pak::pkg_install(c( \
        'bigsnpr', \
        'bigsparser', \
        'bigreadr', \
        'runonce', \
        'argparser', \
        'stringr', \
        'ggplot2', \
        'cowplot', \
        'data.table', \
        'tibble', \
        'tidyr', \
        'dplyr')); \
    stopifnot(packageVersion('bigsnpr') >= '1.12.0'); \
    stopifnot(packageVersion('bigsparser') >= '0.6')"

# =============================================================================
# Final runtime image
# =============================================================================
FROM eclipse-temurin:11

COPY --from=rust_builder /usr/local/cargo/bin/b3sum /usr/bin/b3sum

WORKDIR /pgscalculator

# System libs. Kept identical to the pre-r_builder version of this Dockerfile
# (i.e. -dev variants, not runtime-only) to avoid base-image-codename mismatch:
# `eclipse-temurin:11` is unpinned in this repo and may resolve to Focal or
# Jammy at build time, which would change libssl runtime names (libssl1.1
# vs libssl3). The -dev libs implicitly drag in the correct runtime variants.
# Slimming to runtime-only is tracked as a follow-up.
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    curl \
    gawk \
    graphviz \
    dos2unix \
    pigz \
    libbz2-dev \
    liblzma-dev \
    vim \
    parallel \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgomp1 \
    git \
    --no-install-recommends && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Pull the prebuilt R library tree from the builder. On Debian/Ubuntu, pak
# installs into /usr/local/lib/R/site-library when running as root. We also
# copy /usr/lib/R/site-library in case any package landed there (harmless
# no-op if empty).
COPY --from=r_builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=r_builder /usr/lib/R/site-library       /usr/lib/R/site-library

# Install PRS-CS
WORKDIR /repos
RUN git clone https://github.com/getian107/PRScs.git PRScs

# Install GCTB (sbayesR)
WORKDIR /pgscalculator
ARG GCTB_URL_PRIMARY=https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.05beta_Linux.zip
ARG GCTB_URL_FALLBACK=https://cnsgenomics.com/software/gctb/download/gctb_2.05beta_Linux.zip
RUN (wget "${GCTB_URL_PRIMARY}" --no-check-certificate || wget "${GCTB_URL_FALLBACK}" --no-check-certificate) && \
    unzip gctb_2.05beta_Linux.zip && \
    mv gctb_2.05beta_Linux/gctb /usr/bin/gctb

# Install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip && \
    unzip plink_linux_x86_64_20231211.zip && \
    mv plink /usr/local/bin/

# Install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip && \
    unzip plink2_linux_x86_64_20240105.zip && \
    mv plink2 /usr/local/bin/

# Set up directories with correct permissions
RUN useradd -rm -s /bin/bash -g users pgsuser && \
    chown -R pgsuser:users /pgscalculator && \
    chmod -R 777 /pgscalculator && \
    chmod g+s /pgscalculator

USER pgsuser

# pgscalculator v2 CLI
COPY bin/ /pgscalculator/bin/
COPY assets/ /pgscalculator/assets/
COPY config.template.yaml /pgscalculator/
RUN chmod +x /pgscalculator/bin/pgscalculator
ENV PATH="/pgscalculator/bin:${PATH}"
