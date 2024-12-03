#!/usr/bin/env bash

set -euo pipefail

sumstat=$1
prefix=$2
chr=$3

# Extract SNPs from sumstat
awk 'NR>1{print $3}' ${sumstat} > snps

# Run PLINK2 pruning
if ! plink2 --pfile ${prefix} \
    --indep-pairwise 500kb 1 0.2 \
    --threads 1 \
    --memory 1000 \
    --extract snps \
    --rm-dup force-first \
    --out chr${chr}; then
    # If PLINK fails due to no variants in LD, use all variants
    echo "Warning: No variants in LD found, using all variants" >&2
    cp ${sumstat} chr${chr}_sumstat
    exit 0
fi

# Filter sumstat based on pruned SNPs (only if pruning succeeded)
if [ -f chr${chr}.prune.in ]; then
    awk 'NR==FNR{a[$1]; next} FNR==1 || ($3 in a)' chr${chr}.prune.in ${sumstat}
else
    # Fallback to using all variants if no pruned set is created
    cp ${sumstat} chr${chr}_sumstat
fi