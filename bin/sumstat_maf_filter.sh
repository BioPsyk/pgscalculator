#!/bin/bash

# Usage: ./filter_variants.sh sumstat_file maf_file maf_threshold

sumstat_file="$1"
maf_file="$2"
maf_threshold="$3"

# Ensure the threshold is provided and is a number
if ! [[ "$maf_threshold" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: MAF threshold must be a number." > "/dev/stderr"
    exit 1
fi

# Process MAF file and sumstat file in a single awk command
head -n 1 "$sumstat_file"
awk -v threshold="$maf_threshold" '
    NR==FNR {if ($5 > threshold) variants[$2]; next}
    ($4 in variants)
' "$maf_file" "$sumstat_file"

