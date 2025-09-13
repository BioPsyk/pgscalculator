#!/bin/bash

# Script to process a single PLINK MAF file with robust column extraction
# Uses column names instead of positional indices for better maintainability

set -euo pipefail

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_afreq_file> <output_file>"
    echo "Example: $0 chr1_geno_maf.afreq chr1_processed_maf.txt"
    exit 1
fi

input_file="$1"
output_file="$2"

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file $input_file does not exist"
    exit 1
fi

# Process the file using robust column extraction
awk -vOFS=" " '
    NR==1 {
        # Map column names to positions
        for(i=1; i<=NF; i++) {
            gsub(/^#/, "", $i);  # Remove # prefix if present
            col[$i] = i;
        }
        next;
    }
    NR>1 {
        # Extract columns by name: CHR, SNP, REF, ALT, ALT_FREQS, OBS_CT
        chr = (col["CHROM"] ? $col["CHROM"] : $col["CHR"]);
        snp = $col["ID"];
        a1 = $col["REF"];
        a2 = $col["ALT"];
        maf = $col["ALT_FREQS"];
        obs = $col["OBS_CT"];
        print chr, snp, a1, a2, maf, obs;
    }
' "$input_file" > "$output_file"

echo "Successfully processed $input_file -> $output_file"
