#!/bin/bash

# Script to concatenate PLINK MAF files with robust column extraction
# Uses column names instead of positional indices for better maintainability

set -euo pipefail

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <output_file> <input_file1> [input_file2] ..."
    echo "Example: $0 raw_maf_chrall chr1.afreq chr2.afreq chr3.afreq"
    exit 1
fi

output_file="$1"
shift
input_files=("$@")

# Create header
echo "CHR SNP A1 A2 MAF NCHROBS" > "$output_file"

# Process each input file
for chrfile in "${input_files[@]}"; do
    if [[ ! -f "$chrfile" ]]; then
        echo "Warning: File $chrfile does not exist, skipping..."
        continue
    fi
    
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
    ' "$chrfile" >> "$output_file"
done

echo "Successfully concatenated ${#input_files[@]} files into $output_file"
