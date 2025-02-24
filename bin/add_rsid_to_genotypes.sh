#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pvar> <rsidmap>"
    exit 1
fi

pvar="$1"
ref_file="$2"

# Step 1: Print all header lines from PVAR file
grep "^#" "$pvar"

# Step 2: Create a temporary file with chr:pos as key for non-header lines
# Format: chr:pos chr pos id ref alt qual filter info linenum
#echo "Creating sorted variant file..." >&2
grep -v "^#" "$pvar" | \
    awk 'BEGIN {OFS="\t"} 
         { 
           key = $1":"$2;
           print key, $0, NR 
         }' > "temp_variants"

# Step 3: Sort the variants file for joining
LC_ALL=C sort -k1,1 "temp_variants" > "sorted_variants"

# Step 4: Process reference file to match format
# Input format: b37 b38 ss_SNP ss_A1 ss_A2 bim_SNP bim_A1 bim_A2 ld_SNP ld_A1 ld_A2
#echo "Processing reference file..." >&2
awk 'BEGIN {OFS="\t"} 
     NR > 1 { 
       print $1, $3
     }' "$ref_file" | LC_ALL=C sort -k1,1 > "sorted_ref"

# Step 5: Join variants with reference file and format output
#echo "Joining files and creating output..." >&2
LC_ALL=C join -a 1 -e "." -1 1 -2 1 \
    "sorted_variants" "sorted_ref" | \
    sort -n -k10,10 | \
    awk -vOFS="\t" '
    {
        # Store variant info for duplicate/biallelic checking
        key = $11;  # rsid
        if (length($5)==1 && length($6)==1) {
            biallelic[key]++;
        }
        seen[key]++;
        
        # Print original line but replace ID with rsid
        print $2, $3, $11, $5, $6, $7, $8, $9
    }
    END {
        # Output rsids that are unique and biallelic to tokeep file
        for (k in seen) {
            if (seen[k]==1 && biallelic[k]==1 && k != ".") {
                print k > "tokeep"
            }
        }
    }'

# Step 6: Clean up temporary files
#echo "Cleaning up..." >&2
rm "temp_variants" "sorted_variants" "sorted_ref"


