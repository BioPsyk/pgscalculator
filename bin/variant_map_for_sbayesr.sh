#!/bin/bash

# Usage: ./script.sh <build_version> <ss2_file> <bim2_file> <ld2_file> <output_file>
# <build_version> is either b37 or b38

ss2_file=$1
bim2_file=$2
ld2_file=$3
output_file=$4


# Create temporary files with sorted content for join
sort -k1,1 $bim2_file > bim2_sorted.tmp
sort -k1,1 $ld2_file > ld2_sorted.tmp
sort -k1,1 $ss2_file > ss2_sorted.tmp

# Perform the joins
join -1 1 -2 1 -t $'\t' -o 1.1 1.2 1.5 1.3 1.4 2.4 2.2 2.3 ss2_sorted.tmp bim2_sorted.tmp > join1.tmp
join -1 1 -2 1 -t $'\t' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.4 2.2 2.3 join1.tmp ld2_sorted.tmp > join2.tmp

# Filter rows where neither A1==A1 nor A1==A2 between the files
# Assuming columns for A1/A2 are as follows in the final join: ss_A1=4, ss_A2=5, bim_A1=7, bim_A2=8, ld_A1=10, ld_A2=11
awk -vFS="\t" '{
    if ( (($4==$7 || $4==$8) && ($4==$10 || $4==$11)) && (($5==$7 || $5==$8) && ($5==$10 || $5==$11)) ) print
}' join2.tmp > filtered_join.tmp

# Add a header
echo -e "b37\tb38\tss_SNP\tss_A1\tss_A2\tbim_SNP\tbim_A1\tbim_A2\tld_SNP\tld_A1\tld_A2" > $output_file

# Append the joined data to the output file
cat filtered_join.tmp >> $output_file

## Step 1: Filter tracking after the first join
## Identify rows from ss2 not in join1
#comm -23 ss2_sorted.tmp join1.tmp > filtered_out_after_join1_ss2.tmp
## Identify rows from bim2 not in join1
#comm -23 bim2_sorted.tmp join1.tmp > filtered_out_after_join1_bim2.tmp
#
## Step 2: Filter tracking after the second join
## Identify rows from join1 not in join2 (effectively tracking from both ss2 and bim2 as join1 is their combination)
#comm -23 join1.tmp join2.tmp > filtered_out_after_join2_ld.tmp
#

