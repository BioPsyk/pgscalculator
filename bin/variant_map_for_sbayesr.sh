#!/bin/bash

# Usage: ./script.sh <build_version> <ss2_file> <bim2_file> <ld2_file> <output_file>
# <build_version> is either b37 or b38

ss2_file="${1}"
bim2_file="${2}"
snp2_sorted="${3}"
ld2_file="${4}"
gbuild=${5}
lbuild=${6}
output_file="${7}"

head ${snp2_sorted} > snp2_head

if [ ${gbuild} == "38" ]; then
  # if build 38, switch b37 and b38 columns, and then switch back as last step
  awk -vFS="\t" -vOFS="\t" '{print $2,$1,$3,$4,$5}' ${ss2_file} | LC_ALL=C  sort -k1,1 > ss2_sorted.tmp
else
  # if no switch, just sort
  sort -k1,1 ${ss2_file} > ss2_sorted.tmp
fi

#echo "----" >&2
#cat ss2_sorted.tmp >&2
#echo "----" >&2

# Create temporary files with sorted content for join
LC_ALL=C sort -k4,4 ${bim2_file} > bim2_sorted_1.tmp
LC_ALL=C sort -k1,1 ${ld2_file} > ld2_sorted.tmp

# Perform the joins
LC_ALL=C join -1 1 -2 4 -o 2.1 2.2 2.3 2.4 "${snp2_sorted}" bim2_sorted_1.tmp > join0.tmp
LC_ALL=C sort -k1,1 join0.tmp > bim2_sorted_2.tmp

LC_ALL=C join -1 1 -2 1 -o 1.1 1.2 1.5 1.3 1.4 2.4 2.2 2.3 ss2_sorted.tmp bim2_sorted_2.tmp > join1.tmp


if [ "${lbuild}" == "${gbuild}" ] ; then
#  echo "hej ${gbuild} and ${lbuild}" >&2
  LC_ALL=C  join -1 1 -2 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.4 2.2 2.3 join1.tmp ld2_sorted.tmp > join2.tmp
else
  LC_ALL=C sort -k2,2 join1.tmp > join1.resorted.tmp
  LC_ALL=C join -1 2 -2 1 -a 1 -e 'NA' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.4 2.2 2.3 join1.resorted.tmp ld2_sorted.tmp > join2.tmp
fi

# Filter rows where neither A1==A1 nor A1==A2 between the files
# Assuming columns for A1/A2 are as follows in the final join: ss_A1=4, ss_A2=5, bim_A1=7, bim_A2=8, ld_A1=10, ld_A2=11
#awk -vFS="\t" '{
#    if ( (($4==$7 || $4==$8) && ($4==$10 || $4==$11)) && (($5==$7 || $5==$8) && ($5==$10 || $5==$11)) ) print
#}' join2.tmp > filtered_join2.tmp

awk -vFS=" " -vOFS="\t" '
BEGIN{
  f["A"]="T"
  f["T"]="A"
  f["G"]="C"
  f["C"]="G"
}
{
    # Conditions for matching ss alleles with ld alleles, or ld alleles are NA
    # match_ld_na = $10 == "NA" && $11 == "NA"
    match_ss_ld = ($4 == $10 || $4 == $11) && ($5 == $10 || $5 == $11);
    match_ss_ld_flip = (f[$4] == $10 || f[$4] == $11) && (f[$5]== $10 || f[$5] == $11);
    if (!match_ss_ld && !match_ss_ld_flip) {
      $10 == "NA" 
      $11 == "NA"
    }

    # Conditions for matching ss alleles with bim alleles
    match_ss_bim = ($4 == $7 || $4 == $8) && ($5 == $7 || $5 == $8);
    match_ss_bim_flip = (f[$4] == $7 || f[$4] == $8) && (f[$5] == $7 || f[$5] == $8);
    if (match_ss_bim || match_ss_bim_flip) {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11;
    }
}' join2.tmp > filtered_join2.tmp


# Add a header
echo -e "b37\tb38\tss_SNP\tss_A1\tss_A2\tpvar_SNP\tpvar_A1\tpvar_A2\tld_SNP\tld_A1\tld_A2" > $output_file

if [ ${gbuild} == "38" ]; then
  # if build 38, switch b37 and b38 columns, and then switch back as last step
  awk 'BEGIN {FS="\t"; OFS="\t"} {temp = $1; $1 = $2; $2 = temp; print}' filtered_join2.tmp >> $output_file
else
  # Append the joined data to the output file
  cat filtered_join2.tmp >> $output_file
fi


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

