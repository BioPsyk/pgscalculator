#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bim> <rsidmap>"
    exit 1
fi

bim="$1"
ref_file="$2"


#bim="/home/jesgaaopen/ibp_migration_opengdk/PROJECT_pgscalculator/ibp_pgs_pipelines_2/references/genotypes_test/plink_old/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bim"
#build="37"
#b37="/home/jesgaaopen/ibp_migration_opengdk/PROJECT_pgscalculator/ibp_pgs_pipelines_2/references/cleansumstat_rsid_map/All_20180418_GRCh37_GRCh38.sorted.bed.gz"
#b38="/home/jesgaaopen/ibp_migration_opengdk/PROJECT_pgscalculator/ibp_pgs_pipelines_2/references/cleansumstat_rsid_map/All_20180418_GRCh38_GRCh37.sorted.bed.gz"

#[jesgaaopen@fe-open-01 plink_old]$ head ALL.chr18.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bim
#18      .       0       10083   <CN0>   A
#18      .       0       10644   G       C
#18      .       0       10648   A       C
#18      .       0       10652   A       G

#[jesgaaopen@fe-open-01 cleansumstat_rsid_map]$ zcat All_20180418_GRCh37_GRCh38.sorted.bed.gz | head
#10:10000000 10:9958037 rs1223870629 C T
#10:100000000 10:98240243 rs925217917 T C
#10:100000003 10:98240246 rs537453558 C T
#10:100000004 10:98240247 rs1382475200 G A,C
#10:100000005 10:98240248 rs530933119 G C

# Sort the BIM file using chr:pos for matching
awk '{ print $1":"$4, $1, $2, $3, $4, $5, $6, NR }' "$bim" | LC_ALL=C sort -k1,1 > "sorted_bim"

# Join the sorted BIM file with the reference file on the chr:pos column to fill in the rsids
LC_ALL=C join -a 1 -e "." -1 1 -2 1 -o 1.2,2.3,1.4,1.5,1.6,1.7,1.8 "sorted_bim" <(zcat "$ref_file") | sort -n -k7,7 > "matched_output" || {
    echo "An error occurred during join operation."
    exit 1
}

# skip nonmatching alleles and Prepare the final BIM file with rsids updated
awk -vOFS="\t" '{seen[$2]++; if(length($5)==1 && length($6)==1){biallelic[$2]++}; print $1,$2,$3,$4,$5,$6}; END{for (k in seen){if(seen[k]==1 && biallelic[k]==1){print k > "tokeep"}}}' "matched_output"

# Clean up the temporary files
rm "sorted_bim" "matched_output"


