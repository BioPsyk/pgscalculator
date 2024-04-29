process prepare_sumstat_for_benchmark_scoring {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sumstat_map),path(sumstat_map_noNA),path(map),path(map_noNA)

    output:
    tuple val(chr), path("${chr}_sumstat")

    script:
    """
    colinx="\$(find_col_indices.sh ${sumstat_map_noNA} "RSID,EffectAllele,B")"
    echo "\${colinx}" > colinx
    format_posteriors.sh ${sumstat_map_noNA} "\${colinx}" ${map_noNA} "3,6" "false" > "${chr}_sumstat"
    """
}

process sumstat_maf_filter {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path("sumstat"),path("maffile")
    val(maf)

    output:
    tuple val(chr), path("${chr}_sumstat")

    script:
    """
    sumstat_maf_filter.sh ${sumstat} ${maffile} ${maf} > "${chr}_sumstat"
    """
}


// rename to indep_pairwise_filter_for_benchmark
process indep_pairwise_for_benchmark {
    publishDir "${params.outdir}/intermediates/indep_pairwise_for_benchmark", mode: 'rellink', overwrite: true
    
    input:
        tuple val(chr), path(sumstat), path("geno.pgen"), path("geno.pvar"), path("geno.psam")

    output:
        tuple val(chr), path("chr${chr}_sumstat")
    
    script:
        """
        awk 'NR>1{print \$4}' ${sumstat} > snps
        plink2 --pfile geno \
         --indep-pairwise 500kb 1 0.2 \
         --threads 1 \
         --memory 1000 \
         --extract snps \
         --rm-dup force-first \
         --out "chr${chr}"

        awk 'NR==FNR{a[\$1]; next} FNR==1 || (\$4 in a)' "chr${chr}.prune.in" ${sumstat} > "chr${chr}_sumstat" 
        """
}

