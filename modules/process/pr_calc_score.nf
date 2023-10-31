#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process calc_score_prscs {
    publishDir "${params.outdir}/scores", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'mod_mem'
    
    input:
        tuple val(chr), path(snp_posteriors), path("geno.bed"), path("geno.bim"), path("geno.fam")

    output:
        path("chr${chr}.score")
    
    script:
        """
        ./plink2 --bfile geno \
        --out tmp \
        --score ${snp_posteriors} 2 4 6 header cols=+scoresums ignore-dup-ids

        awk '{gsub(/^#/, ""); print}' tmp.score > chr${chr}.score
        """
}

