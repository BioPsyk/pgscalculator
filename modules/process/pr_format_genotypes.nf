#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process add_rsid_to_genotypes {
    publishDir "${params.outdir}/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'mod_mem'
    
    input:
        tuple val(chr), path("geno.bed"), path("geno.bim"), path("geno.fam"), path(rsid_ref)

    output:
        tuple val(chr), path("geno.bed"), path("${chr}_geno_rsid.bim"), path("geno.fam")
    
    script:
        """
        
        add_rsid_to_genotypes.sh "geno.bim" ${rsid_ref} > "${chr}_geno_rsid.bim"
        """
}

