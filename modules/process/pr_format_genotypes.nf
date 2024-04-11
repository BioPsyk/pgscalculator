#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process add_rsid_to_genotypes {
    //publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'mod_mem'
    
    input:
        tuple val(chr), path("geno.bed"), path("genoX.bim"), path("geno.fam"), path(rsid_ref)

    output:
        tuple val(chr), path("geno2.pgen"), path("geno2.pvar"), path("geno2.psam")
    
    script:
        """
        
        # add rsid and make list of tokeep
        add_rsid_to_genotypes.sh "genoX.bim" ${rsid_ref} > "geno.bim"

        # only keep variants with rsid, remove dups and indels
        plink2 --bfile geno --make-pgen --extract tokeep --out geno2
        """
}

process concat_genotypes {
    publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'mod_mem'
    
    input:
       path("allfiles.txt")

    output:
        tuple val("all"), path("allgeno.pgen"), path("allgeno.pvar"), path("allgeno.psam")
    
    script:
        """
        plink2 --pmerge-list allfiles.txt --make-pgen --multiallelics-already-joined --out allgeno
        """

}

