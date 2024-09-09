#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process make_geno_bim_snpid_unique_bim {
    publishDir "${params.outdir}/intermediates/make_geno_bim_snpid_unique_bim", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(bim)

    output:
        tuple val(chr), path("${chr}_geno.bim")
    
    script:
        """
        make_geno_bim_snpid_unique.sh ${bim} "${chr}_geno.bim"
        """
}

process make_geno_bim_snpid_unique_bim_fam_bed {
    publishDir "${params.outdir}/intermediates/make_geno_bim_snpid_unique_bim_fam_bed", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(bed), path(bim), path(fam)

    output:
        tuple val(chr), path(bed), path("${chr}_geno.bim"), path(fam)
    
    script:
        """
        make_geno_bim_snpid_unique.sh ${bim} "${chr}_geno.bim"
        """
}

process add_rsid_to_genotypes {
    //publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path("geno.bed"), path("genoX.bim"), path("geno.fam"), path(rsid_ref)

    output:
        tuple val(chr), path("geno2.pgen"), path("geno2.pvar"), path("geno2.psam")
    
    script:
        """
        
        # add rsid and make list of tokeep
        add_rsid_to_genotypes.sh "genoX.bim" ${rsid_ref} > "geno.bim"

        # only keep variants with rsid, remove dups and indels
        plink2 --bfile geno --make-pgen --memory ${params.memory.plink.add_rsid_to_genotypes} --threads 1 --extract tokeep --out geno2
        """
}

process concat_genotypes {
    publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
       path("allfiles.txt")

    output:
        tuple val("all"), path("allgeno.pgen"), path("allgeno.pvar"), path("allgeno.psam")
    
    script:
        """
        plink2 --pmerge-list allfiles.txt --make-pgen --memory ${params.memory.plink.concat_genotypes} --threads 1 --multiallelics-already-joined --out allgeno
        """

}

