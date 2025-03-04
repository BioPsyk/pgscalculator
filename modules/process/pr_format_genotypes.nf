#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process make_geno_pvar_snpid_unique_pvar {
    publishDir "${params.outdir}/intermediates/make_geno_pvar_snpid_unique_pvar", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(pvar)

    output:
        tuple val(chr), path("${chr}_geno.pvar")
    
    script:
        """
        make_geno_pvar_snpid_unique.sh ${pvar} "${chr}_geno.pvar"
        """
}

process make_geno_pvar_snpid_unique_pvar_psam_pgen {
    publishDir "${params.outdir}/intermediates/make_geno_pvar_snpid_unique_pvar_psam_pgen", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(pgen), path(pvar), path(psam)

    output:
        tuple val(chr), path(pgen), path("${chr}_geno.pvar"), path(psam)
    
    script:
        """
        make_geno_pvar_snpid_unique.sh ${pvar} "${chr}_geno.pvar"
        """
}

process add_rsid_to_genotypes {
    //publishDir "${params.outdir}/intermediates/add_rsid_to_genotypes", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path("geno.pgen"), path("geno.pvar"), path("geno.psam"), path(rsid_ref)

    output:
        tuple val(chr), path("geno2.pgen"), path("geno2.pvar"), path("geno2.psam")
    
    script:
        """
        # add rsid and make list of tokeep
        add_rsid_to_genotypes.sh "geno.pvar" ${rsid_ref} > "geno_with_rsid.pvar"
        mv geno_with_rsid.pvar geno.pvar

        # only keep variants with rsid, remove dups and indels
        plink2 --pfile geno --make-pgen --memory ${params.memory.plink.add_rsid_to_genotypes} --threads 1 --extract tokeep --out geno2
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

