#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process convert_plink1_to_plink2 {
    publishDir "${params.outdir}/intermediates/convert_plink1_to_plink2", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(bed), path(bim), path(fam)

    output:
        tuple val(chr), path("${chr}_converted.pgen"), path("${chr}_converted.pvar"), path("${chr}_converted.psam")
    
    script:
        """
        # Convert PLINK1 format to PLINK2 format
        plink2 --bed ${bed} --bim ${bim} --fam ${fam} \\
               --make-pgen \\
               --memory ${params.memory.plink.extract_maf_from_genotypes} \\
               --threads 1 \\
               --out ${chr}_converted
        """
}

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

process standardize_psam_to_iid_only {
    publishDir "${params.outdir}/intermediates/standardize_psam", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(pgen), path(pvar), path(psam)

    output:
        tuple val(chr), path(pgen), path(pvar), path("${chr}_standardized.psam")
    
    script:
        """
        # Standardize .psam to IID-only format
        if head -1 ${psam} | grep -q "^#FID"; then
            # Has FID column - remove it, keep IID and other columns
            awk 'NR==1{gsub(/^#FID[[:space:]]+/, "#"); print} NR>1{\$1=""; gsub(/^[[:space:]]+/, ""); print}' ${psam} > "${chr}_standardized.psam"
        else
            # Already IID-only format - just copy
            cp ${psam} "${chr}_standardized.psam"
        fi
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

