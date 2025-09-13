// Extract data from genotypes, such as maf, etc.
process extract_maf_from_genotypes {  
    publishDir "${params.outdir}/intermediates/maf_from_genotypes", mode: 'rellink', overwrite: true  

    input:  
        tuple val(chr), path("geno.bed"), path("geno.bim"), path("geno.fam"), path("map"), path("map_noNA")

    output:                                                                                 
        tuple val(chr), path("${chr}_geno_maf.afreq")    

    script:                                                                                 
        """
        cut -f 6 $map > bimIDs
        plink2 --bfile geno --extract bimIDs --threads 1 --memory ${params.memory.plink.extract_maf_from_genotypes} --freq --out ${chr}_geno_maf
        """                                                                                 
}

process concatenate_plink_maf {
    publishDir "${params.outdir}/extra", mode: 'copy', overwrite: true

    cpus 6

    input:
        path(chrplinkmaf)
    
    output:
        tuple val("all"), path("raw_maf_chrall")
    script:
        """
        concatenate_plink_maf.sh "raw_maf_chrall" ${chrplinkmaf}
        """
}

