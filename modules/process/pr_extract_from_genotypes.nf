// Extract data from genotypes, such as maf, etc.
process extract_maf_from_genotypes {  
    publishDir "${params.outdir}/intermediates/maf_from_genotypes", mode: 'rellink', overwrite: true  

    input:  
        tuple val(chr), path("geno.bed"), path("geno.bim"), path("geno.fam") 

    output:                                                                                 
        tuple val(chr), path("${chr}_geno_maf.frq")    

    script:                                                                                 
        """
        plink --bfile geno --freq --out ${chr}_geno_maf
        """                                                                                 
}

process concatenate_plink_maf {
    publishDir "${params.outdir}/extra", mode: 'copy', overwrite: true

    cpus 22

    input:
        path(chrplinkmaf)
    
    output:
        tuple val("all"), path("raw_maf_chrall")
    script:
        """
        echo " CHR           SNP   A1   A2          MAF  NCHROBS"  > "raw_maf_chrall"
        parallel tail -n+2 ::: ${chrplinkmaf} >> "raw_maf_chrall"

       # for chrfile in ${chrplinkmaf}
       # do
       #   tail -n+2 \$chrfile >> "raw_maf_chrall"
       # done
        """
}


