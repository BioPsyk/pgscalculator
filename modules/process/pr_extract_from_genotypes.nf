// Extract data from genotypes, such as maf, etc.
process extract_maf_from_genotypes {  
    publishDir "${params.outdir}/intermediates/maf_from_genotypes", mode: 'rellink', overwrite: true  
    label 'low_mem'

    input:  
        tuple val(chr), path("geno.bed"), path("geno.bim"), path("geno.fam"), path("map"), path("map_noNA")

    output:                                                                                 
        tuple val(chr), path("${chr}_geno_maf.afreq")    

    script:                                                                                 
        """
        cut -f 6 $map > bimIDs
        plink2 --bfile geno --extract bimIDs --threads 1 --memory 1000 --freq --out ${chr}_geno_maf
        """                                                                                 
}

process concatenate_plink_maf {
    publishDir "${params.outdir}/extra", mode: 'copy', overwrite: true

    cpus 6
    label 'low_mem'

    input:
        path(chrplinkmaf)
    
    output:
        tuple val("all"), path("raw_maf_chrall")
    script:
        """
        echo "CHR SNP A1 A2 MAF NCHROBS"  > "raw_maf_chrall"
        for chrfile in ${chrplinkmaf}
        do
          awk -vOFS=" " 'NR>1{print \$1, \$2, \$3, \$4, \$6, \$7}' \$chrfile >> "raw_maf_chrall"
        done
        """
}

