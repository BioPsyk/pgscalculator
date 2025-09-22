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
        
        # Process the .afreq file using robust column extraction script
        process_plink_maf_single.sh ${chr}_geno_maf.afreq ${chr}_geno_maf.afreq.processed
        mv ${chr}_geno_maf.afreq.processed ${chr}_geno_maf.afreq
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
        # Create header
        echo "CHR SNP A1 A2 MAF NCHROBS" > "raw_maf_chrall"
        
        # Simply concatenate the already-processed files
        for chrfile in ${chrplinkmaf}
        do
          cat \$chrfile >> "raw_maf_chrall"
        done
        """
}

