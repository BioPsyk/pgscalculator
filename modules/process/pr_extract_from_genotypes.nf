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
        
        # Post-process the .afreq file to extract correct columns using robust method
        awk -vOFS=" " '
            NR==1 {
                # Map column names to positions
                for(i=1; i<=NF; i++) {
                    gsub(/^#/, "", \\$i);  # Remove # prefix if present
                    col[\\$i] = i;
                }
                next;
            }
            NR>1 {
                # Extract columns by name: CHR, SNP, REF, ALT, ALT_FREQS, OBS_CT
                chr = (col["CHROM"] ? \\$col["CHROM"] : \\$col["CHR"]);
                snp = \\$col["ID"];
                a1 = \\$col["REF"];
                a2 = \\$col["ALT"];
                maf = \\$col["ALT_FREQS"];
                obs = \\$col["OBS_CT"];
                print chr, snp, a1, a2, maf, obs;
            }
        ' ${chr}_geno_maf.afreq > ${chr}_geno_maf.afreq.tmp
        mv ${chr}_geno_maf.afreq.tmp ${chr}_geno_maf.afreq
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

