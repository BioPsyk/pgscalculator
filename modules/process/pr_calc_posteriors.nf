// calculate per chromosome posterior SNP effects for sBayesR

process calc_posteriors_prscs {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true
    label 'big_mem'

    input:
        tuple val(chr), path(gwas), path(lddir), path("geno.bed"), path("geno.bim"), path("geno.fam"), val(N)
            
    output:
        tuple val(chr), path("chr${chr}.posteriors")

    script:
        def integerN = Math.round(N as float)
        """
        mkdir prscs

        # Count the number of lines in the file
        num_lines=\$(head -n20 "$gwas" | wc -l)

        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          # PRScs requires python2, which is best accessed through conda
          #source /opt/micromamba/bin/activate py27

          micromamba run -n py27 python /repos/PRScs/PRScs.py \
            --ref_dir=$lddir \
            --sst_file=${gwas} \
            --bim_prefix="geno" \
            --n_gwas=${integerN} \
            --chrom=${chr} \
            --out_dir=prscs

          mv prscs_pst_eff_a1_b0.5_phiauto_chr${chr}.txt chr${chr}.posteriors
        else
          touch chr${chr}.posteriors
        fi
        """ 
}

// calculate per chromosome posterior SNP effects for sBayesR
process calc_posteriors_sbayesr {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true

    label 'big_mem'
    cpus 6
    
    input:
        tuple val(chr), path(gwas_chr), path(ldbin), path(ldinfo)
    
    output:
        tuple val(chr), path("chr${chr}.posteriors")

    script:
        ld_prefix="band_chr${chr}.ldm.sparse"

        def options = params.calc_posteriors_sbayesr.options.collect { key, value ->
            "--${key.replace('_', '-')} ${value}"
        }.join(' ')
        
        def flags = params.calc_posteriors_sbayesr.flags.findAll { it.value }
                      .collect { key, _ -> "--${key.replace('_', '-')}" }
                      .join(' ')

        """
        # Count the number of lines in the file
        num_lines=\$(head -n20 "$gwas_chr" | wc -l)
        
        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          gctb --sbayes R \
            --gwas-summary ${gwas_chr} \
            --ldm ${ld_prefix} \
            $options \
            $flags \
            --out chr${chr}
          mv "chr${chr}.snpRes" "chr${chr}.posteriors"
        else
          touch "chr${chr}.posteriors"
        fi
        """
}


// Concatenate per chromosome posterior SNP effects for sBayesR
process concatenate_sbayes_posteriors {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true

    input:
        path(chrposteriors)
    
    output:
        tuple val("all"), path("allchr.posteriors")
    script:
        """
        echo "    Id                 Name  Chrom     Position     A1     A2        A1Frq     A1Effect           SE            PIP  LastSampleEff"  > "allchr.posteriors"
        for chrfile in ${chrposteriors}
        do
          tail -n+2 \$chrfile >> "allchr.posteriors"
        done
        """
}

// calculate per chromosome posterior SNP effects for sBayesR
process qc_posteriors {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true

    label 'big_mem'
    cpus 6
    
    input:
        tuple val(chr), path("chr${chr}.posteriors")
    
    output:
        tuple val(chr), path("chr${chr}.png")

    script:

        """
        qc_posteriors.sh "/pgscalculator/bin/R/qc_posteriors.R" "chr${chr}.posteriors"
        """
}

