// calculate per chromosome posterior SNP effects for prscs
process calc_posteriors_prscs {
    publishDir "${params.outdir}/intermediates/calc_posteriors", mode: 'rellink', overwrite: true
    label 'big_mem'

    input:
        tuple val(chr), path(gwas), path(lddir), path("geno.bim"), val(N)
            
    output:
        tuple val(chr), path("chr${chr}.posteriors")

    script:
        def integerN = Math.round(N as float)
        def options = params.calc_posteriors_prscs.options.collect { key, value ->
            "--${key}=${value}"
        }.join(' ')
        """
        mkdir prscs

        # Count the number of lines in the file
        num_lines=\$(head -n20 "$gwas" | wc -l)

        # Check if the file has more than one line (more than just the header)
        if [ "\$num_lines" -gt 1 ]; then
          # PRScs requires python2, which is best accessed through conda
          micromamba run -n py27 python /repos/PRScs/PRScs.py \
            --ref_dir=$lddir \
            --sst_file=${gwas} \
            --bim_prefix="geno" \
            --n_gwas=${integerN} \
            --chrom=${chr} \
            $options \
            --out_dir=prscs

          mv prscs_pst_eff_a1_b0.5_phiauto_chr${chr}.txt chr${chr}.posteriors
        else
          touch chr${chr}.posteriors
        fi
        """ 
}

// calculate per chromosome posterior SNP effects for sBayesR
process calc_posteriors_sbayesr {
    publishDir "${params.outdir}/intermediates/calc_posteriors", mode: 'rellink', overwrite: true

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

// Concatenate per chromosome posterior SNP effects for PRSCS
process concatenate_sbayes_PRSCS {
    publishDir "${params.outdir}/calc_posteriors", mode: 'copy', overwrite: true

    input:
        path(chrposteriors)
    
    output:
        tuple val("all"), path("raw_posteriors_chrall")
    script:
        """
        echo "Chrom	Id	Position	A1	A2"  > "raw_posteriors_chrall"
        for chrfile in ${chrposteriors}
        do
          tail -n+2 \$chrfile >> "raw_posteriors_chrall"
        done
        """
}

// Concatenate per chromosome posterior SNP effects for sBayesR
process concatenate_sbayes_posteriors {
    publishDir "${params.outdir}/extra", mode: 'copy', overwrite: true

    input:
        path(chrposteriors)
    
    output:
        tuple val("all"), path("raw_posteriors_chrall")
    script:
        """
        echo "    Id                 Name  Chrom     Position     A1     A2        A1Frq     A1Effect           SE            PIP  LastSampleEff"  > "raw_posteriors_chrall"
        for chrfile in ${chrposteriors}
        do
          tail -n+2 \$chrfile >> "raw_posteriors_chrall"
        done
        """
}

// calculate per chromosome posterior SNP effects for sBayesR
process qc_posteriors {
    publishDir "${params.outdir}/qc_posteriors", mode: 'copy', overwrite: true

    label 'big_mem'
    cpus 6
    
    input:
        tuple val(chr), path("posteriors"), path("input")
    
    output:
        tuple val(chr), path('*.png')

    script:

        """
        if [ "\$(head -n2 ${posteriors} | wc -l)" -gt 1 ]; then
          qc_posteriors.sh "/pgscalculator/bin/R/qc_posteriors.R" "$posteriors" "$input" "chr${chr}"
        else
          touch ${chr}_noplot.png 
        fi
        """
}

// Prepare a to-score format using variant map to align with the genotype snp ids
process format_sbayesr_posteriors {
    publishDir "${params.outdir}/intermediates", mode: 'copy', overwrite: true

    input:
        tuple val(chr), path(posterior), path(map), path(map_noNA)
    
    output:
        tuple val(chr), path("${chr}_posterior_scoreformat")
    script:
        snp_posteriors_cols="2,5,8"
        // select ss_ID and bim_ID
        map_from_to="3,6"
        """
        format_posteriors.sh ${posterior} ${snp_posteriors_cols} ${map} ${map_from_to} "true" > "${chr}_posterior_scoreformat"
        """
}


