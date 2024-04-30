process make_augmented_gwas {
    publishDir "${params.outdir}/intermediates/augmented_sumstat", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'low_mem'

    input:
    tuple val(chr), path(sumstat_map), path(sumstat_map_noNA), path(maffile), path(posteriors), path(benchmark), path(map),path(map_noNA)

    output:
    tuple val(chr), path("${chr}_sumstat_augmented")

    script:
    """
  #  cat ${sumstat_map} > ss2
  #  cat ${map} > map2
  #  cat ${maffile} > maf2
  #  cat ${posteriors} > post2
  #  cat ${benchmark} > bench2
    make_augmented_output.sh ${sumstat_map} ${map} ${maffile} ${posteriors} ${benchmark} > "${chr}_sumstat_augmented"
    """
}

// Concatenate per chromosome posterior SNP effects for sBayesR
process concatenate_augmented_sumstat {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    label 'low_mem'

    input:
        path(chrfiles)
    
    output:
        path("augmented_sumstat.gz")
    script:
        """
        echo "RSID	CHR	POS	EffectAllele	EffectAllele	B	SE	Z	P	genoID	MAF	postEffect	benchEffect"  > "augmented_sumstat"
        for chrfile in ${chrfiles}
        do
          tail -n+2 \$chrfile >> "augmented_sumstat"
        done

        gzip "augmented_sumstat"
        """
}
