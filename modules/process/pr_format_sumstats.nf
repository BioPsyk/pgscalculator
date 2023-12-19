// Nextflow processes format of a gwas cleansumstats default output 
process change_build_sumstats {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
        path(input_file)
        path(input_map)

    output:
        path('clean_sumstat_grch37')

    script:
        """
        change_build_sumstats.sh ${input_file} ${input_map} "clean_sumstat_grch37"
        """
}

process format_sumstats {
    publishDir "${params.outdir}/intermediates/format_sumstats", mode: 'rellink', overwrite: true, enabled: params.dev
     
    input:
    tuple val(chr), path(input_file)
    file(mapfile)
    val(method)

    output:
    tuple val(chr), path('formatted')

    script:
        """
        format_sumstats.sh ${input_file} ${mapfile} ${method} > formatted
        """
}

process add_N_effective {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    
    tuple val(chr), path(sfile), path(metafile)
    val(whichN)

    output:
    tuple val(chr), path("${chr}_sfile_added_N"), path(metafile)

    script:
    """
    add_N_to_sumstat.sh ${sfile} ${metafile} ${whichN} > ${chr}_sfile_added_N
    """
}

process force_EAF_to_sumstat {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile), path(metafile)

    output:
    tuple val(chr), path("${chr}_sfile_forced_EAF")

    script:
    """
    force_EAF_to_sumstat.sh ${sfile} ${metafile} > ${chr}_sfile_forced_EAF
    """
}

process add_B_and_SE {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile)

    output:
    tuple val(chr), path("${chr}_sfile_added_B_SE")

    script:
    """
    fill_in_beta_and_se.sh ${sfile} > ${chr}_sfile_added_B_SE
    """
}

process filter_bad_values {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile)

    output:
    tuple val(chr), path("${chr}_sfile_filter_bad_values")

    script:
    """
    filter_bad_values.sh ${sfile} > ${chr}_sfile_filter_bad_values
    """
}

process filter_on_ldref_rsids {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)
    path(rsids)

    output:
    path("sfile_filter_on_rsids")

    script:
    """
    filter_on_ldref_rsids.sh ${sfile} ${rsids} > sfile_filter_on_rsids
    """
}

process split_on_chromosome {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)

    output:
    path('split*')

    script:
    """
    split_on_chromosome.sh ${sfile} "split"
    """
}

