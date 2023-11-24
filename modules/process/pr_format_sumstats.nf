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
     
    label 'low_mem'

    input:
        path(input_file)
        file(mapfile)
        val(method)

    output:
        path('formatted_*')

    script:
        """
        format_sumstats.sh ${input_file} ${mapfile} ${method} formatted
        """
}

process add_N_effective {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)
    path(metafile)
    val(whichN)

    output:
    path("sfile_added_N.gz")

    script:
    """
    add_N_to_sumstat.sh ${sfile} ${metafile} ${whichN} | gzip -c > sfile_added_N.gz
    """
}

process force_EAF_to_sumstat {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)
    path(metafile)

    output:
    path("sfile_forced_EAF.gz")

    script:
    """
    force_EAF_to_sumstat.sh ${sfile} ${metafile} | gzip -c > sfile_forced_EAF.gz
    """
}

process add_B_and_SE {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)

    output:
    path("sfile_added_N.gz")

    script:
    """
    add_B_and_SE.sh ${sfile} | gzip -c > sfile_added_N.gz
    """
}




