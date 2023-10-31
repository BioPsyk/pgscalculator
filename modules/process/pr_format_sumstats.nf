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
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     
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

