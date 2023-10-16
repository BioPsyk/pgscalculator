// Nextflow process to format a gwas cleansumstats default output 

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
        format_sumstats.sh $input_file ${mapfile} ${method} formatted
        """
}

