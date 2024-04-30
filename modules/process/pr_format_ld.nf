#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process split_on_chromosome_prscs_ld {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    label 'low_mem'

    input:
    path(ldfile)

    output:
    path('split*')

    script:
    """
    cp ${ldfile} ld
    split_on_chromosome.sh ${ldfile} "CHR" "splitld" "cat"
    """
}
