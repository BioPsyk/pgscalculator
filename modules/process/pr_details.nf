process copyConfigFiles {
    publishDir "${params.outdir}/details", mode: 'copy', overwrite: true
    
    input:
    path("nextflow.config")

    output:
    path("nextflow.config")
    path("*.config")

    script:
    """
    mkdir -p ${params.outdir}/details
    cp nextflow.config ${params.outdir}/details/
    cp ${params.conffile} ${params.outdir}/details/
    # Also copy files to output for Nextflow
    cp nextflow.config .
    cp ${params.conffile} .
    """
}

process copyVersionFile {
    publishDir "${params.outdir}/details", mode: 'copy', overwrite: true

    input:
    path("VERSION")

    output:
    path("VERSION")

    script:
    """
    mkdir -p ${params.outdir}/details
    cp VERSION ${params.outdir}/details/
    # Also copy to output for Nextflow
    cp VERSION .
    """
}


