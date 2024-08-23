process copyConfigFiles {
    
    input:
    path("nextflow.config")

    script:
    """
    mkdir -p ${params.outdir}/details
    cp nextflow.config ${params.outdir}/details/
    cp ${params.conffile} ${params.outdir}/details/

    """
}


