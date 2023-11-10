process copyConfigFiles {
    
    input:
    path("nextflow.config")

    script:
    """
    mkdir -p ${params.outdir}/details
    cp nextflow.config ${params.outdir}/details/
    cp /pgscalculator/confdir/* ${params.outdir}/details/
    """
}

