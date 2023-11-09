process copyConfigFiles {
    
    input:
    path 'nextflow.config'
    path 'base.config'

    script:
    """
    mkdir -p ${params.outdir}/details
    cp nextflow.config ${params.outdir}/details/
    cp base.config ${params.outdir}/details/base.config
    """
}

