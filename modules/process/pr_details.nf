process copyConfigFiles {
    publishDir "${params.outdir}/details", mode: 'copy', overwrite: true
    
    input:
    path(config_file)

    output:
    path("pipeline_nextflow.config")
    path("user_pipeline.config")

    script:
    """
    # Copy input config file and user config to output with different names
    cp ${config_file} pipeline_nextflow.config
    cp ${params.conffile} user_pipeline.config
    """
}

process copyVersionFile {
    publishDir "${params.outdir}/details", mode: 'copy', overwrite: true

    input:
    path(version_file)

    output:
    path("pipeline_VERSION")

    script:
    """
    # Copy version file to output with different name
    cp ${version_file} pipeline_VERSION
    """
}


