process copyConfigFiles {
    publishDir "${params.outdir}/details", mode: 'copy', overwrite: true
    
    input:
    path(config_file)

    output:
    path("nextflow.config")
    path("user.config")

    script:
    """
    # Copy input config file and user config to output
    cp ${config_file} nextflow.config
    cp ${params.conffile} user.config
    """
}

process copyVersionFile {
    publishDir "${params.outdir}/details", mode: 'copy', overwrite: true

    input:
    path(version_file)

    output:
    path("VERSION")

    script:
    """
    # Copy version file to output
    cp ${version_file} VERSION
    """
}


