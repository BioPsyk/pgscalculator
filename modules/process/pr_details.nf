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
    if [ "${config_file}" != "nextflow.config" ]; then
        cp ${config_file} nextflow.config
    else
        # File is already named correctly, just touch to ensure it exists
        touch nextflow.config
    fi
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


