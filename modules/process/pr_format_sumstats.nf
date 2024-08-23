// Nextflow processes format of a gwas cleansumstats default output 
process add_build_sumstats {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
        path(input_file)
        path(input_map)

    output:
        path('added_sumstat_grch37')

    script:
        """
        # Add grch37 chr and pos as new column 1 and 2
        add_build_sumstats.sh ${input_file} ${input_map} "added_sumstat_grch37"
        """
}
process filter_NA_coordinates {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
        tuple val(chr), path(file1)

    output:
        tuple val(chr), path("${chr}_filtered_on_NA")

    script:
        """
        filter_NA_coordinates.sh ${file1} > "${chr}_filtered_on_NA"
        """
}
process rmcol_build_sumstats {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     

    input:
        tuple val(chr), path(file1), path(file2)
        val(torm)

    output:
        tuple val(chr), path("${chr}_rmcol_sumstat_grch37_map"), path("${chr}_rmcol_sumstat_grch37_map_noNA")

    script:
        """
        # Remove one of two builds 1 (col 1 and 2) or 2 (col 3 and 4)
        rmcol_build_sumstats.sh ${file1} ${torm} "${chr}_rmcol_sumstat_grch37_map"
        rmcol_build_sumstats.sh ${file2} ${torm} "${chr}_rmcol_sumstat_grch37_map_noNA"
        """
}

process change_build_sumstats {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     

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
     
    input:
    tuple val(chr), path(input_file), path(input_file_noNA)
    file(mapfile)
    val(method)

    output:
    tuple val(chr), path("${chr}_formatted")

    script:
        """
        format_sumstats.sh ${input_file_noNA} ${mapfile} ${method} > "${chr}_formatted"
        """
}

process add_N_effective {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    
    tuple val(chr), path(sfile), path(metafile)
    val(whichN)

    output:
    tuple val(chr), path("${chr}_sfile_added_N"), path(metafile)

    script:
    """
    add_N_to_sumstat.sh ${sfile} ${metafile} ${whichN} > ${chr}_sfile_added_N
    """
}

process force_EAF_to_sumstat {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile), path(metafile)

    output:
    tuple val(chr), path("${chr}_sfile_forced_EAF")

    script:
    """
    force_EAF_to_sumstat.sh ${sfile} ${metafile} > ${chr}_sfile_forced_EAF
    """
}

process add_B_and_SE {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile)

    output:
    tuple val(chr), path("${chr}_sfile_added_B_SE")

    script:
    """
    fill_in_beta_and_se.sh ${sfile} > ${chr}_sfile_added_B_SE
    """
}

process filter_bad_values_1 {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile)

    output:
    tuple val(chr), path("${chr}_sfile_filter_bad_values_1")

    script:
    """
    filter_bad_values.sh ${sfile} > ${chr}_sfile_filter_bad_values_1
    """
}
process filter_bad_values_2 {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    tuple val(chr), path(sfile)

    output:
    tuple val(chr), path("${chr}_sfile_filter_bad_values_2")

    script:
    """
    filter_bad_values.sh ${sfile} > ${chr}_sfile_filter_bad_values_2
    """
}

process filter_on_ldref_rsids {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)
    path(rsids)

    output:
    path("sfile_filter_on_rsids")

    script:
    """
    filter_on_ldref_rsids.sh ${sfile} ${rsids} > sfile_filter_on_rsids
    """
}

process split_on_chromosome {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
    input:
    path(sfile)

    output:
    path('split*')

    script:
    """
    split_on_chromosome.sh ${sfile} "CHR" "splitss" "zcat"
    """
}

process concatenate_sumstat_input {
    publishDir "${params.outdir}/intermediates", mode: 'copy', overwrite: true

    input:
        path(chrinput)

    output:
        tuple val("all"), path("allchr_input")
    script:
        """
        arr=(${chrinput})
        head -n1 "\${arr[1]}" > "allchr_input"
        for chrfile in ${chrinput}
        do
          tail -n+2 \$chrfile >> "allchr_input"
        done
        """
}

