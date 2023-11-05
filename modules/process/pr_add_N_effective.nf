
process pr_add_N_effective {
    input:
    path(sfile)
    path(metafile)

    output:
    path("sfile_added_N.gz")

    script:
    """
    add_N_to_sumstat.sh ${sfile} ${metafile} | gzip -c > sfile_added_N.gz
    """
}

