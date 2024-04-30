
process extract_metadata_from_sumstat {
    label 'low_mem'

    input:
    path(metafile)

    output:
    env(neff) 

    script:
    """
    neff=\$(awk -F ': ' '\$1=="stats_EffectiveN"{print \$2}' ${metafile} )
    """
}

