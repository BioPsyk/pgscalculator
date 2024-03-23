// Nextflow processes to make initial variant map

process variant_map_for_prscs {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
    tuple val(chr), path(ss), path(bim), path(ldmeta)

    output:
    tuple val(chr), path("${chr}_variants_mapfile")

    script:
        """

        # b37 b38 a1 a2 markername
        # chr:pos chr:pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" 'NR>1{print \$1":"\$2, \$3":"\$4, \$7, \$8, \$6 }' ${ss} > ss2
        # chr pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" '{print \$1":"\$4, \$5, \$6, \$2 }' ${bim} > bim2
        # chr pos a1 a2 markername
        awk -vFS=" " -vOFS="\t" 'NR>1{print \$1":"\$3, \$4, \$5, \$2 }' ${ldmeta} > ld2

        # use same variant map as for sbayesr
        variant_map_for_sbayesr.sh ss2 bim2 ld2 "${chr}_variants_mapfile"
        """
}

process variant_map_for_sbayesr {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
    tuple val(chr), path(ss), path(bim), path(ldbin), path(ldinfo)

    output:
    tuple val(chr), path("${chr}_variants_mapfile")

    script:
        """

        # b37 b38 a1 a2 markername
        # chr:pos chr:pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" 'NR>1{print \$1":"\$2, \$3":"\$4, \$7, \$8, \$6 }' ${ss} > ss2
        # chr pos a1 a2 markername
        awk -vFS="\t" -vOFS="\t" '{print \$1":"\$4, \$5, \$6, \$2 }' ${bim} > bim2
        # chr pos a1 a2 markername
        awk -vFS=" " -vOFS="\t" 'NR>1{print \$1":"\$4, \$5, \$6, \$2 }' ${ldinfo} > ld2

        variant_map_for_sbayesr.sh ss2 bim2 ld2 "${chr}_variants_mapfile"
        """
}

process filter_sumstat_variants_on_map_file {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
    tuple val(chr), path(ss), path(map)

    output:
    tuple val(chr), path("${chr}_subset_sumstats")

    script:
        """
        # subset on chr:pos:a1:a2
        subset_sumstat_on_mapfile.sh ${ss} ${map} "${chr}_subset_sumstats"
        """
}

// Deprecated
process concat_bim_files {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
        path(chrinput)

    output:
        path('allchr_bim')

    script:
        """
        arr=(${chrinput})
        head -n1 "\${arr[1]}" > "allchr_bim"
        for chrfile in ${chrinput}
        do
          tail -n+2 \$chrfile >> "allchr_bim"
        done
        """
}

// Deprecated
process subset_on_bim_file {
    publishDir "${params.outdir}/intermediates", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
        tuple val(chr), path(sumstat), path(bim)
      

    output:
        tuple val(chr), path("${chr}_sumstat_subset")

    script:
        """
        subset_sumstat_on_bim.sh sumstat bim ${chr}_sumstat_subset
        """
}

