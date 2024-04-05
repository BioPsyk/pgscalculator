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
    tuple val(chr), path("${chr}_variants_mapfile"), emit: map
    tuple val(chr), path("${chr}_variants_mapfile_no_NA"), emit: map_noNA

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
        awk -vFS="\t" -vOFS="\t" '{if (\$9 != "NA") { print \$0 }}' "${chr}_variants_mapfile" > "${chr}_variants_mapfile_no_NA"
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

process prepare_augmentad_sumstat_files {
    publishDir "${params.outdir}/intermediates/mapgeneration", mode: 'rellink', overwrite: true, enabled: params.dev
     
    label 'low_mem'

    input:
    tuple val(chr), path(ss), path(map)

    output:
    tuple val(chr), path("${chr}_augmented_ss")

    script:
        """
        # subset on chr:pos:a1:a2
        subset_sumstat_on_mapfile.sh ${ss} ${map} "${chr}_subset_sumstats"
        """
}

