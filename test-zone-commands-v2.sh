#!/bin/bash
# pgscalculator v2.0.0 - Test commands script
# Compatible with existing commands.sh structure

# Configuration (same as original commands.sh)
inlist="short_list.txt"
outdir="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/test-zone/out_test_v2"
indir="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.6.7"

pgsfold="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator"

mkdir -p ${outdir}

while read id ; do
  infold="${indir}/sumstat_${id}"
  sleep 0.1
  
  # Example: Run all steps
  sbatch --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
    --account ibp_pipeline_cleansumstats \
    --job-name="pgs_v2_${id}" \
    --output="pgs_v2_${id}.out" \
    --error="pgs_v2_${id}.err" \
    --wrap="
  echo ${id}
  date
  /bin/bash ${pgsfold}/pgscalculator-v2.sh \
    -i ${infold} \
    -l ${pgsfold}/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
    -c ${pgsfold}/conf/sbayesr.config \
    -o ${outdir}/sumstat_${id} \
    -g ${pgsfold}/references/genotypes_test/plink \
    -f ${pgsfold}/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt
  echo ${id}
  date
  "
  
  # Example: Run only prep steps (once per project)
  # Uncomment to run prep separately:
  # sbatch --mem=10g --ntasks 1 --cpus-per-task 6 --time=0:30:00 \
  #   --account ibp_pipeline_cleansumstats \
  #   --job-name="prep_${id}" \
  #   --output="prep_${id}.out" \
  #   --error="prep_${id}.err" \
  #   --wrap="
  # /bin/bash ${pgsfold}/pgscalculator-v2.sh \
  #   -i ${infold} \
  #   -l ${pgsfold}/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  #   -c ${pgsfold}/conf/sbayesr.config \
  #   -o ${outdir}/sumstat_${id} \
  #   -g ${pgsfold}/references/genotypes_test/plink \
  #   -f ${pgsfold}/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
  #   --steps prep
  # "
  
  # Example: Run only posteriors (skip prep if already done)
  # Uncomment to run posteriors separately:
  # sbatch --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  #   --account ibp_pipeline_cleansumstats \
  #   --job-name="post_${id}" \
  #   --output="post_${id}.out" \
  #   --error="post_${id}.err" \
  #   --wrap="
  # /bin/bash ${pgsfold}/pgscalculator-v2.sh \
  #   -i ${infold} \
  #   -l ${pgsfold}/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
  #   -c ${pgsfold}/conf/sbayesr.config \
  #   -o ${outdir}/sumstat_${id} \
  #   --steps posteriors --skip-prep
  # "

done < ${inlist}


