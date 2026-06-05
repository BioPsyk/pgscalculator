#!/bin/bash
# pgscalculator v2.0.0 - Test Commands
# Adapted from test-zone/2024-11-28-commands.sh for v2 testing
#
# This file contains different ways of running tests for pgscalculator v2.
# Copy the relevant section and run manually, or use as a reference.
#
# Created: 2024-12-11
# Based on: test-zone/2024-11-28-commands.sh

# =============================================================================
# Common paths (update these for your environment)
# =============================================================================

# Project location
pgsfold="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator"

# Sumstats (cleansumstats output) - use version_1.12.0 or latest
indir="/faststorage/project/ibp_pipeline_cleansumstats/raw_library/sumstat_clean_library/version_1.12.0"

# Reduced sumstat for quick testing
indir_reduced="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/test-zone/reduced_sumstat"

# LD reference
ldref="${pgsfold}/references/ld-sbayesr/ukb/band_ukb_10k_hm3"

# Test genotypes (plink format)
genodir="${pgsfold}/references/genotypes_test/plink"
genofile="${pgsfold}/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt"

# Test genotypes (plink2 format)
genodir_plink2="${pgsfold}/references/genotypes_test/plink2"
genofile_plink2="${pgsfold}/references/genotypes_test/mapfiles/plink2_genodir_genofiles.txt"

# Config file
configfile="${pgsfold}/config.template.yaml"

# Test sumstat ID (example)
id=814


# =============================================================================
# TEST 1: Normal Run (Full Pipeline)
# =============================================================================
# This runs all steps: prep, format-sumstat, posteriors, score

outdir="${pgsfold}/test-v2/out_test_normal"
mkdir -p ${outdir}

infold="${indir}/sumstat_${id}"

sbatch --mem=40g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --job-name="pgs_v2_${id}" \
  --output="${outdir}/pgs_v2_${id}.out" \
  --error="${outdir}/pgs_v2_${id}.err" \
  --wrap="
echo 'Starting pgscalculator v2 - Full Pipeline'
echo 'Sumstat ID: ${id}'
date

/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -g ${genodir} \
  -f ${genofile} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  -d

echo 'Completed: ${id}'
date
"


# =============================================================================
# TEST 2: Reduced Sumstat (Quick Test)
# =============================================================================
# Run with reduced set of variants in sumstat
# (only 30 on each chromosome, at least 1 million bp apart)
# This is faster for debugging

outdir="${pgsfold}/test-v2/out_test_reduced"
mkdir -p ${outdir}

infold="${indir_reduced}/sumstat_${id}"

sbatch --mem=40g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --job-name="pgs_v2_red_${id}" \
  --output="${outdir}/pgs_v2_red_${id}.out" \
  --error="${outdir}/pgs_v2_red_${id}.err" \
  --wrap="
echo 'Starting pgscalculator v2 - Reduced Sumstat'
echo 'Sumstat ID: ${id}'
date

/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -g ${genodir} \
  -f ${genofile} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  -d

echo 'Completed: ${id}'
date
"


# =============================================================================
# TEST 3: Posteriors Only (Skip Scoring)
# =============================================================================
# Run only posterior calculation, skip scoring step

outdir="${pgsfold}/test-v2/out_test_posteriors_only"
mkdir -p ${outdir}

infold="${indir}/sumstat_${id}"

sbatch --mem=40g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --job-name="pgs_v2_post_${id}" \
  --output="${outdir}/pgs_v2_post_${id}.out" \
  --error="${outdir}/pgs_v2_post_${id}.err" \
  --wrap="
echo 'Starting pgscalculator v2 - Posteriors Only'
echo 'Sumstat ID: ${id}'
date

/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -g ${genodir} \
  -f ${genofile} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  -d \
  -2

echo 'Completed: ${id}'
date
"


# =============================================================================
# TEST 4: Step-by-Step Execution (NEW in v2)
# =============================================================================
# Run specific steps using --steps flag

outdir="${pgsfold}/test-v2/out_test_stepwise"
mkdir -p ${outdir}

infold="${indir}/sumstat_${id}"

# Step 1: Prep only (reusable across sumstats)
sbatch --mem=10g --ntasks 1 --cpus-per-task 6 --time=0:30:00 \
  --account ibp_pipeline_cleansumstats \
  --job-name="pgs_v2_prep_${id}" \
  --output="${outdir}/pgs_v2_prep_${id}.out" \
  --error="${outdir}/pgs_v2_prep_${id}.err" \
  --wrap="
echo 'Starting pgscalculator v2 - Prep Steps Only'
date

/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -g ${genodir} \
  -f ${genofile} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  --steps prep

echo 'Prep completed'
date
"

# Step 2: Posteriors only (after prep is done)
sbatch --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
  --account ibp_pipeline_cleansumstats \
  --job-name="pgs_v2_post2_${id}" \
  --output="${outdir}/pgs_v2_post2_${id}.out" \
  --error="${outdir}/pgs_v2_post2_${id}.err" \
  --dependency=singleton \
  --wrap="
echo 'Starting pgscalculator v2 - Posteriors (skip prep)'
echo 'Sumstat ID: ${id}'
date

/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  --steps posteriors --skip-prep

echo 'Posteriors completed'
date
"

# Step 3: Scoring only (after posteriors is done)
sbatch --mem=10g --ntasks 1 --cpus-per-task 6 --time=0:30:00 \
  --account ibp_pipeline_cleansumstats \
  --job-name="pgs_v2_score_${id}" \
  --output="${outdir}/pgs_v2_score_${id}.out" \
  --error="${outdir}/pgs_v2_score_${id}.err" \
  --dependency=singleton \
  --wrap="
echo 'Starting pgscalculator v2 - Scoring (skip prep)'
echo 'Sumstat ID: ${id}'
date

/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -g ${genodir} \
  -f ${genofile} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  --steps score --skip-prep

echo 'Scoring completed'
date
"


# =============================================================================
# TEST 5: Interactive Development (Inside Container)
# =============================================================================
# For debugging and development - run interactively inside container

# First, request an interactive node:
srun --mem=40g --ntasks 1 --cpus-per-task 8 --time=2:00:00 \
  --account ibp_pipeline_cleansumstats \
  --pty /bin/bash

# Then set up variables:
outdir="${pgsfold}/test-v2/out_test_dev"
mkdir -p ${outdir}
id=814
infold="${indir}/sumstat_${id}"

# Run wrapper script directly (for debugging):
/bin/bash ${pgsfold}/pgscalculator.sh \
  -i ${infold} \
  -l ${ldref} \
  -g ${genodir_plink2} \
  -f ${genofile_plink2} \
  -c ${configfile} \
  -o ${outdir}/sumstat_${id} \
  -d

# Or enter container for interactive debugging:
singularity shell --contain --cleanenv \
  -B /faststorage:/faststorage \
  ${pgsfold}/sif/ibp-pgscalculator-base_version-2.0.0.sif

# Inside container, use the CLI directly:
# pgscalculator status --config /path/to/config.yaml
# pgscalculator prep-genotypes --config /path/to/config.yaml
# pgscalculator run --all --sumstat TRAIT --config /path/to/config.yaml


# =============================================================================
# TEST 6: Batch Processing (Multiple Sumstats)
# =============================================================================
# Process multiple sumstats from a list

outdir="${pgsfold}/test-v2/out_test_batch"
mkdir -p ${outdir}

# Create or use an existing list of sumstat IDs
inlist="${pgsfold}/test-v2/short_list.txt"

# Example list (create this file):
# echo "814" > ${inlist}
# echo "815" >> ${inlist}
# echo "816" >> ${inlist}

while read id ; do
  infold="${indir}/sumstat_${id}"
  sleep 0.1
  
  sbatch --mem=20g --ntasks 1 --cpus-per-task 22 --time=1:00:00 \
    --account ibp_pipeline_cleansumstats \
    --job-name="pgs_v2_batch_${id}" \
    --output="${outdir}/pgs_v2_batch_${id}.out" \
    --error="${outdir}/pgs_v2_batch_${id}.err" \
    --wrap="
  echo 'Batch: ${id}'
  date
  /bin/bash ${pgsfold}/pgscalculator.sh \
    -i ${infold} \
    -l ${ldref} \
    -g ${genodir} \
    -f ${genofile} \
    -c ${configfile} \
    -o ${outdir}/sumstat_${id}
  echo 'Completed: ${id}'
  date
  "
done < ${inlist}


# =============================================================================
# Original v1 command (for reference)
# =============================================================================
# This is what the old pgscalculator.sh looked like:
#
# /bin/bash ${pgsfold}/pgscalculator.sh \
#   -i ${infold} \
#   -l ${pgsfold}/references/ld-sbayesr/ukb/band_ukb_10k_hm3 \
#   -g ${pgsfold}/references/genotypes_test/plink \
#   -f ${pgsfold}/references/genotypes_test/mapfiles/plink_genodir_genofiles.txt \
#   -c ${pgsfold}/config.template.yaml \
#   -o ${outdir}/sumstat_${id} \
#   -d




