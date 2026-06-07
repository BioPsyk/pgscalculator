#!/bin/bash
# Test scenarios for pgscalculator v2.1 output implementation
#
# Usage:
#   ./test_scenarios.sh [scenario]
#   
#   scenario: 1, 2, or 3 (or "all" to run all)
#
# Scenarios:
#   1. No INFO/MAF files - compute MAF from genotypes
#   2. With mock INFO file - apply INFO filter
#   3. With mock MAF file - use provided MAF

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGS_DIR="$(dirname "$SCRIPT_DIR")"
WRAPPER="${PGS_DIR}/pgscalculator.sh"

# Base paths (adjust as needed)
INPUT_DIR="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/sumstats/5668"
LDDIR="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/ld-sbayesr/ukb/band_ukb_10k_hm3"
GENODIR="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/plink2"
GENOFILE="/faststorage/project/ibp_migration_opengdk/PROJECT_pgscalculator/pgscalculator/references/genotypes_test/mapfiles/plink2_genodir_genofiles.txt"

SCENARIO="${1:-1}"

echo "=========================================="
echo "pgscalculator v2.1 - Test Scenarios"
echo "=========================================="

run_test_1() {
    echo ""
    echo "=== Scenario 1: No INFO/MAF files ==="
    echo "Expected: MAF computed from genotypes, INFO filter skipped"
    
    local OUTDIR="${SCRIPT_DIR}/out_test_scenario1"
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"
    
    # Create config WITHOUT info_file or maf_file
    cat > "${OUTDIR}/config.yaml" << EOF
# Scenario 1: No reference files
# MAF will be computed from genotypes
# INFO filter will be skipped

ld_reference: ${LDDIR}
genotypes: ${GENODIR}
genotype_manifest: ${GENOFILE}
outdir: ${OUTDIR}

filters:
  info_threshold: 0.8
  maf_threshold: 0.01

whichn: totalN
chromosomes: 21-22

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 4
  seed: 12345
  exclude_mhc: true

score_columns: 1 2 5

slurm:
  account: ibp_pipeline_pgscalculator
  partition: normal
  prep: { mem: 16g, cpus: 4, time: '02:00:00' }
  posteriors: { mem: 20g, cpus: 6, time: '04:00:00' }
  score: { mem: 10g, cpus: 4, time: '01:00:00' }
EOF
    
    echo "Config created: ${OUTDIR}/config.yaml"
    echo ""
    echo "To run prep step:"
    echo "  ${WRAPPER} --config ${OUTDIR}/config.yaml --steps prep --sbatch"
    echo ""
    echo "Or without SLURM (interactive):"
    echo "  ${WRAPPER} --config ${OUTDIR}/config.yaml --steps prep"
}

run_test_2() {
    echo ""
    echo "=== Scenario 2: With INFO file ==="
    echo "Expected: INFO filter applied, MAF computed from genotypes"
    
    local OUTDIR="${SCRIPT_DIR}/out_test_scenario2"
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"
    
    # Create mock INFO file
    echo "Creating mock INFO file..."
    
    mkdir -p "${OUTDIR}/references"
    cat > "${OUTDIR}/references/info_scores.tsv" << 'EOF'
GENO_ID	INFO
rs7287144	0.95
rs4010554	0.85
rs4010558	0.75
rs2379981	0.92
rs9605903	0.65
rs16980739	0.88
EOF
    # Note: rs4010558 and rs9605903 have INFO < 0.8, should be filtered
    
    # Create config WITH info_file
    cat > "${OUTDIR}/config.yaml" << EOF
# Scenario 2: With INFO file
# INFO filter will be applied
# MAF computed from genotypes

ld_reference: ${LDDIR}
genotypes: ${GENODIR}
genotype_manifest: ${GENOFILE}
outdir: ${OUTDIR}

references:
  info_file: ${OUTDIR}/references/info_scores.tsv

filters:
  info_threshold: 0.8
  maf_threshold: 0.01

whichn: totalN
chromosomes: 21-22

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 4
  seed: 12345
  exclude_mhc: true

score_columns: 1 2 5

slurm:
  account: ibp_pipeline_pgscalculator
  partition: normal
  prep: { mem: 16g, cpus: 4, time: '02:00:00' }
EOF
    
    echo "Config created: ${OUTDIR}/config.yaml"
    echo "INFO file: ${OUTDIR}/references/info_scores.tsv"
    echo ""
    echo "To run:"
    echo "  ${WRAPPER} --config ${OUTDIR}/config.yaml --steps prep --sbatch"
}

run_test_3() {
    echo ""
    echo "=== Scenario 3: With MAF file ==="
    echo "Expected: MAF filter uses provided file (not computed)"
    
    local OUTDIR="${SCRIPT_DIR}/out_test_scenario3"
    rm -rf "$OUTDIR"
    mkdir -p "$OUTDIR"
    
    # Create mock MAF file
    echo "Creating mock MAF file..."
    
    mkdir -p "${OUTDIR}/references"
    cat > "${OUTDIR}/references/maf.tsv" << 'EOF'
GENO_ID	MAF
rs7287144	0.15
rs4010554	0.02
rs4010558	0.005
rs2379981	0.25
rs9605903	0.18
rs16980739	0.03
EOF
    # Note: rs4010558 has MAF < 0.01, should be filtered
    
    # Create config WITH maf_file
    cat > "${OUTDIR}/config.yaml" << EOF
# Scenario 3: With MAF file
# MAF filter uses provided file (not computed from genotypes)

ld_reference: ${LDDIR}
genotypes: ${GENODIR}
genotype_manifest: ${GENOFILE}
outdir: ${OUTDIR}

references:
  maf_file: ${OUTDIR}/references/maf.tsv

filters:
  info_threshold: 0.8
  maf_threshold: 0.01

whichn: totalN
chromosomes: 21-22

sbayesr:
  gamma: "0.0,0.01,0.1,1"
  pi: "0.95,0.02,0.02,0.01"
  burn_in: 2000
  chain_length: 10000
  threads: 4
  seed: 12345
  exclude_mhc: true

score_columns: 1 2 5

slurm:
  account: ibp_pipeline_pgscalculator
  partition: normal
  prep: { mem: 16g, cpus: 4, time: '02:00:00' }
EOF
    
    echo "Config created: ${OUTDIR}/config.yaml"
    echo "MAF file: ${OUTDIR}/references/maf.tsv"
    echo ""
    echo "To run:"
    echo "  ${WRAPPER} --config ${OUTDIR}/config.yaml --steps prep --sbatch"
}

case "$SCENARIO" in
    1) run_test_1 ;;
    2) run_test_2 ;;
    3) run_test_3 ;;
    all)
        run_test_1
        run_test_2
        run_test_3
        ;;
    *)
        echo "Usage: $0 [1|2|3|all]"
        echo ""
        echo "Scenarios:"
        echo "  1 - No INFO/MAF files (compute MAF from genotypes)"
        echo "  2 - With INFO file"
        echo "  3 - With MAF file"
        echo "  all - Set up all scenarios"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Test setup complete!"
echo ""
echo "Next steps:"
echo "  1. Submit prep job for your scenario"
echo "  2. After prep completes, run: --steps sumstat,posteriors,score"
echo "  3. Check outputs in out_test_scenario*/sumstats/*/"
echo "=========================================="

