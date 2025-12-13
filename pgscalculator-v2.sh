#!/usr/bin/env bash

# pgscalculator v2.1.0 wrapper script
# Config-first approach: paths in config.yaml, minimal CLI

################################################################################
# Help page
################################################################################

function general_usage(){
  echo "Usage:"
  echo "  ./pgscalculator-v2.sh --config <file> [options]"
  echo ""
  echo "Required:"
  echo "  --config <file>   Path to config.yaml with all settings"
  echo ""
  echo "Optional:"
  echo "  -i <dir>          Path to sumstats folder (overrides config)"
  echo "  -o <dir>          Path to output directory (overrides config)"
  echo "  --sumstat <name>  Sumstat name/ID (extracted from -i path if not provided)"
  echo "  --steps <list>    Steps to run: prep, sumstat, posteriors, score (default: all)"
  echo "  --skip-prep       Skip prep steps if already completed"
  echo "  --chr <range>     Chromosomes to process (e.g., '21-22', default: 1-22)"
  echo "  -d                Dev mode (verbose output)"
  echo "  -v                Show version"
  echo "  -h                Show this help"
  echo ""
  echo "Config file (config.yaml) should contain:"
  echo "  ld_reference: /path/to/band_ukb_10k_hm3"
  echo "  genotypes: /path/to/genotypes"
  echo "  genotype_manifest: /path/to/manifest.txt"
  echo "  outdir: /path/to/output  # optional, can use -o instead"
  echo "  # Plus sbayesr parameters, thresholds, etc."
  echo ""
  echo "Examples:"
  echo "  # Run all steps for a sumstat"
  echo "  ./pgscalculator-v2.sh --config config.yaml -i /path/to/sumstat_814"
  echo ""
  echo "  # Run only prep (reusable across sumstats)"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep"
  echo ""
  echo "  # Run per-sumstat steps (after prep is done)"
  echo "  ./pgscalculator-v2.sh --config config.yaml -i /path/to/sumstat_814 --skip-prep"
  echo ""
  echo "  # Run specific steps only"
  echo "  ./pgscalculator-v2.sh --config config.yaml -i /path/to/sumstat_814 --steps posteriors,score --skip-prep"
}

################################################################################
# Prepare path parsing
################################################################################
present_dir="${PWD}"
project_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

################################################################################
# Parameter parsing
################################################################################
paramarray=($@)

# Parse all arguments
config_file=""
infold=""
outdir=""
steps_arg=""
sumstat_name=""
skip_prep=false
chromosomes=""
devmode=""

i=0
while [ $i -lt ${#paramarray[@]} ]; do
  case "${paramarray[$i]}" in
    --config)
      config_file="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    --steps)
      steps_arg="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    --sumstat)
      sumstat_name="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    --skip-prep)
      skip_prep=true
      i=$((i+1))
      ;;
    --chr)
      chromosomes="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    -i)
      infold="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    -o)
      outdir="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    -d)
      devmode="--verbose"
      i=$((i+1))
      ;;
    -v)
      cat ${project_dir}/VERSION 1>&2
      exit 0
      ;;
    -h|--help)
      general_usage 1>&2
      exit 0
      ;;
    *)
      echo "Unknown option: ${paramarray[$i]}" 1>&2
      general_usage 1>&2
      exit 1
      ;;
  esac
done

################################################################################
# Validate config file
################################################################################
if [[ -z "$config_file" ]]; then
  >&2 echo "Error: --config is required"
  >&2 echo "Run with -h for help"
  exit 1
fi

if [[ ! -f "$config_file" ]]; then
  >&2 echo "Error: Config file not found: $config_file"
  exit 1
fi

config_file_host=$(realpath "$config_file")

################################################################################
# Parse config file (simple YAML parsing with awk)
################################################################################
parse_yaml_value() {
  local key="$1"
  local file="$2"
  awk -F': ' -v key="$key" '$1 == key {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$file"
}

# Read paths from config
cfg_ld_reference=$(parse_yaml_value "ld_reference" "$config_file_host")
cfg_genotypes=$(parse_yaml_value "genotypes" "$config_file_host")
cfg_genotype_manifest=$(parse_yaml_value "genotype_manifest" "$config_file_host")
cfg_outdir=$(parse_yaml_value "outdir" "$config_file_host")
cfg_chromosomes=$(parse_yaml_value "chromosomes" "$config_file_host")

# CLI overrides config
if [[ -n "$outdir" ]]; then
  cfg_outdir="$outdir"
fi

if [[ -n "$chromosomes" ]]; then
  cfg_chromosomes="$chromosomes"
fi

################################################################################
# Validate required paths
################################################################################
# LD reference is required
if [[ -z "$cfg_ld_reference" ]]; then
  >&2 echo "Error: ld_reference not found in config file"
  exit 1
fi

if [[ ! -d "$cfg_ld_reference" ]]; then
  >&2 echo "Error: LD reference directory not found: $cfg_ld_reference"
  exit 1
fi

# Genotypes required for scoring
if [[ -z "$cfg_genotypes" ]] || [[ -z "$cfg_genotype_manifest" ]]; then
  >&2 echo "Warning: genotypes and genotype_manifest should be in config for scoring"
fi

# Output directory required (from config or CLI)
if [[ -z "$cfg_outdir" ]]; then
  >&2 echo "Error: outdir not found in config file and -o not provided"
  exit 1
fi

# Input sumstat: required for non-prep steps
if [[ -z "$infold" ]] && [[ "$steps_arg" != "prep" ]]; then
  # Check if we're only running prep
  if [[ -z "$steps_arg" ]] || [[ "$steps_arg" == *"sumstat"* ]] || [[ "$steps_arg" == *"posteriors"* ]] || [[ "$steps_arg" == *"score"* ]]; then
    >&2 echo "Error: -i (sumstat input) is required for non-prep steps"
    exit 1
  fi
fi

################################################################################
# Resolve paths
################################################################################
mkdir -p "${cfg_outdir}"
outdir_host=$(realpath "${cfg_outdir}")
lddir_host=$(realpath "${cfg_ld_reference}")

if [[ -n "$cfg_genotypes" ]] && [[ -d "$cfg_genotypes" ]]; then
  genodir_host=$(realpath "${cfg_genotypes}")
else
  genodir_host=""
fi

if [[ -n "$cfg_genotype_manifest" ]] && [[ -f "$cfg_genotype_manifest" ]]; then
  genofile_host=$(realpath "${cfg_genotype_manifest}")
else
  genofile_host=""
fi

if [[ -n "$infold" ]]; then
  infold_host=$(realpath "${infold}")
  if [[ ! -d "$infold_host" ]]; then
    >&2 echo "Error: Input directory doesn't exist: $infold_host"
    exit 1
  fi
else
  # For prep-only, use a placeholder
  infold_host="${outdir_host}"
fi

# Extract sumstat name from input path if not provided
if [[ -z "$sumstat_name" ]] && [[ -n "$infold" ]]; then
  sumstat_name=$(basename "$infold_host" | sed 's/^sumstat_//')
  if [[ "$sumstat_name" == "$(basename "$infold_host")" ]]; then
    sumstat_name=$(basename "$infold_host")
  fi
fi

################################################################################
# Determine steps to run
################################################################################
if [[ -n "$steps_arg" ]]; then
  run_all=false
else
  run_all=true
fi

################################################################################
# Prepare container variables
################################################################################
source "${project_dir}/scripts/init-containerization.sh"

# Default to singularity
mountflag="-B"

# Container paths
indir_container="/pgscalculator/input"
foldername=$(basename "$lddir_host")
lddir_container="/pgscalculator/$foldername"
genodir_container="/pgscalculator/genodir"
outdir_container="/pgscalculator/outdir"
config_container="/pgscalculator/config"

if [[ -n "$genofile_host" ]]; then
  genodir2_host=$(dirname "${genofile_host}")
  genofile_name=$(basename "${genofile_host}")
  genodir2_container="/pgscalculator/genodir2"
  genofile_container="${genodir2_container}/${genofile_name}"
else
  genodir2_host=""
  genofile_container=""
fi

################################################################################
# Generate container config.yaml
################################################################################
config_yaml_host="${outdir_host}/config.yaml"
config_yaml_container="${outdir_container}/config.yaml"

# Copy original config and update paths for container
cat > "${config_yaml_host}" << EOF
# pgscalculator v2.1.0 - Auto-generated config for container
# Generated from: ${config_file_host}
input: ${indir_container}
outdir: ${outdir_container}
lddir: ${lddir_container}
genodir: ${genodir_container}
genofile: ${genofile_container}
EOF

# Copy parameters from user's config (skip path keys we've already set)
awk '
  !/^(ld_reference|genotypes|genotype_manifest|outdir|input|lddir|genodir|genofile):/ {
    print
  }
' "$config_file_host" >> "${config_yaml_host}"

# Add/override chromosome range if specified
if [[ -n "$cfg_chromosomes" ]]; then
  # Remove existing chromosomes line and add new one
  sed -i '/^chromosomes:/d' "${config_yaml_host}"
  echo "chromosomes: ${cfg_chromosomes}" >> "${config_yaml_host}"
fi

################################################################################
# Determine image
################################################################################
runimage="sif/${singularity_image_tag}"
if [[ ! -f "${project_dir}/${runimage}" ]]; then
  # Try to find any available sif file
  runimage=$(ls -1 "${project_dir}"/sif/*.sif 2>/dev/null | head -1)
  if [[ -z "$runimage" ]]; then
    >&2 echo "Error: No singularity image found in ${project_dir}/sif/"
    exit 1
  fi
fi

source "${project_dir}/conf/init-docker-config.sh"

################################################################################
# Build command
################################################################################
# Build mount flags
mount_flags=$(format_mount_flags "${mountflag}")

# Build CLI command
if [[ "$run_all" == true ]]; then
  cli_cmd="/pgscalculator/bin/pgscalculator run --all --sumstat ${sumstat_name} --config ${config_yaml_container}"
else
  cli_cmd="/pgscalculator/bin/pgscalculator run --steps ${steps_arg} --sumstat ${sumstat_name} --config ${config_yaml_container}"
  if [[ "$skip_prep" == true ]]; then
    cli_cmd="${cli_cmd} --skip-prep"
  fi
fi

if [[ -n "$devmode" ]]; then
  cli_cmd="${cli_cmd} ${devmode}"
fi

################################################################################
# Build mount options
################################################################################
mount_opts=""
mount_opts="${mount_opts} ${mountflag} ${outdir_host}:${outdir_container}"
mount_opts="${mount_opts} ${mountflag} ${lddir_host}:${lddir_container}"

if [[ -n "$infold" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${infold_host}:${indir_container}"
fi

if [[ -n "$genodir_host" ]] && [[ -d "$genodir_host" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${genodir_host}:${genodir_container}"
fi

if [[ -n "$genodir2_host" ]] && [[ -d "$genodir2_host" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${genodir2_host}:${genodir2_container}"
fi

################################################################################
# Execute
################################################################################
echo "Running pgscalculator v2.1.0 in Singularity"
echo "Config: ${config_file_host}"
echo "Output: ${outdir_host}"
echo "Command: ${cli_cmd}"

singularity run \
  --contain \
  --cleanenv \
  ${mount_flags} \
  ${mount_opts} \
  "${runimage}" \
  ${cli_cmd}
