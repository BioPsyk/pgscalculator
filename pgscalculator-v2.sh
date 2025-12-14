#!/usr/bin/env bash

# pgscalculator v2.1.0 wrapper script
# Config-first approach: paths in config.yaml, minimal CLI

################################################################################
# Help page
################################################################################

function general_usage(){
  echo "Usage:"
  echo "  ./pgscalculator-v2.sh --config <file> --steps <steps> [options]"
  echo ""
  echo "Required:"
  echo "  --config <file>   Path to config.yaml with all settings"
  echo "  --steps <list>    Steps to run: prep, sumstat, posteriors, score"
  echo ""
  echo "Optional:"
  echo "  -i <dir>          Path to sumstats folder (required for non-prep steps)"
  echo "  -o <dir>          Path to output directory (overrides config)"
  echo "  --chr <range>     Chromosomes to process (e.g., '21-22', default: 1-22)"
  echo "  --sbatch          Submit as SLURM job using sbatch settings from config"
  echo "  -d                Dev mode (verbose output)"
  echo "  -v                Show version"
  echo "  -h                Show this help"
  echo ""
  echo "Step groups:"
  echo "  prep        Run prep steps (genotypes, ldref, inclusion-list)"
  echo "  sumstat     Format and filter sumstat"
  echo "  posteriors  Calculate posteriors with sbayesR"
  echo "  score       Calculate PGS scores"
  echo ""
  echo "Config file (config.yaml) should contain:"
  echo "  ld_reference: /path/to/band_ukb_10k_hm3"
  echo "  genotypes: /path/to/genotypes"
  echo "  genotype_manifest: /path/to/manifest.txt"
  echo "  outdir: /path/to/output"
  echo ""
  echo "  # Optional: SLURM settings for --sbatch"
  echo "  slurm:"
  echo "    account: my_account"
  echo "    prep:       { mem: 10g, cpus: 6, time: '1:00:00' }"
  echo "    posteriors: { mem: 20g, cpus: 8, time: '2:00:00' }"
  echo "    score:      { mem: 10g, cpus: 4, time: '0:30:00' }"
  echo ""
  echo "Examples:"
  echo "  # Step 1: Run prep (once per project)"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep"
  echo ""
  echo "  # Step 2: Run per-sumstat steps"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_814"
  echo ""
  echo "  # Or submit as SLURM jobs"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat_814 --sbatch"
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
chromosomes=""
devmode=""
use_sbatch=false

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
    --chr)
      chromosomes="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    --sbatch)
      use_sbatch=true
      i=$((i+1))
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
# Validate required arguments
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

if [[ -z "$steps_arg" ]]; then
  >&2 echo "Error: --steps is required"
  >&2 echo ""
  >&2 echo "Available steps:"
  >&2 echo "  prep        - Prepare genotypes and LD reference (run once)"
  >&2 echo "  sumstat     - Format and filter sumstat"
  >&2 echo "  posteriors  - Calculate posteriors with sbayesR"
  >&2 echo "  score       - Calculate PGS scores"
  >&2 echo ""
  >&2 echo "Examples:"
  >&2 echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep"
  >&2 echo "  ./pgscalculator-v2.sh --config config.yaml --steps sumstat,posteriors,score -i /path/to/sumstat"
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

# Parse nested YAML values (e.g., slurm.account, slurm.prep.mem)
parse_yaml_nested() {
  local section="$1"
  local key="$2"
  local file="$3"
  awk -v section="$section" -v key="$key" '
    BEGIN { in_section = 0 }
    /^[a-zA-Z]/ { in_section = 0 }
    $0 ~ "^"section":" { in_section = 1; next }
    in_section && $0 ~ "^  "key":" {
      gsub(/^  [a-zA-Z_]+: */, "")
      gsub(/[{}]/, "")
      print
      exit
    }
  ' "$file"
}

# Parse inline YAML dict (e.g., "{ mem: 10g, cpus: 6, time: '1:00:00' }")
parse_inline_dict() {
  local dict="$1"
  local key="$2"
  echo "$dict" | sed 's/[{}]//g' | tr ',' '\n' | awk -F': ' -v key="$key" '
    $1 ~ key { gsub(/^[ \t]+|[ \t]+$|'"'"'/, "", $2); print $2 }
  '
}

# Read paths from config
cfg_ld_reference=$(parse_yaml_value "ld_reference" "$config_file_host")
cfg_genotypes=$(parse_yaml_value "genotypes" "$config_file_host")
cfg_genotype_manifest=$(parse_yaml_value "genotype_manifest" "$config_file_host")
cfg_outdir=$(parse_yaml_value "outdir" "$config_file_host")
cfg_chromosomes=$(parse_yaml_value "chromosomes" "$config_file_host")

# Read optional reference files (for INFO/MAF filtering)
cfg_info_file=$(parse_yaml_nested "references" "info_file" "$config_file_host")
cfg_maf_file=$(parse_yaml_nested "references" "maf_file" "$config_file_host")

# CLI overrides config
if [[ -n "$outdir" ]]; then
  cfg_outdir="$outdir"
fi

if [[ -n "$chromosomes" ]]; then
  cfg_chromosomes="$chromosomes"
fi

################################################################################
# Handle --sbatch: submit as SLURM job
################################################################################
if [[ "$use_sbatch" == true ]]; then
  # Read SLURM settings from config
  slurm_account=$(parse_yaml_nested "slurm" "account" "$config_file_host")
  slurm_partition=$(parse_yaml_nested "slurm" "partition" "$config_file_host")
  
  # Determine which step profile to use
  step_profile="default"
  if [[ -n "$steps_arg" ]]; then
    # Use first step as profile (prep, posteriors, score, etc.)
    step_profile=$(echo "$steps_arg" | cut -d',' -f1)
  fi
  
  # Get step-specific settings
  step_settings=$(parse_yaml_nested "slurm" "$step_profile" "$config_file_host")
  
  # Parse settings or use defaults
  if [[ -n "$step_settings" ]]; then
    slurm_mem=$(parse_inline_dict "$step_settings" "mem")
    slurm_cpus=$(parse_inline_dict "$step_settings" "cpus")
    slurm_time=$(parse_inline_dict "$step_settings" "time")
  fi
  
  # Apply defaults if not set
  slurm_mem="${slurm_mem:-20g}"
  slurm_cpus="${slurm_cpus:-8}"
  slurm_time="${slurm_time:-2:00:00}"
  
  # Build job name
  if [[ -n "$infold" ]]; then
    job_name="pgs_$(basename "$infold" | sed 's/^sumstat_//')"
  else
    job_name="pgs_${step_profile}"
  fi
  
  # Build the command to run (same command without --sbatch)
  run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps ${steps_arg}"
  [[ -n "$infold" ]] && run_cmd="${run_cmd} -i ${infold}"
  [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
  [[ -n "$chromosomes" ]] && run_cmd="${run_cmd} --chr ${chromosomes}"
  [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"
  
  # Determine output directory for logs
  log_dir="${cfg_outdir:-./}"
  mkdir -p "$log_dir"
  
  # Build sbatch command
  sbatch_cmd="sbatch"
  sbatch_cmd="${sbatch_cmd} --mem=${slurm_mem}"
  sbatch_cmd="${sbatch_cmd} --cpus-per-task=${slurm_cpus}"
  sbatch_cmd="${sbatch_cmd} --time=${slurm_time}"
  sbatch_cmd="${sbatch_cmd} --job-name=${job_name}"
  sbatch_cmd="${sbatch_cmd} --output=${log_dir}/${job_name}.out"
  sbatch_cmd="${sbatch_cmd} --error=${log_dir}/${job_name}.err"
  [[ -n "$slurm_account" ]] && sbatch_cmd="${sbatch_cmd} --account=${slurm_account}"
  [[ -n "$slurm_partition" ]] && sbatch_cmd="${sbatch_cmd} --partition=${slurm_partition}"
  sbatch_cmd="${sbatch_cmd} --wrap=\"${run_cmd}\""
  
  echo "Submitting SLURM job..."
  echo "  Job name: ${job_name}"
  echo "  Resources: mem=${slurm_mem}, cpus=${slurm_cpus}, time=${slurm_time}"
  echo "  Logs: ${log_dir}/${job_name}.out"
  echo "  Command: ${run_cmd}"
  echo ""
  
  # Submit and exit
  eval ${sbatch_cmd}
  exit $?
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

# Extract sumstat name from input path
sumstat_name=""
if [[ -n "$infold" ]]; then
  sumstat_name=$(basename "$infold_host" | sed 's/^sumstat_//')
fi

################################################################################
# Prerequisite checking
################################################################################
check_prep_exists() {
  local outdir="$1"
  local missing=""
  
  # Check for genotype prep outputs
  if [[ ! -d "${outdir}/prep/genotypes" ]] || [[ -z "$(ls -A "${outdir}/prep/genotypes" 2>/dev/null)" ]]; then
    missing="${missing}  - Genotype prep: ${outdir}/prep/genotypes/\n"
  fi
  
  # Check for inclusion list
  if [[ ! -f "${outdir}/prep/inclusion-list/inclusion_list.txt" ]]; then
    missing="${missing}  - Inclusion list: ${outdir}/prep/inclusion-list/inclusion_list.txt\n"
  fi
  
  if [[ -n "$missing" ]]; then
    echo "$missing"
    return 1
  fi
  return 0
}

check_sumstat_exists() {
  local outdir="$1"
  local sumstat="$2"
  
  if [[ ! -f "${outdir}/sumstats/${sumstat}/formatted/cleaned_sumstat.tsv" ]]; then
    echo "  - Formatted sumstat: ${outdir}/sumstats/${sumstat}/formatted/cleaned_sumstat.tsv"
    return 1
  fi
  return 0
}

check_posteriors_exists() {
  local outdir="$1"
  local sumstat="$2"
  
  if [[ ! -d "${outdir}/sumstats/${sumstat}/posteriors" ]] || [[ -z "$(ls -A "${outdir}/sumstats/${sumstat}/posteriors" 2>/dev/null)" ]]; then
    echo "  - Posteriors: ${outdir}/sumstats/${sumstat}/posteriors/"
    return 1
  fi
  return 0
}

# Check prerequisites based on requested steps
if [[ "$steps_arg" != "prep" ]]; then
  # Non-prep steps require prep to be completed
  missing_prep=$(check_prep_exists "$outdir_host")
  if [[ $? -ne 0 ]]; then
    >&2 echo "Error: Prep outputs not found."
    >&2 echo ""
    >&2 echo "Missing:"
    >&2 echo -e "$missing_prep"
    >&2 echo ""
    >&2 echo "Run prep first:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps prep"
    >&2 echo ""
    >&2 echo "Then retry your command."
    exit 1
  fi
fi

# Check sumstat prerequisite for posteriors
if [[ "$steps_arg" == *"posteriors"* ]] && [[ "$steps_arg" != *"sumstat"* ]]; then
  missing_sumstat=$(check_sumstat_exists "$outdir_host" "$sumstat_name")
  if [[ $? -ne 0 ]]; then
    >&2 echo "Error: Sumstat formatting not completed."
    >&2 echo ""
    >&2 echo "Missing:"
    >&2 echo "$missing_sumstat"
    >&2 echo ""
    >&2 echo "Run sumstat step first:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps sumstat -i ${infold}"
    >&2 echo ""
    >&2 echo "Or include sumstat in your steps:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps sumstat,posteriors -i ${infold}"
    exit 1
  fi
fi

# Check posteriors prerequisite for score
if [[ "$steps_arg" == *"score"* ]] && [[ "$steps_arg" != *"posteriors"* ]]; then
  missing_posteriors=$(check_posteriors_exists "$outdir_host" "$sumstat_name")
  if [[ $? -ne 0 ]]; then
    >&2 echo "Error: Posteriors calculation not completed."
    >&2 echo ""
    >&2 echo "Missing:"
    >&2 echo "$missing_posteriors"
    >&2 echo ""
    >&2 echo "Run posteriors step first:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps posteriors -i ${infold}"
    >&2 echo ""
    >&2 echo "Or include posteriors in your steps:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps posteriors,score -i ${infold}"
    exit 1
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
# Use a different filename to avoid overwriting user's original config
config_yaml_host="${outdir_host}/config_container.yaml"
config_yaml_container="${outdir_container}/config_container.yaml"

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

# Copy parameters from user's config (skip path keys and references section - handled separately)
awk '
  BEGIN { in_references = 0 }
  /^references:/ { in_references = 1; next }
  /^[a-zA-Z]/ && in_references { in_references = 0 }
  in_references { next }
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
  cli_cmd="/pgscalculator/bin/pgscalculator run --all --config ${config_yaml_container}"
  if [[ -n "$sumstat_name" ]]; then
    cli_cmd="${cli_cmd} --sumstat ${sumstat_name}"
  fi
else
  cli_cmd="/pgscalculator/bin/pgscalculator run --steps ${steps_arg} --config ${config_yaml_container}"
  if [[ -n "$sumstat_name" ]]; then
    cli_cmd="${cli_cmd} --sumstat ${sumstat_name}"
  fi
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

# Mount INFO file if provided
if [[ -n "$cfg_info_file" ]] && [[ -f "$cfg_info_file" ]]; then
  info_file_host=$(realpath "$cfg_info_file")
  info_file_container="/pgscalculator/references/info_scores.tsv"
  mount_opts="${mount_opts} ${mountflag} ${info_file_host}:${info_file_container}"
  # Add to config
  echo "info_file: ${info_file_container}" >> "${config_yaml_host}"
fi

# Mount MAF file if provided
if [[ -n "$cfg_maf_file" ]] && [[ -f "$cfg_maf_file" ]]; then
  maf_file_host=$(realpath "$cfg_maf_file")
  maf_file_container="/pgscalculator/references/maf.tsv"
  mount_opts="${mount_opts} ${mountflag} ${maf_file_host}:${maf_file_container}"
  # Add to config
  echo "maf_file: ${maf_file_container}" >> "${config_yaml_host}"
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
