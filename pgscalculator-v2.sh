#!/usr/bin/env bash

# pgscalculator v2.0.0 wrapper script
# Compatible with v1 command-line interface, runs new modular CLI

################################################################################
# Help page
################################################################################

function general_usage(){
 echo "Usage:"
 echo " ./pgscalculator-v2.sh -i <dir> -o <dir> -l <dir> -c <file> [options]"
 echo ""
 echo "options:"
 echo "-h          Display help message"
 echo "-i <dir>    Path to sumstats folder (cleansumstats output)"
 echo "-l <dir>    LD map dir, absolute paths"
 echo "-g <dir>    Target genotypes directory"
 echo "-f <file>   Genotype manifest file"
 echo "-s <file>   SNP list filtering (default: none)"
 echo "-c <file>   Config file (sbayesr.config or prscs.config)"
 echo "-o <dir>    Path to output directory"
 echo "-b <dir>    Path to system tmp or scratch (default: /tmp)"
 echo "-w <dir>    Path to workdir/intermediate files (default: work)"
 echo "-j <mode>   Image mode: docker, dockerhub_biopsyk, or singularity (default: singularity)"
 echo "-d          Dev mode, keep intermediates"
 echo "-v          Get version number"
 echo "-1          Disable step1 (calc posteriors) - only format sumstat"
 echo "-2          Disable step2 (calc score) - only calc posteriors"
 echo "--steps <steps>  Specify steps to run (e.g., prep,posteriors,score or all)"
 echo "--sumstat <name> Sumstat name/ID (extracted from -i path if not provided)"
 echo ""
 echo "Examples:"
 echo "  # Run all steps"
 echo "  ./pgscalculator-v2.sh -i /path/to/sumstat_123 -o /path/to/out -l /path/to/ld -c conf/sbayesr.config"
 echo ""
 echo "  # Run only prep and posteriors"
 echo "  ./pgscalculator-v2.sh -i /path/to/sumstat_123 -o /path/to/out -l /path/to/ld -c conf/sbayesr.config --steps prep,posteriors"
 echo ""
 echo "  # Run only scoring (skip prep if already done)"
 echo "  ./pgscalculator-v2.sh -i /path/to/sumstat_123 -o /path/to/out -l /path/to/ld -c conf/sbayesr.config --steps score --skip-prep"
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

# Parse new-style arguments first
steps_arg=""
sumstat_name=""
skip_prep=false

# Remove --steps, --sumstat, --skip-prep from array for getopts
new_paramarray=()
i=0
while [ $i -lt ${#paramarray[@]} ]; do
    if [[ "${paramarray[$i]}" == "--steps" ]]; then
        steps_arg="${paramarray[$((i+1))]}"
        i=$((i+2))
    elif [[ "${paramarray[$i]}" == "--sumstat" ]]; then
        sumstat_name="${paramarray[$((i+1))]}"
        i=$((i+2))
    elif [[ "${paramarray[$i]}" == "--skip-prep" ]]; then
        skip_prep=true
        i=$((i+1))
    else
        new_paramarray+=("${paramarray[$i]}")
        i=$((i+1))
    fi
done

# Defaults
infold=""
lddir=""
genodir=""
genofile=""
snpfile=""
conffile=""
outdir="out"
container_image=""
calc_posterior=true
calc_score=true

infold_given=false
lddir_given=false
genodir_given=false
genofile_given=false
snpfile_given=false
conffile_given=false
outdir_given=false
tmpdir_given=false
devmode_given=false
calc_posterior_given=false
calc_score_given=false
container_image_given=false

tmpdir="/tmp"
workdir="${present_dir}/work"
devmode=""

getoptsstring=":hvi:o:b:w:l:g:f:s:m:c:db:12j:"

while getopts "${getoptsstring}" opt "${new_paramarray[@]}"; do
  case ${opt} in
    h )
      general_usage 1>&2
      exit 0
      ;;
    v )
      cat ${project_dir}/VERSION 1>&2
      exit 0
      ;;
    i )
      infold="$OPTARG"
      infold_given=true
      ;;
    l )
      lddir="$OPTARG"
      lddir_given=true
      ;;
    g )
      genodir="$OPTARG"
      genodir_given=true
      ;;
    f )
      genofile="$OPTARG"
      genofile_given=true
      ;;
    s )
      snpfile="$OPTARG"
      snpfile_given=true
      ;;
    c )
      conffile="$OPTARG"
      conffile_given=true
      ;;
    o )
      outdir="$OPTARG"
      outdir_given=true
      ;;
    j )
      container_image="$OPTARG"
      container_image_given=true
      ;;
    b )
      tmpdir="$OPTARG"
      tmpdir_given=true
      ;;
    w )
      workdir="$OPTARG"
      workdir_given=true
      ;;
    d )
      devmode="--dev"
      devmode_given=true
      ;;
    1 )
      calc_posterior=false
      calc_posterior_given=true
      ;;
    2 )
      calc_score=false
      calc_score_given=true
      ;;
    \? )
      echo "Invalid Option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Invalid Option: -$OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done

################################################################################
# Validate required arguments
################################################################################
if ! ${infold_given}; then
  >&2 echo "Error: -i (input directory) is required"
  exit 1
fi

if ! ${outdir_given}; then
  >&2 echo "Error: -o (output directory) is required"
  exit 1
fi

if ! ${lddir_given}; then
  >&2 echo "Error: -l (LD directory) is required"
  exit 1
fi

if ! ${conffile_given}; then
  >&2 echo "Error: -c (config file) is required"
  exit 1
fi

if ! ${genodir_given} || ! ${genofile_given}; then
  >&2 echo "Warning: -g and -f (genotype directory and manifest) are required for scoring"
fi

################################################################################
# Resolve paths
################################################################################
mkdir -p ${outdir}
mkdir -p ${workdir}
mkdir -p ${tmpdir}

infold_host=$(realpath "${infold}")
outdir_host=$(realpath "${outdir}")
tmpdir_host=$(realpath "${tmpdir}")
workdir_host=$(realpath "${workdir}")

if [ ! -d "$infold_host" ]; then
  >&2 echo "Error: Input directory doesn't exist: $infold_host"
  exit 1
fi

lddir_host=$(realpath "${lddir}")
if [ ! -d "$lddir_host" ]; then
  >&2 echo "Error: LD directory doesn't exist: $lddir_host"
  exit 1
fi

if ${genodir_given}; then
  genodir_host=$(realpath "${genodir}")
  if [ ! -d "$genodir_host" ]; then
    >&2 echo "Error: Genotype directory doesn't exist: $genodir_host"
    exit 1
  fi
fi

if ${genofile_given}; then
  genofile_host=$(realpath "${genofile}")
  if [ ! -f "$genofile_host" ]; then
    >&2 echo "Error: Genotype manifest file doesn't exist: $genofile_host"
    exit 1
  fi
fi

conffile_host=$(realpath "${conffile}")
if [ ! -f "$conffile_host" ]; then
  >&2 echo "Error: Config file doesn't exist: $conffile_host"
  exit 1
fi

# Extract sumstat name from input path if not provided
if [[ -z "$sumstat_name" ]]; then
  sumstat_name=$(basename "$infold_host" | sed 's/^sumstat_//')
  if [[ "$sumstat_name" == "$(basename "$infold_host")" ]]; then
    # If no sumstat_ prefix, use the directory name as-is
    sumstat_name=$(basename "$infold_host")
  fi
fi

################################################################################
# Determine steps to run
################################################################################
if [[ -n "$steps_arg" ]]; then
  # User specified steps explicitly
  run_all=false
elif [[ "$calc_posterior" == false ]] && [[ "$calc_score" == false ]]; then
  # Both disabled - only format sumstat
  steps_arg="sumstat"
  run_all=false
elif [[ "$calc_posterior" == false ]]; then
  # Only posteriors disabled - run prep, sumstat, score
  steps_arg="prep,sumstat,score"
  run_all=false
elif [[ "$calc_score" == false ]]; then
  # Only score disabled - run prep, sumstat, posteriors
  steps_arg="prep,sumstat,posteriors"
  run_all=false
else
  # Default: run all steps
  run_all=true
fi

################################################################################
# Prepare container variables
################################################################################
source "${project_dir}/scripts/init-containerization.sh"

# Which mount symbol to use
if [ "${container_image}" == "docker" ] || [ "${container_image}" == "dockerhub_biopsyk" ]; then
  mountflag="-v"
else
  mountflag="-B"
fi

# Container paths (same as v1)
indir_container="/pgscalculator/input"
foldername=$(basename "$lddir_host")
lddir_container="/pgscalculator/$foldername"
genodir_container="/pgscalculator/genodir"
genodir2_host=$(dirname "${genofile_host}")
genofile_name=$(basename "${genofile_host}")
genodir2_container="/pgscalculator/genodir2"
genofile_container="${genodir2_container}/${genofile_name}"
confdir_host=$(dirname "${conffile_host}")
conffile_name=$(basename "${conffile_host}")
confdir_container="/pgscalculator/confdir"
conffile_container="${confdir_container}/${conffile_name}"
outdir_container="/pgscalculator/outdir"
tmpdir_container="/tmp"
workdir_container="/pgscalculator/work"

if ${snpfile_given}; then
  snpdir_host=$(dirname "${snpfile_host}")
  snpfile_name=$(basename "${snpfile_host}")
  snpdir_container="/pgscalculator/snpdir"
  snpfile_container="${snpdir_container}/${snpfile_name}"
  snplist_host_container="${mountflag} ${snpdir_host}:${snpdir_container}"
else
  snplist_host_container=""
fi

# Create config.yaml in output directory
config_yaml_host="${outdir_host}/config.yaml"
config_yaml_container="${outdir_container}/config.yaml"

# Generate config.yaml
cat > "${config_yaml_host}" << EOF
# pgscalculator v2.0.0 - Auto-generated config
input: ${indir_container}
outdir: ${outdir_container}
lddir: ${lddir_container}
genodir: ${genodir_container}
genofile: ${genofile_container}
EOF

# Add optional parameters from config file if it's a yaml file
if [[ "$conffile_name" == *.yaml ]] || [[ "$conffile_name" == *.yml ]]; then
  # If config is already yaml, we can source some values
  # For now, just add basic sbayesr defaults
  cat >> "${config_yaml_host}" << EOF
info_threshold: 0.8
maf_threshold: 0.01
whichn: totalN
sbayesr:
  gamma: 0.0,0.01,0.1,1
  pi: 0.95,0.02,0.02,0.01
  burn_in: 2000
  chain_length: 10000
  threads: 6
  seed: 80851
  exclude_mhc: true
score_columns: 2 5 9
EOF
else
  # For .config files, add defaults (user can override later)
  cat >> "${config_yaml_host}" << EOF
info_threshold: 0.8
maf_threshold: 0.01
whichn: totalN
sbayesr:
  gamma: 0.0,0.01,0.1,1
  pi: 0.95,0.02,0.02,0.01
  burn_in: 2000
  chain_length: 10000
  threads: 6
  seed: 80851
  exclude_mhc: true
score_columns: 2 5 9
EOF
fi

################################################################################
# Determine image
################################################################################
if [ "${container_image}" == "docker" ]; then
  runimage="${image_tag}" 
elif [ "${container_image}" == "dockerhub_biopsyk" ]; then
  runimage="${deploy_image_tag_docker_hub}" 
elif [ "${container_image}" == "" ]; then
  runimage="sif/${singularity_image_tag}" 
else
  runimage="${container_image}" 
fi

source "${project_dir}/conf/init-docker-config.sh"

################################################################################
# Build command
################################################################################
# Build mount flags
mount_flags=$(format_mount_flags "${mountflag}")

# Build CLI command
if [[ "$run_all" == true ]]; then
  cli_cmd="pgscalculator run --all --sumstat ${sumstat_name} --config ${config_yaml_container}"
else
  cli_cmd="pgscalculator run --steps ${steps_arg} --sumstat ${sumstat_name} --config ${config_yaml_container}"
  if [[ "$skip_prep" == true ]]; then
    cli_cmd="${cli_cmd} --skip-prep"
  fi
fi

if [[ -n "$devmode" ]]; then
  cli_cmd="${cli_cmd} --verbose"
fi

################################################################################
# Execute
################################################################################
if [ "${container_image}" == "docker" ] || [ "${container_image}" == "dockerhub_biopsyk" ]; then
  echo "Running pgscalculator v2.0.0 in Docker"
  echo "Command: ${cli_cmd}"
  exec docker run \
     --rm \
     ${docker_run_args} \
     ${mount_flags} \
     ${mountflag} "${infold_host}:${indir_container}" \
     ${mountflag} "${outdir_host}:${outdir_container}" \
     ${mountflag} "${lddir_host}:${lddir_container}" \
     ${mountflag} "${genodir_host}:${genodir_container}" \
     ${mountflag} "${genodir2_host}:${genodir2_container}" \
     ${mountflag} "${confdir_host}:${confdir_container}" \
     ${mountflag} "${tmpdir_host}:${tmpdir_container}" \
     ${mountflag} "${workdir_host}:${workdir_container}" \
     ${snplist_host_container} \
     "${runimage}" \
     ${cli_cmd}
else
  echo "Running pgscalculator v2.0.0 in Singularity"
  echo "Command: ${cli_cmd}"
  singularity run \
     --contain \
     --cleanenv \
     ${mount_flags} \
     ${mountflag} "${infold_host}:${indir_container}" \
     ${mountflag} "${outdir_host}:${outdir_container}" \
     ${mountflag} "${lddir_host}:${lddir_container}" \
     ${mountflag} "${genodir_host}:${genodir_container}" \
     ${mountflag} "${genodir2_host}:${genodir2_container}" \
     ${mountflag} "${confdir_host}:${confdir_container}" \
     ${mountflag} "${tmpdir_host}:${tmpdir_container}" \
     ${mountflag} "${workdir_host}:${workdir_container}" \
     ${snplist_host_container} \
     "${runimage}" \
     ${cli_cmd}
fi


