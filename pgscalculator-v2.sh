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
  echo "  --steps <list>    Steps to run: prep, sumstat, weights, score, finalize"
  echo ""
  echo "Optional:"
  echo "  -i <dir>          Path to sumstats folder (required for non-prep steps)"
  echo "  -o <dir>          Path to output directory (overrides config)"
  echo "  --methods <list>  Active posterior methods: sbayesr, ldpred2, or both"
  echo "                   (comma-separated; overrides config methods:)"
  echo "  --force           Force re-run of steps even if already completed"
  echo "  --sbatch          Submit as SLURM job using sbatch settings from config"
  echo "                   For per-sumstat steps (sumstat/weights/score), this submits a"
  echo "                   lightweight *driver job* that runs sumstat and launches/monitors"
  echo "                   chromosome-parallel arrays per active method (see --methods)."
  echo "  -d                Dev mode (verbose output)"
  echo "  --cleanup         After a successful run/driver job, remove per-run work/ and tmp/ folders"
  echo "                   (Default for now: keep work/tmp, which is useful during development)"
  echo "  -v                Show version"
  echo "  -h                Show this help"
  echo ""
  echo "Step groups:"
  echo "  prep        Run prep steps (genotypes, ldref, inclusion-list)"
  echo "  sumstat     Format and filter sumstat"
  echo "  weights     Posterior weights per active method (sBayesR array + LDpred2 single job + benchmark)"
  echo "  score       Calculate PGS scores (one array per active method)"
  echo ""
 echo "Config file (config.yaml) should contain:"
 echo "  lddir: /path/to/band_ukb_10k_hm3"
 echo "  genodir: /path/to/genotypes"
 echo "  genofile: /path/to/manifest.txt"
 echo "  outdir: /path/to/output"
 echo ""
  echo "  # Optional: SLURM settings for --sbatch (weights = two array jobs)"
  echo "  slurm:"
  echo "    account: my_account"
  echo "    partition: normal"
  echo "    driver:            { mem: 1g, cpus: 1, time: '2:00:00' }"
  echo "    prep:              { mem: 10g, cpus: 1, time: '1:00:00', max_parallel: 22 }"
  echo "    sumstat:          { mem: 1g, cpus: 1, time: '0:30:00', max_parallel: 22 }"
    echo "    weights_sbayesr:   { mem: 20g, cpus: 6, time: '2:00:00', max_parallel: 22 }"
    echo "    weights_ldpred2:   { mem: 64g, cpus: 16, time: '4:00:00' }"
    echo "    weights_benchmark: { mem: 2g, cpus: 2, time: '0:30:00', max_parallel: 22 }"
    echo "    score_sbayesr:     { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }"
    echo "    score_ldpred2:     { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }"
    echo "    score:             { mem: 10g, cpus: 4, time: '0:30:00', max_parallel: 22 }  # inherits to score_*"
  echo "  benchmark:"
  echo "    maf_threshold: 0.05"
  echo "    indep_pairwise: [250, 50, 0.25]"
 echo ""
 echo "Examples:"
  echo "  # Step 1: Run prep (once per project)"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep"
 echo ""
  echo "  # Step 2: Run per-sumstat steps"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_814"
 echo ""
  echo "  # Or submit as SLURM jobs"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep --sbatch"
  echo "  ./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat_814 --sbatch"
}

################################################################################
# Prepare path parsing
################################################################################
present_dir="${PWD}"
project_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${project_dir}/bin/lib/common.sh"

################################################################################
# Parameter parsing
################################################################################
paramarray=($@)

# Parse all arguments
config_file=""
infold=""
outdir=""
steps_arg=""
chromosomes_override=""
methods_cli=""
method_override=""
devmode=""
force_mode=""
use_sbatch=false
driver_run=false
do_cleanup=false

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
    --_chr)
        chromosomes_override="${paramarray[$((i+1))]}"
        i=$((i+2))
      ;;
    --sbatch)
      use_sbatch=true
        i=$((i+1))
      ;;
    --_driver-run)
      # Internal flag: run inside a driver job; orchestrates arrays, does not submit itself.
      driver_run=true
        i=$((i+1))
      ;;
    --cleanup)
      do_cleanup=true
      i=$((i+1))
      ;;
    --force|-f)
      force_mode="--force"
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
    --methods)
      methods_cli="${paramarray[$((i+1))]}"
      i=$((i+2))
      ;;
    --method)
      method_override="${paramarray[$((i+1))]}"
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
  >&2 echo "  weights     - sBayesR posteriors + benchmark weights (two array jobs)"
  >&2 echo "  score       - Calculate PGS scores"
  >&2 echo ""
  >&2 echo "Examples:"
  >&2 echo "  ./pgscalculator-v2.sh --config config.yaml --steps prep"
  >&2 echo "  ./pgscalculator-v2.sh --config config.yaml --steps sumstat,weights,score,finalize -i /path/to/sumstat"
  exit 1
fi

# (legacy) no-op: --sbatch-array removed; keep block absent

config_file_host=$(realpath "$config_file")

# Load config + resolve active methods (CLI --methods overrides config methods:)
parse_config "$config_file_host"
load_active_methods "$methods_cli" "$config_file_host"
echo "Active methods: $(methods_to_csv "$CFG_METHODS")"

################################################################################
# Parse config file (simple YAML parsing with awk)
################################################################################
parse_yaml_value() {
  local key="$1"
  local file="$2"
  awk -F': ' -v key="$key" '$1 == key {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$file"
}

# Parse nested YAML values (supports one or two levels under a section).
# Defined in bin/lib/common.sh as parse_yaml_nested().

# Parse a YAML list (e.g., modules:\n  - tools\n  - singularity/4.1.2)
# Defined in bin/lib/common.sh as parse_yaml_list().

# Parse inline YAML dict (e.g., "{ mem: 10g, cpus: 6, time: '1:00:00' }")
parse_inline_dict() {
  local dict="$1"
  local key="$2"
  echo "$dict" | sed 's/[{}]//g' | tr ',' '\n' | awk -F': ' -v key="$key" '
    $1 ~ key { gsub(/^[ \t]+|[ \t]+$|'"'"'/, "", $2); print $2 }
  '
}

# Format SBATCH settings for logging.
format_sbatch_settings() {
  local array_spec="${1:-none}"
  local max_parallel="${2:-}"
  local account="${slurm_account:-none}"
  local partition="${slurm_partition:-none}"
  local mem="${slurm_mem:-}"
  local cpus="${slurm_cpus:-}"
  local time="${slurm_time:-}"
  if [[ -n "$max_parallel" ]]; then
    echo "SBATCH: account=${account}, partition=${partition}, mem=${mem}, cpus=${cpus}, time=${time}, array=${array_spec}, max_parallel=${max_parallel}"
  else
    echo "SBATCH: account=${account}, partition=${partition}, mem=${mem}, cpus=${cpus}, time=${time}, array=${array_spec}"
  fi
}

# Validate sumstat metadata has required N fields before processing.
# This provides early feedback for missing sample size information.
validate_sumstat_metadata() {
  local sumstat_dir="$1"
  local config_file="$2"
  
  local metadata_file="${sumstat_dir}/cleaned_metadata.yaml"
  if [[ ! -f "$metadata_file" ]]; then
    >&2 echo ""
    >&2 echo "Error: Metadata file not found: ${metadata_file}"
    >&2 echo ""
    >&2 echo "The sumstat directory must contain a cleaned_metadata.yaml file"
    >&2 echo "(produced by cleansumstats)."
    >&2 echo ""
    exit 1
  fi
  
  # Determine which N field is needed based on config
  local which_n
  which_n=$(parse_yaml_value "whichn" "$config_file")
  which_n="${which_n:-totalN}"
  
  local n_value=""
  local n_field=""
  local alt_fields=""
  
  if [[ "$which_n" == "effectiveN" ]]; then
    n_field="stats_EffectiveN"
    n_value=$(awk -F': ' '$1=="stats_EffectiveN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
    alt_fields="stats_CaseN + stats_ControlN"
    
    # If effectiveN is missing, check if case/control counts are available (can derive effectiveN)
    if [[ -z "$n_value" ]]; then
      local case_n ctrl_n
      case_n=$(awk -F': ' '$1=="stats_CaseN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
      ctrl_n=$(awk -F': ' '$1=="stats_ControlN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
      if [[ -n "$case_n" ]] && [[ -n "$ctrl_n" ]] && [[ "$case_n" =~ ^[0-9]+$ ]] && [[ "$ctrl_n" =~ ^[0-9]+$ ]]; then
        # Can derive effectiveN from case/control
        n_value="derivable"
      fi
    fi
  else
    # totalN (default)
    n_field="stats_TotalN"
    n_value=$(awk -F': ' '$1=="stats_TotalN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
    alt_fields="stats_CaseN + stats_ControlN (for effectiveN)"
  fi
  
  if [[ -z "$n_value" ]]; then
    >&2 echo ""
    >&2 echo "Error: Required sample size field is missing or empty in metadata."
    >&2 echo ""
    >&2 echo "  Config:    whichn: ${which_n}"
    >&2 echo "  Expected:  ${n_field} (or ${alt_fields})"
    >&2 echo "  Metadata:  ${metadata_file}"
    >&2 echo ""
    >&2 echo "The posteriors calculation requires sample size (N) to be present."
    >&2 echo "Please ensure the cleansumstats output includes this field, or"
    >&2 echo "manually add '${n_field}: <value>' to the metadata file."
    >&2 echo ""
    exit 1
  fi
  
  echo "Metadata validated: ${n_field} = ${n_value}"
}

# Load modules from config (e.g., singularity on HPC systems that use module load)
while IFS= read -r mod; do
  [[ -z "$mod" ]] && continue
  echo "Loading module: $mod"
  module load "$mod"
done < <(parse_yaml_list "modules" "$config_file_host")

# Read paths from config (support new keys + legacy aliases)
cfg_ld_reference=$(parse_yaml_value "ld_reference" "$config_file_host")
cfg_lddir=$(parse_yaml_value "lddir" "$config_file_host")
if [[ -z "$cfg_ld_reference" && -n "$cfg_lddir" ]]; then
  cfg_ld_reference="$cfg_lddir"
fi

cfg_genotypes=$(parse_yaml_value "genotypes" "$config_file_host")
cfg_genodir=$(parse_yaml_value "genodir" "$config_file_host")
if [[ -z "$cfg_genotypes" && -n "$cfg_genodir" ]]; then
  cfg_genotypes="$cfg_genodir"
fi

cfg_genotype_manifest=$(parse_yaml_value "genotype_manifest" "$config_file_host")
cfg_genofile=$(parse_yaml_value "genofile" "$config_file_host")
if [[ -z "$cfg_genotype_manifest" && -n "$cfg_genofile" ]]; then
  cfg_genotype_manifest="$cfg_genofile"
fi

cfg_outdir=$(parse_yaml_value "outdir" "$config_file_host")
cfg_chromosomes=$(parse_yaml_value "chromosomes" "$config_file_host")

# Read optional reference files (legacy INFO/MAF filtering)
cfg_info_file=$(parse_yaml_nested "references" "info_file" "$config_file_host")
cfg_maf_file=$(parse_yaml_nested "references" "maf_file" "$config_file_host")

# Read optional user inclusion lists
cfg_inclusion_gt=$(parse_yaml_nested "filters" "inclusion_list.gt" "$config_file_host")
cfg_inclusion_ss=$(parse_yaml_nested "filters" "inclusion_list.ss" "$config_file_host")
cfg_inclusion_ld=$(parse_yaml_nested "filters" "inclusion_list.ld" "$config_file_host")

# Read genome build and liftover reference (for dual-position mapfile)
cfg_genotype_build=$(parse_yaml_value "genotype_build" "$config_file_host")
cfg_genotype_build="${cfg_genotype_build:-GRCh37}"  # Default to GRCh37
cfg_liftover_reference=$(parse_yaml_value "liftover_reference" "$config_file_host")

# CLI overrides config
if [[ -n "$outdir" ]]; then
  cfg_outdir="$outdir"
fi

if [[ -n "$chromosomes_override" ]]; then
  cfg_chromosomes="$chromosomes_override"
fi

################################################################################
# Handle --sbatch-array: submit chromosome-parallel steps as a SLURM job array
################################################################################
expand_chromosome_list() {
  local chr_spec="$1"
  if [[ -z "$chr_spec" ]]; then
    chr_spec="1-22"
  fi
  # Normalize separators
  chr_spec=$(echo "$chr_spec" | tr ',' ' ')
  local out=""
  local tok
  for tok in $chr_spec; do
    if [[ "$tok" =~ ^[0-9]+-[0-9]+$ ]]; then
      local start="${tok%-*}"
      local end="${tok#*-}"
      if [[ "$start" -le "$end" ]]; then
        local c
        for ((c=start; c<=end; c++)); do out="${out} ${c}"; done
      else
        local c
        for ((c=start; c>=end; c--)); do out="${out} ${c}"; done
      fi
    else
      out="${out} ${tok}"
    fi
  done
  echo "$out" | awk '{$1=$1;print}'
}

format_elapsed() {
  local secs="$1"
  if [[ -z "$secs" ]] || ! [[ "$secs" =~ ^[0-9]+$ ]]; then
    echo "NA"
    return
  fi
  local h=$((secs/3600))
  local m=$(((secs%3600)/60))
  local s=$((secs%60))
  printf "%02d:%02d:%02d" "$h" "$m" "$s"
}

cleanup_work_and_tmp() {
  local outdir="$1"
  local sumstat="${2:-}"

  # Always safe to remove prep tmp if present.
  rm -rf "${outdir}/prep/tmp" 2>/dev/null || true

  if [[ -n "$sumstat" ]]; then
    rm -rf "${outdir}/sumstats/${sumstat}/tmp" 2>/dev/null || true
    rm -rf "${outdir}/sumstats/${sumstat}/work" 2>/dev/null || true
  fi
}

################################################################################
# Prerequisite checking helpers (must be defined before driver mode runs)
################################################################################
check_prep_exists() {
  local outdir="$1"
  # Optional: validate that prep outputs exist for the chromosomes requested by the current config.
  # This prevents re-using a previously "successful" prep run from a different chromosome subset.
  local chr_spec="${2:-${cfg_chromosomes:-1-22}}"
  local chr_list=""
  chr_list=$(expand_chromosome_list "$chr_spec")

  local missing=""

  # Check for genotype prep outputs
  if [[ ! -d "${outdir}/prep/genotypes" ]] || [[ -z "$(ls -A "${outdir}/prep/genotypes" 2>/dev/null)" ]]; then
    missing="${missing}  - Genotype prep: ${outdir}/prep/genotypes/\n"
  fi

  # Check for inclusion list (v2.1)
  if [[ ! -f "${outdir}/prep/inclusion_list/variant_inclusion_list.tsv" ]]; then
    missing="${missing}  - Inclusion list: ${outdir}/prep/inclusion_list/variant_inclusion_list.tsv\n"
  fi
  if [[ ! -f "${outdir}/prep/inclusion_list/.rsid_index" ]]; then
    missing="${missing}  - Inclusion index: ${outdir}/prep/inclusion_list/.rsid_index\n"
  fi
  if [[ ! -f "${outdir}/prep/variant_map.tsv" ]]; then
    missing="${missing}  - Variant map: ${outdir}/prep/variant_map.tsv\n"
  fi

  # Chromosome-specific prep artifacts (must exist for all requested chromosomes)
  if [[ -n "$chr_list" ]]; then
    local chr=""
    for chr in $chr_list; do
      if [[ ! -f "${outdir}/prep/genotypes/chr${chr}_pvar_fmt" ]]; then
        missing="${missing}  - Genotype pvar (chr${chr}): ${outdir}/prep/genotypes/chr${chr}_pvar_fmt\n"
      fi
      if [[ ! -f "${outdir}/prep/ldref/chr${chr}_ld_rsids" ]]; then
        missing="${missing}  - LD ref rsids (chr${chr}): ${outdir}/prep/ldref/chr${chr}_ld_rsids\n"
      fi
      if [[ ! -f "${outdir}/prep/inclusion_list/chr${chr}_variant_map" ]]; then
        missing="${missing}  - Inclusion variant map (chr${chr}): ${outdir}/prep/inclusion_list/chr${chr}_variant_map\n"
      fi
    done
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
  local methods_list="${3:-${CFG_METHODS:-sbayesr}}"

  local formatted_new="${outdir}/sumstats/${sumstat}/work/formatted/sumstat_formatted.tsv.gz"
  local formatted_legacy="${outdir}/sumstats/${sumstat}/formatted/sumstat_formatted.tsv.gz"
  local formatted_chr_new_glob="${outdir}/sumstats/${sumstat}/work/formatted/chr*.tsv"
  local work_base="${outdir}/sumstats/${sumstat}/work"
  local missing=""

  if [[ ! -f "$formatted_new" && ! -f "$formatted_legacy" ]]; then
    if ! compgen -G "$formatted_chr_new_glob" >/dev/null 2>&1; then
      missing="${missing}  - Formatted sumstat (expected): ${formatted_new}\n"
      missing="${missing}    (alt): ${formatted_chr_new_glob}\n"
    fi
  fi

  local m
  for m in $methods_list; do
    local filtered_dir="${work_base}/filtered_${m}"
    local filtered_chr_glob="${filtered_dir}/chr*_filtered.tsv"
    if compgen -G "$filtered_chr_glob" >/dev/null 2>&1; then
      continue
    fi
    if [[ "$m" == "sbayesr" ]]; then
      local legacy_filtered="${work_base}/filtered/chr*_filtered.tsv"
      local legacy_root="${outdir}/sumstats/${sumstat}/filtered/chr*_filtered.tsv"
      if compgen -G "$legacy_filtered" >/dev/null 2>&1 || compgen -G "$legacy_root" >/dev/null 2>&1; then
        continue
      fi
    fi
    missing="${missing}  - Filtered sumstat (${m}): ${filtered_chr_glob}\n"
  done

  if [[ -n "$missing" ]]; then
    echo -e "$missing"
    return 1
  fi
  return 0
}

check_posteriors_exists() {
  local outdir="$1"
  local sumstat="$2"
  local methods_list="${3:-${CFG_METHODS:-sbayesr}}"
  local work_base="${outdir}/sumstats/${sumstat}/work"
  local missing=""
  local m

  for m in $methods_list; do
    case "$m" in
      sbayesr)
        local post_dir="${work_base}/posteriors"
        local mapped_dir="${work_base}/posteriors_mapped"
        local post_legacy="${outdir}/sumstats/${sumstat}/posteriors"
        local mapped_legacy="${outdir}/sumstats/${sumstat}/posteriors_mapped"
        if [[ ( ! -d "$post_dir" || -z "$(ls -A "$post_dir" 2>/dev/null)" ) && ( ! -d "$post_legacy" || -z "$(ls -A "$post_legacy" 2>/dev/null)" ) ]]; then
          missing="${missing}  - Posteriors (sbayesr): ${post_dir}/\n"
        fi
        if [[ ( ! -d "$mapped_dir" || -z "$(ls -A "$mapped_dir" 2>/dev/null)" ) && ( ! -d "$mapped_legacy" || -z "$(ls -A "$mapped_legacy" 2>/dev/null)" ) ]]; then
          missing="${missing}  - Posteriors mapped (sbayesr): ${mapped_dir}/\n"
        fi
        ;;
      ldpred2)
        local ldp_post="${work_base}/posteriors_ldpred2"
        local ldp_mapped="${work_base}/posteriors_mapped_ldpred2"
        if [[ ! -d "$ldp_post" || -z "$(ls -A "$ldp_post" 2>/dev/null)" ]]; then
          missing="${missing}  - Posteriors (ldpred2): ${ldp_post}/\n"
        fi
        if [[ ! -d "$ldp_mapped" || -z "$(ls -A "$ldp_mapped" 2>/dev/null)" ]]; then
          missing="${missing}  - Posteriors mapped (ldpred2): ${ldp_mapped}/\n"
        fi
        ;;
    esac
  done

  if [[ -n "$missing" ]]; then
    echo -e "$missing"
    return 1
  fi
  return 0
}

# Build filter-variants invocations for all active methods (sumstat array tasks).
build_sumstat_filter_chr_cmds() {
  local run_cmd="$1"
  local body=""
  local m
  for m in ${CFG_METHODS:-sbayesr}; do
    body="${body}${run_cmd} --method ${m} --_chr \"\$CHR\"; "
  done
  echo "$body"
}

# Resolve mapped-posteriors directory for a score step profile.
score_profile_mapped_dir() {
  local outdir_host="$1"
  local sumstat_name="$2"
  local step_profile="$3"
  local method=""
  case "$step_profile" in
    score_sbayesr|weights_sbayesr) method="sbayesr" ;;
    score_ldpred2) method="ldpred2" ;;
    score) method="sbayesr" ;;
    *) echo ""; return 0 ;;
  esac
  local base="${outdir_host}/sumstats/${sumstat_name}/work"
  if [[ "$method" == "sbayesr" ]]; then
    if [[ -d "${base}/posteriors_mapped" ]]; then
      echo "${base}/posteriors_mapped"
    elif [[ -d "${outdir_host}/sumstats/${sumstat_name}/posteriors_mapped" ]]; then
      echo "${outdir_host}/sumstats/${sumstat_name}/posteriors_mapped"
    fi
  else
    echo "${base}/posteriors_mapped_ldpred2"
  fi
}

if [[ "$driver_run" == true ]]; then
  # Support: any combination of:
  #   --steps sumstat
  #   --steps posteriors
  #   --steps score
  #   --steps sumstat,posteriors,score
  # Order is enforced: sumstat -> posteriors -> score
  if [[ -z "${steps_arg:-}" ]]; then
    >&2 echo "Error: --_driver-run requires --steps (prep, sumstat, weights, score, finalize, or combinations thereof)"
    exit 1
  fi
  has_prep=false
  has_sumstat=false
  has_weights=false
  has_score=false
  has_finalize=false
  IFS=',' read -r -a _sbatch_steps <<< "$steps_arg"
  for _s in "${_sbatch_steps[@]}"; do
    _s="$(echo "$_s" | awk '{$1=$1;print}')"
    [[ -z "$_s" ]] && continue
    if [[ "$_s" == "prep" ]]; then
      has_prep=true
    elif [[ "$_s" == "sumstat" ]]; then
      has_sumstat=true
    elif [[ "$_s" == "weights" ]]; then
      has_weights=true
    elif [[ "$_s" == "score" ]]; then
      has_score=true
    elif [[ "$_s" == "finalize" ]]; then
      has_finalize=true
    else
      >&2 echo "Error: driver mode only supports --steps prep, sumstat, weights, score, finalize (or combinations) (got: '${steps_arg}')"
      exit 1
    fi
  done
  if [[ "$has_prep" == true && ( "$has_sumstat" == true || "$has_weights" == true || "$has_score" == true || "$has_finalize" == true ) ]]; then
    >&2 echo "Error: prep must be run on its own (do not include prep with other steps in driver jobs)"
    exit 1
  fi
  if [[ "$has_prep" != true && "$has_sumstat" != true && "$has_weights" != true && "$has_score" != true && "$has_finalize" != true ]]; then
    >&2 echo "Error: --_driver-run requires --steps to include prep and/or sumstat/weights/score/finalize"
    exit 1
  fi

  # Need outdir for logs and chromosome list file
  if [[ -z "$cfg_outdir" ]]; then
    >&2 echo "Error: outdir not found in config file and -o not provided"
    exit 1
  fi
  mkdir -p "${cfg_outdir}"
  outdir_host=$(realpath "${cfg_outdir}")

  # In driver mode we exit before the later "Resolve paths" block.
  # For non-prep steps, compute sumstat_name here so watcher sanity-checks can run.
  if [[ "$has_prep" != true ]]; then
    if [[ -z "${infold:-}" ]]; then
      >&2 echo "Error: --_driver-run requires -i <sumstat_dir> for sumstat/weights/score"
      exit 1
    fi
infold_host=$(realpath "${infold}")
    if [[ ! -d "$infold_host" ]]; then
      >&2 echo "Error: Input directory doesn't exist: $infold_host"
      exit 1
    fi
    sumstat_name=$(basename "$infold_host")
  fi

  # Read SLURM settings from config
  slurm_account=$(parse_yaml_nested "slurm" "account" "$config_file_host")
  slurm_partition=$(parse_yaml_nested "slurm" "partition" "$config_file_host")

  # Build chromosome list file (robust to non-contiguous ranges)
  chr_list=$(expand_chromosome_list "${cfg_chromosomes:-1-22}")
  if [[ -z "$chr_list" ]]; then
    >&2 echo "Error: failed to determine chromosomes list (cfg_chromosomes='${cfg_chromosomes}')"
    exit 1
  fi
  chr_count=$(echo "$chr_list" | wc -w | awk '{print $1}')

  # Keep SLURM logs contained within the appropriate folder.
  # (prep is shared; sumstat/weights/score are per-sumstat)
  if [[ "$has_prep" == true ]]; then
    log_dir="${outdir_host}/prep/logs/slurm"
  else
    log_dir="${outdir_host}/sumstats/${sumstat_name}/logs/slurm"
  fi
  mkdir -p "$log_dir"

  watch_array() {
    local array_jobid="$1"
    local chr_file="$2"
    local chr_count="$3"

    is_terminal_state() {
      local st="$1"
      case "$st" in
        COMPLETED*|FAILED*|CANCELLED*|TIMEOUT*|OUT_OF_MEMORY*|NODE_FAIL*|PREEMPTED*|BOOT_FAIL*|DEADLINE*|SPECIAL_EXIT*)
          return 0
          ;;
        *)
          return 1
          ;;
      esac
    }

    # Best-effort watch: report task starts/finishes + elapsed time.
    # Important: must terminate when the array is done, even if sacct is delayed.
    if ! command -v squeue >/dev/null 2>&1; then
      echo "Note: squeue not available; not watching job progress."
      return 0
    fi

    echo ""
    echo "Watching array progress (Ctrl+C stops watching; jobs continue)..."

    declare -A started
    declare -A finished
    declare -A task_state
    declare -A task_failed
    start_ts=$(date +%s)

    while [[ ${#finished[@]} -lt $chr_count ]]; do
      # Snapshot current tasks in queue.
      # IMPORTANT: squeue can show arrays as a single aggregated line (e.g. 123_[1-22]).
      # Use %A (array job id) and %a (task id) to get per-task rows when available.
      mapfile -t sq_lines < <(squeue -h --array -j "${array_jobid}" -o "%A|%a|%T" 2>/dev/null || squeue -h -j "${array_jobid}" -o "%A|%a|%T" 2>/dev/null || true)

      declare -A in_queue
      for line in "${sq_lines[@]}"; do
        jid="${line%%|*}"
        rest="${line#*|}"
        tid="${rest%%|*}"
        state="${rest##*|}"

        # tid can be numeric (single task) or a range/list (aggregated). We only act on numeric rows.
        if [[ "$jid" == "$array_jobid" ]] && [[ "$tid" =~ ^[0-9]+$ ]]; then
          idx="$tid"
          in_queue["$idx"]="$state"
          if [[ -z "${started[$idx]:-}" ]] && [[ "$state" == "RUNNING" || "$state" == "COMPLETING" ]]; then
            started["$idx"]=1
            chr=$(sed -n "${idx}p" "$chr_file" 2>/dev/null || true)
            echo "  Started: task=${idx} chr=${chr} state=${state} time=$(date)"
          fi
        fi
      done

      # Detect finished tasks (not in queue anymore)
      for ((idx=1; idx<=chr_count; idx++)); do
        if [[ -n "${finished[$idx]:-}" ]]; then
          continue
        fi
        if [[ -z "${in_queue[$idx]:-}" ]]; then
          # Try to fetch accounting info
          chr=$(sed -n "${idx}p" "$chr_file" 2>/dev/null || true)
          if command -v sacct >/dev/null 2>&1; then
            acct_line=$(sacct -j "${array_jobid}_${idx}" --format=JobIDRaw,State,ElapsedRaw -n -P 2>/dev/null | head -n 1 || true)
            if [[ -n "$acct_line" ]]; then
              rest="${acct_line#*|}"
              st="${rest%%|*}"
              elapsed_raw="${rest##*|}"
              elapsed_fmt=$(format_elapsed "$elapsed_raw")
              task_state["$idx"]="$st"
              if is_terminal_state "$st"; then
                finished["$idx"]=1
                if [[ "$st" != COMPLETED* ]]; then
                  task_failed["$idx"]=1
                fi
                echo "  Finished: task=${idx} chr=${chr} state=${st} elapsed=${elapsed_fmt} time=$(date)"
              fi
            fi
          else
            # No sacct; mark finished when it leaves queue
            finished["$idx"]=1
            task_state["$idx"]="UNKNOWN"
            echo "  Finished: task=${idx} chr=${chr} (no sacct available) time=$(date)"
          fi
        fi
      done

      # If the entire array disappeared from squeue, finish promptly.
      # This prevents hangs if sacct is delayed for some tasks.
      if [[ ${#sq_lines[@]} -eq 0 ]]; then
        sleep 2
        for ((idx=1; idx<=chr_count; idx++)); do
          if [[ -n "${finished[$idx]:-}" ]]; then
            continue
          fi
          chr=$(sed -n "${idx}p" "$chr_file" 2>/dev/null || true)
          if command -v sacct >/dev/null 2>&1; then
            acct_line=$(sacct -j "${array_jobid}_${idx}" --format=State,ElapsedRaw -n -P 2>/dev/null | head -n 1 || true)
            if [[ -n "$acct_line" ]]; then
              st="${acct_line%%|*}"
              elapsed_raw="${acct_line##*|}"
              elapsed_fmt=$(format_elapsed "$elapsed_raw")
              task_state["$idx"]="$st"
              if is_terminal_state "$st"; then
                finished["$idx"]=1
                if [[ "$st" != COMPLETED* ]]; then
                  task_failed["$idx"]=1
                fi
                echo "  Finished: task=${idx} chr=${chr} state=${st} elapsed=${elapsed_fmt} time=$(date)"
              fi
            fi
          else
            finished["$idx"]=1
            task_state["$idx"]="UNKNOWN"
            echo "  Finished: task=${idx} chr=${chr} (array ended) time=$(date)"
          fi
        done
      fi

      sleep 15
    done

    # sacct can lag; do a final poll to resolve UNKNOWN tasks and capture failures.
    if command -v sacct >/dev/null 2>&1; then
      local tries=0
      while [[ $tries -lt 12 ]]; do
        local unknown_left=0
        for ((idx=1; idx<=chr_count; idx++)); do
          if [[ "${task_state[$idx]:-}" != "UNKNOWN" ]]; then
            continue
          fi
          acct_line=$(sacct -j "${array_jobid}_${idx}" --format=State -n -P 2>/dev/null | head -n 1 || true)
          if [[ -n "$acct_line" ]]; then
            st="${acct_line%%|*}"
            task_state["$idx"]="$st"
            if [[ "$st" != COMPLETED* ]]; then
              task_failed["$idx"]=1
            fi
          else
            unknown_left=$((unknown_left + 1))
          fi
        done
        if [[ $unknown_left -eq 0 ]]; then
          break
        fi
        tries=$((tries + 1))
        sleep 5
      done
    fi

    end_ts=$(date +%s)
    total_elapsed=$((end_ts - start_ts))
    echo ""
    echo "Array completed: job=${array_jobid} total_elapsed=$(format_elapsed "$total_elapsed")"
    local failed_count=0
    for ((idx=1; idx<=chr_count; idx++)); do
      if [[ -n "${task_failed[$idx]:-}" ]]; then
        failed_count=$((failed_count + 1))
      fi
    done
    if [[ $failed_count -gt 0 ]]; then
      echo "Array result: FAILED_TASKS=${failed_count}/${chr_count}"
      return 1
    fi

    # If some tasks are still UNKNOWN (no sacct), fall back to parent array job state.
    local unknown_count=0
    for ((idx=1; idx<=chr_count; idx++)); do
      if [[ "${task_state[$idx]:-}" == "UNKNOWN" ]]; then
        unknown_count=$((unknown_count + 1))
      fi
    done
    if [[ $unknown_count -gt 0 ]]; then
      parent_state=$(sacct -j "${array_jobid}" --format=State -n -P 2>/dev/null | head -n 1 || true)
      parent_state="${parent_state%%|*}"
      echo "Array result: UNKNOWN_TASKS=${unknown_count}/${chr_count} (parent_state=${parent_state:-UNKNOWN})"
      if [[ "$parent_state" == COMPLETED* ]]; then
        >&2 echo "Warning: some task states are UNKNOWN due to accounting lag; treating array as successful because parent job is COMPLETED."
        return 0
      fi
      return 1
    fi
    return 0
  }

  watch_job() {
    local jobid="$1"
    local label="${2:-job}"

    if ! command -v squeue >/dev/null 2>&1; then
      echo "Note: squeue not available; not watching ${label}."
      return 0
    fi

    echo ""
    echo "Watching ${label} (Ctrl+C stops watching; job continues)..."
    local started=0
    local start_ts
    start_ts=$(date +%s)

    while true; do
      st=$(squeue -h -j "${jobid}" -o "%T" 2>/dev/null | head -n 1 || true)
      if [[ -n "$st" ]]; then
        if [[ $started -eq 0 ]]; then
          started=1
          echo "  Started: job=${jobid} state=${st} time=$(date)"
        fi
        sleep 10
        continue
      fi
      break
    done

    # Final accounting (best effort)
    local final_state="UNKNOWN"
    local elapsed_fmt="NA"
    if command -v sacct >/dev/null 2>&1; then
      # sacct can lag; try a few times
      local tries=0
      while [[ $tries -lt 12 ]]; do
        acct_line=$(sacct -j "${jobid}" --format=State,ElapsedRaw -n -P 2>/dev/null | head -n 1 || true)
        if [[ -n "$acct_line" ]]; then
          final_state="${acct_line%%|*}"
          elapsed_raw="${acct_line##*|}"
          elapsed_fmt=$(format_elapsed "$elapsed_raw")
          break
        fi
        tries=$((tries + 1))
        sleep 5
      done
    fi

    end_ts=$(date +%s)
    total_elapsed=$((end_ts - start_ts))
    echo "  Finished: job=${jobid} state=${final_state} elapsed=${elapsed_fmt} total_elapsed=$(format_elapsed "$total_elapsed") time=$(date)"

    if [[ "$final_state" == COMPLETED* ]]; then
      return 0
    fi
    return 1
  }

  submit_array_for_step() {
    local step_profile="$1"

    step_settings=$(resolve_slurm_step_settings "$step_profile" "$config_file_host")
    slurm_mem=""
    slurm_cpus=""
    slurm_time=""
    slurm_max_parallel_step=""
    if [[ -n "$step_settings" ]]; then
      slurm_mem=$(parse_inline_dict "$step_settings" "mem")
      slurm_cpus=$(parse_inline_dict "$step_settings" "cpus")
      slurm_time=$(parse_inline_dict "$step_settings" "time")
      slurm_max_parallel_step=$(parse_inline_dict "$step_settings" "max_parallel")
    fi
    # Step-specific defaults
    if [[ "$step_profile" == "sumstat" ]]; then
      slurm_mem="${slurm_mem:-1g}"
      slurm_cpus="${slurm_cpus:-1}"
      slurm_time="${slurm_time:-0:30:00}"
    elif [[ "$step_profile" == "prep" ]]; then
      slurm_mem="${slurm_mem:-10g}"
      slurm_cpus="${slurm_cpus:-2}"
      slurm_time="${slurm_time:-2:00:00}"
    elif [[ "$step_profile" == "weights_sbayesr" ]]; then
      slurm_mem="${slurm_mem:-20g}"
      slurm_cpus="${slurm_cpus:-8}"
      slurm_time="${slurm_time:-2:00:00}"
    elif [[ "$step_profile" == "weights_benchmark" ]]; then
      slurm_mem="${slurm_mem:-2g}"
      slurm_cpus="${slurm_cpus:-2}"
      slurm_time="${slurm_time:-0:30:00}"
    elif [[ "$step_profile" == "weights_ldpred2" ]]; then
      slurm_mem="${slurm_mem:-64g}"
      slurm_cpus="${slurm_cpus:-16}"
      slurm_time="${slurm_time:-4:00:00}"
    elif [[ "$step_profile" == "score_sbayesr" || "$step_profile" == "score_ldpred2" || "$step_profile" == "score" ]]; then
      slurm_mem="${slurm_mem:-10g}"
      slurm_cpus="${slurm_cpus:-4}"
      slurm_time="${slurm_time:-0:30:00}"
    else
      slurm_mem="${slurm_mem:-20g}"
      slurm_cpus="${slurm_cpus:-8}"
      slurm_time="${slurm_time:-2:00:00}"
    fi

    # Determine max parallel tasks (step > default 22)
    max_parallel="${slurm_max_parallel_step:-22}"
    if ! [[ "$max_parallel" =~ ^[0-9]+$ ]] || [[ "$max_parallel" -lt 1 ]]; then
      >&2 echo "Error: invalid max_parallel for ${step_profile}: $max_parallel"
      exit 1
    fi

    # Build job name
    if [[ -n "$infold" ]]; then
      job_name="pgs_$(basename "$infold")_${step_profile}"
    else
      job_name="pgs_${step_profile}"
    fi

    # Determine chromosomes for this array.
    # Default: config chromosomes (chr_list). For scoring, only run chromosomes that produced
    # non-empty mapped posteriors (so we don't waste time scoring empty chromosomes).
    step_chr_list="$chr_list"
    step_chr_count="$chr_count"
    if [[ "$step_profile" == score* && -n "${sumstat_name:-}" ]]; then
      mapped_dir=$(score_profile_mapped_dir "$outdir_host" "$sumstat_name" "$step_profile")

      if [[ -n "$mapped_dir" && -d "$mapped_dir" ]]; then
        map_chrs=()
        for f in "${mapped_dir}"/chr*.snpRes; do
          [[ -f "$f" ]] || continue
          nlines=$(wc -l < "$f" 2>/dev/null || echo 0)
          if [[ "$nlines" -gt 1 ]]; then
            chr="${f##*/}"
            chr="${chr#chr}"
            chr="${chr%.snpRes}"
            map_chrs+=("$chr")
          fi
        done
        if [[ ${#map_chrs[@]} -gt 0 ]]; then
          step_chr_list=$(printf '%s\n' "${map_chrs[@]}" | sort -n | tr '\n' ' ' | awk '{$1=$1;print}')
          step_chr_count=$(echo "$step_chr_list" | wc -w | awk '{print $1}')
        else
          >&2 echo "Warning: no non-empty mapped posteriors found; skipping score array."
          return 0
        fi
      fi
    fi

    chr_file="${log_dir}/${job_name}.chromosomes.txt"
    echo "$step_chr_list" | tr ' ' '\n' > "$chr_file"

    # In array mode, we must only run chromosome-parallel work.
    # - sumstat: filter-variants per-chromosome (format-sumstat already ran in driver)
    # - weights_sbayesr: calc-posteriors + format-posteriors per chr
    # - weights_benchmark: calc-benchmark per chr
    # - score: only calc-score is chr-parallel; combine-scores/finalize-output must be run once after
    steps_arg_for_task="${step_profile}"
    if [[ "$step_profile" == "sumstat" ]]; then
      steps_arg_for_task="filter-variants"
    elif [[ "$step_profile" == "score" || "$step_profile" == "score_sbayesr" || "$step_profile" == "score_ldpred2" ]]; then
      steps_arg_for_task="calc-score"
    elif [[ "$step_profile" == "prep" ]]; then
      steps_arg_for_task="prep-inclusion-list"
    elif [[ "$step_profile" == "weights_sbayesr" ]]; then
      steps_arg_for_task="calc-posteriors,format-posteriors"
    elif [[ "$step_profile" == "weights_benchmark" ]]; then
      steps_arg_for_task="calc-benchmark"
    fi

    run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps ${steps_arg_for_task}"
    [[ -n "$infold" ]] && run_cmd="${run_cmd} -i ${infold}"
    [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"
    [[ -n "$force_mode" ]] && run_cmd="${run_cmd} --force"
    run_cmd="${run_cmd} --methods $(methods_to_csv "$CFG_METHODS")"

    local method_flag=""
    case "$step_profile" in
      weights_sbayesr|score_sbayesr|score)
        method_flag="--method sbayesr"
        ;;
      score_ldpred2)
        method_flag="--method ldpred2"
        ;;
    esac
    [[ -n "$method_flag" ]] && run_cmd="${run_cmd} ${method_flag}"

    local chr_body=""
    if [[ "$step_profile" == "sumstat" ]]; then
      chr_body=$(build_sumstat_filter_chr_cmds "$run_cmd")
    else
      chr_body="${run_cmd} --_chr \"\$CHR\"; "
    fi

    task_wrap="CHR=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" \"${chr_file}\"); \
if [[ -z \"\$CHR\" ]]; then echo \"Error: could not resolve chromosome for task \$SLURM_ARRAY_TASK_ID\" >&2; exit 1; fi; \
echo \"[INFO] Starting ${step_profile} chr\${CHR} at \$(date)\"; \
${chr_body}\
rc=\$?; echo \"[INFO] Finished ${step_profile} chr\${CHR} at \$(date) (exit=\$rc)\"; exit \$rc"

    # Build sbatch args as an array to avoid brittle quoting + eval issues.
    # NOTE: task_wrap intentionally contains escaped '$' so it is evaluated on the compute node, not here.
    sbatch_args=(--parsable)
    sbatch_args+=(--mem="${slurm_mem}")
    sbatch_args+=(--cpus-per-task="${slurm_cpus}")
    sbatch_args+=(--time="${slurm_time}")
    sbatch_args+=(--job-name="${job_name}")
    sbatch_args+=(--output="${log_dir}/${job_name}_%A_%a.out")
    sbatch_args+=(--error="${log_dir}/${job_name}_%A_%a.err")
    sbatch_args+=(--array="1-${step_chr_count}%${max_parallel}")
    [[ -n "$slurm_account" ]] && sbatch_args+=(--account="${slurm_account}")
    [[ -n "$slurm_partition" ]] && sbatch_args+=(--partition="${slurm_partition}")
    sbatch_args+=(--wrap="${task_wrap}")

    echo "Submitting SLURM job array..."
    echo "  Job name: ${job_name}"
    echo "  Step: ${step_profile}"
    if [[ "$step_profile" == score* ]]; then
      echo "  Note: array runs 'calc-score' only; run finalize step (separate job) for combine-scores + finalize-output."
    fi
    echo "  $(format_sbatch_settings "1-${step_chr_count}%${max_parallel}" "${max_parallel}")"
    echo "  Chromosomes: ${step_chr_list}"
    echo "  Array: 1-${step_chr_count}%${max_parallel} (max_parallel=${max_parallel})"
    echo "  Resources per task: mem=${slurm_mem}, cpus=${slurm_cpus}, time=${slurm_time}"
    echo "  Logs: ${log_dir}/${job_name}_<jobid>_<taskid>.out/.err"
    echo "  Chromosome file: ${chr_file}"
    echo ""

    array_jobid=$(sbatch "${sbatch_args[@]}")
    if [[ -z "$array_jobid" ]]; then
      >&2 echo "Error: failed to submit SLURM array job"
      exit 1
    fi
    echo "Submitted: ${array_jobid}"

    if ! watch_array "$array_jobid" "$chr_file" "$step_chr_count"; then
      >&2 echo "Warning: SLURM array job ${array_jobid} for step '${step_profile}' had failed/unknown task(s)."
      >&2 echo "Continuing (placeholders will propagate); check logs under: ${log_dir}/"
    fi

    # Sanity-check expected outputs exist after a successful array.
    # This catches cases where tasks return 0 but accidentally write nothing.
    if [[ -n "$sumstat_name" ]]; then
      base_sumstat_out="${outdir_host}/sumstats/${sumstat_name}/work"
      if [[ "$step_profile" == "sumstat" ]]; then
        local m n_filtered sumstat_warn=0
        for m in ${CFG_METHODS:-sbayesr}; do
          n_filtered=$(ls "${base_sumstat_out}/filtered_${m}"/chr*_filtered.tsv 2>/dev/null | wc -l | awk '{print $1}')
          if [[ "$m" == "sbayesr" && "$n_filtered" -lt "$step_chr_count" ]]; then
            n_filtered=$(ls "${base_sumstat_out}/filtered"/chr*_filtered.tsv 2>/dev/null | wc -l | awk '{print $1}')
          fi
          if [[ "$n_filtered" -lt "$step_chr_count" ]]; then
            sumstat_warn=1
            >&2 echo "  filtered_${m}: expected >=${step_chr_count} chr*_filtered.tsv (found ${n_filtered})"
          fi
        done
        if [[ "$sumstat_warn" -eq 1 ]]; then
          >&2 echo "Warning: sumstat array finished but outputs are missing (continuing)."
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      elif [[ "$step_profile" == "weights_sbayesr" ]]; then
        n_post=$(ls "${base_sumstat_out}/posteriors"/chr*.snpRes 2>/dev/null | wc -l | awk '{print $1}')
        n_mapped=$(ls "${base_sumstat_out}/posteriors_mapped"/chr*.snpRes 2>/dev/null | wc -l | awk '{print $1}')
        if [[ "$n_post" -lt "$step_chr_count" || "$n_mapped" -lt "$step_chr_count" ]]; then
          >&2 echo "Warning: weights_sbayesr array finished but outputs are missing (continuing)."
          >&2 echo "  Expected >=${step_chr_count} files in:"
          >&2 echo "    - ${base_sumstat_out}/posteriors/chr*.snpRes   (found ${n_post})"
          >&2 echo "    - ${base_sumstat_out}/posteriors_mapped/chr*.snpRes (found ${n_mapped})"
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      elif [[ "$step_profile" == score* ]]; then
        local scores_subdir="scores"
        [[ "$step_profile" == "score_ldpred2" ]] && scores_subdir="scores_ldpred2"
        n_scores=$(ls "${base_sumstat_out}/${scores_subdir}"/chr*.sscore 2>/dev/null | wc -l | awk '{print $1}')
        if [[ "$n_scores" -lt "$step_chr_count" ]]; then
          >&2 echo "Warning: ${step_profile} array finished but outputs are missing (continuing)."
          >&2 echo "  Expected >=${step_chr_count} files in: ${base_sumstat_out}/${scores_subdir}/chr*.sscore (found ${n_scores})"
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      fi
    fi
    echo ""
  }

  submit_single_job_for_step() {
    local step_profile="$1"

    step_settings=$(resolve_slurm_step_settings "$step_profile" "$config_file_host")
    slurm_mem=""
    slurm_cpus=""
    slurm_time=""
    if [[ -n "$step_settings" ]]; then
      slurm_mem=$(parse_inline_dict "$step_settings" "mem")
      slurm_cpus=$(parse_inline_dict "$step_settings" "cpus")
      slurm_time=$(parse_inline_dict "$step_settings" "time")
    fi
    if [[ "$step_profile" == "weights_ldpred2" ]]; then
      slurm_mem="${slurm_mem:-64g}"
      slurm_cpus="${slurm_cpus:-16}"
      slurm_time="${slurm_time:-4:00:00}"
    else
      slurm_mem="${slurm_mem:-20g}"
      slurm_cpus="${slurm_cpus:-8}"
      slurm_time="${slurm_time:-2:00:00}"
    fi

    if [[ -n "$infold" ]]; then
      job_name="pgs_$(basename "$infold")_${step_profile}"
    else
      job_name="pgs_${step_profile}"
    fi

    local steps_arg_for_task="calc-ldpred2,format-posteriors"
    run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps ${steps_arg_for_task}"
    [[ -n "$infold" ]] && run_cmd="${run_cmd} -i ${infold}"
    [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"
    [[ -n "$force_mode" ]] && run_cmd="${run_cmd} --force"
    run_cmd="${run_cmd} --methods $(methods_to_csv "$CFG_METHODS") --method ldpred2"

    sbatch_args=(--parsable)
    sbatch_args+=(--mem="${slurm_mem}")
    sbatch_args+=(--cpus-per-task="${slurm_cpus}")
    sbatch_args+=(--time="${slurm_time}")
    sbatch_args+=(--job-name="${job_name}")
    sbatch_args+=(--output="${log_dir}/${job_name}_%j.out")
    sbatch_args+=(--error="${log_dir}/${job_name}_%j.err")
    [[ -n "$slurm_account" ]] && sbatch_args+=(--account="${slurm_account}")
    [[ -n "$slurm_partition" ]] && sbatch_args+=(--partition="${slurm_partition}")
    sbatch_args+=(--wrap="${run_cmd}")

    echo "Submitting SLURM single job (${step_profile})..."
    echo "  Job name: ${job_name}"
    echo "  Steps: ${steps_arg_for_task} (--method ldpred2)"
    echo "  Resources: mem=${slurm_mem}, cpus=${slurm_cpus}, time=${slurm_time}"
    echo "  Logs: ${log_dir}/${job_name}_<jobid>.out/.err"
    echo ""

    local single_jobid
    single_jobid=$(sbatch "${sbatch_args[@]}")
    if [[ -z "$single_jobid" ]]; then
      >&2 echo "Error: failed to submit SLURM job for ${step_profile}"
      exit 1
    fi
    echo "Submitted: ${single_jobid}"

    if command -v squeue >/dev/null 2>&1; then
      echo "Waiting for ${step_profile} job to finish..."
      while squeue -j "$single_jobid" -h 2>/dev/null | grep -q .; do
        sleep 10
      done
      if command -v sacct >/dev/null 2>&1; then
        st=$(sacct -j "$single_jobid" --format=State -n -P 2>/dev/null | head -n 1 || true)
        if [[ -n "$st" && "$st" != COMPLETED* ]]; then
          >&2 echo "Warning: ${step_profile} job ${single_jobid} finished with state: ${st}"
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      fi
    fi

    if [[ -n "$sumstat_name" && "$step_profile" == "weights_ldpred2" ]]; then
      base_sumstat_out="${outdir_host}/sumstats/${sumstat_name}/work"
      n_post=$(ls "${base_sumstat_out}/posteriors_ldpred2"/chr*.snpRes 2>/dev/null | wc -l | awk '{print $1}')
      n_mapped=$(ls "${base_sumstat_out}/posteriors_mapped_ldpred2"/chr*.snpRes 2>/dev/null | wc -l | awk '{print $1}')
      if [[ "$n_post" -lt 1 || "$n_mapped" -lt 1 ]]; then
        >&2 echo "Warning: weights_ldpred2 finished but posteriors outputs look sparse (continuing)."
        >&2 echo "  posteriors_ldpred2 chr*.snpRes: ${n_post}"
        >&2 echo "  posteriors_mapped_ldpred2 chr*.snpRes: ${n_mapped}"
      fi
    fi
    echo ""
  }

  submit_finalize_job() {
    base_sumstat_out="${outdir_host}/sumstats/${sumstat_name}/work"
    n_scores=0
    for scores_dir in "${base_sumstat_out}/scores" "${base_sumstat_out}/scores_ldpred2"; do
      [[ -d "$scores_dir" ]] || continue
      n=$(ls "${scores_dir}"/chr*.sscore 2>/dev/null | wc -l | awk '{print $1}')
      n_scores=$((n_scores + n))
    done
    if [[ "$n_scores" -eq 0 ]]; then
      >&2 echo "Error: cannot run finalize: no score files under ${base_sumstat_out}/scores*"
      >&2 echo "Run the score step first (e.g. --steps score,finalize or run score then --steps finalize)."
      exit 1
    fi

    step_settings=$(parse_yaml_nested "slurm" "finalize" "$config_file_host")
    slurm_mem=""
    slurm_cpus=""
    slurm_time=""
    if [[ -n "$step_settings" ]]; then
      slurm_mem=$(parse_inline_dict "$step_settings" "mem")
      slurm_cpus=$(parse_inline_dict "$step_settings" "cpus")
      slurm_time=$(parse_inline_dict "$step_settings" "time")
    fi
    slurm_mem="${slurm_mem:-16g}"
    slurm_cpus="${slurm_cpus:-1}"
    slurm_time="${slurm_time:-1:00:00}"

    if [[ -n "$infold" ]]; then
      job_name="pgs_$(basename "$infold")_finalize"
    else
      job_name="pgs_finalize"
    fi
    log_dir="${outdir_host}/sumstats/${sumstat_name}/logs/slurm"
    mkdir -p "$log_dir"

    run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps combine-scores,finalize-output"
    [[ -n "$infold" ]] && run_cmd="${run_cmd} -i ${infold}"
    [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"
    [[ -n "$force_mode" ]] && run_cmd="${run_cmd} --force"
    run_cmd="${run_cmd} --methods $(methods_to_csv "$CFG_METHODS")"

    sbatch_args=(--parsable)
    sbatch_args+=(--mem="${slurm_mem}")
    sbatch_args+=(--cpus-per-task="${slurm_cpus}")
    sbatch_args+=(--time="${slurm_time}")
    sbatch_args+=(--job-name="${job_name}")
    sbatch_args+=(--output="${log_dir}/${job_name}_%j.out")
    sbatch_args+=(--error="${log_dir}/${job_name}_%j.err")
    [[ -n "$slurm_account" ]] && sbatch_args+=(--account="${slurm_account}")
    [[ -n "$slurm_partition" ]] && sbatch_args+=(--partition="${slurm_partition}")
    sbatch_args+=(--wrap="${run_cmd}")

    echo "Submitting SLURM finalize job (combine-scores + finalize-output)..."
    echo "  Job name: ${job_name}"
    echo "  Resources: mem=${slurm_mem}, cpus=${slurm_cpus}, time=${slurm_time}"
    echo "  Logs: ${log_dir}/${job_name}_<jobid>.out/.err"
    echo ""

    finalize_jobid=$(sbatch "${sbatch_args[@]}")
    if [[ -z "$finalize_jobid" ]]; then
      >&2 echo "Error: failed to submit SLURM finalize job"
      exit 1
    fi
    echo "Submitted finalize job: ${finalize_jobid}"

    if command -v squeue >/dev/null 2>&1; then
      echo "Waiting for finalize job to finish..."
      while squeue -j "$finalize_jobid" -h 2>/dev/null | grep -q .; do
        sleep 10
      done
      if command -v sacct >/dev/null 2>&1; then
        st=$(sacct -j "$finalize_jobid" --format=State -n -P 2>/dev/null | head -n 1 || true)
        if [[ -n "$st" && "$st" != COMPLETED* ]]; then
          >&2 echo "Warning: finalize job ${finalize_jobid} finished with state: ${st}"
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      fi
    fi
    echo ""
  }

  # Always run weights before score if both requested
  # In driver mode:
  # - format-sumstat runs directly (GRCh37 paste + chromosome split)
  # - filter-variants runs as a SLURM array (per-chromosome)
  # - weights (sBayesR + benchmark) and score run as SLURM arrays (per-chromosome)
  run_base="${project_dir}/pgscalculator-v2.sh --config ${config_file_host}"
  [[ -n "$infold" ]] && run_base="${run_base} -i ${infold}"
  [[ -n "$outdir" ]] && run_base="${run_base} -o ${outdir}"
  [[ -n "$devmode" ]] && run_base="${run_base} -d"
  [[ -n "$force_mode" ]] && run_base="${run_base} --force"
  run_base="${run_base} --methods $(methods_to_csv "$CFG_METHODS")"

  if [[ "$has_prep" == true ]]; then
    echo "Running prep-genotypes and prep-ldref inside driver job..."
    eval "${run_base} --steps prep-genotypes,prep-ldref"

    # Run prep-inclusion-list as a SLURM array (per-chromosome)
    submit_array_for_step "prep"

    # Combine per-chromosome maps + inclusion list after array completes
    echo "Running prep-inclusion-list combine..."
    eval "${run_base} --steps prep-inclusion-combine"

    if [[ "$do_cleanup" == true ]]; then
      echo "Cleanup enabled: removing prep tmp/"
      cleanup_work_and_tmp "$outdir_host" ""
    fi

    exit 0
  fi

  if [[ "$has_sumstat" == true ]]; then
    # Step 1: Run format-sumstat directly (produces per-chr files)
    echo "Running format-sumstat inside driver job..."
    eval "${run_base} --steps format-sumstat"

    # Check that format-sumstat produced per-chromosome files
    format_dir="${outdir_host}/sumstats/${sumstat_name}/work/formatted"
    if [[ ! -f "${format_dir}/chr1.tsv" ]]; then
      >&2 echo "Error: format-sumstat did not produce per-chromosome files."
      >&2 echo "Expected: ${format_dir}/chr*.tsv"
      exit 1
    fi

    # Step 2: Submit filter-variants as a SLURM array (per-chromosome)
    submit_array_for_step "sumstat"

    # Gate downstream arrays on the expected sumstat outputs existing.
    # This prevents submitting 22 tasks that all fail immediately due to missing inputs.
    if ! check_sumstat_exists "$outdir_host" "$sumstat_name" >/dev/null; then
      >&2 echo "Error: sumstat step finished but expected outputs are missing."
      >&2 echo "Missing:"
      check_sumstat_exists "$outdir_host" "$sumstat_name" || true
      exit 1
    fi
  fi
  if [[ "$has_weights" == true ]]; then
    if has_method sbayesr "$CFG_METHODS"; then
      submit_array_for_step "weights_sbayesr"
    fi
    if has_method ldpred2 "$CFG_METHODS"; then
      submit_single_job_for_step "weights_ldpred2"
    fi
    submit_array_for_step "weights_benchmark"
  fi
  if [[ "$has_score" == true ]]; then
    if has_method sbayesr "$CFG_METHODS"; then
      submit_array_for_step "score_sbayesr"
    fi
    if has_method ldpred2 "$CFG_METHODS"; then
      submit_array_for_step "score_ldpred2"
    fi
  fi
  if [[ "$has_finalize" == true ]]; then
    submit_finalize_job
  fi

  # Optional cleanup (default is to keep work/tmp during development)
  if [[ "$do_cleanup" == true ]]; then
    echo "Cleanup enabled: removing work/ and tmp/ for ${sumstat_name}"
    cleanup_work_and_tmp "$outdir_host" "$sumstat_name"
  fi

  exit 0
fi

################################################################################
# Handle --sbatch: submit as SLURM job
################################################################################
if [[ "$use_sbatch" == true ]]; then
  # Read SLURM settings from config
  slurm_account=$(parse_yaml_nested "slurm" "account" "$config_file_host")
  slurm_partition=$(parse_yaml_nested "slurm" "partition" "$config_file_host")

  # Enforce: prep must be run on its own
  if [[ -z "${steps_arg:-}" ]]; then
    >&2 echo "Error: --sbatch requires --steps"
    exit 1
  fi
  IFS=',' read -r -a _steps_list <<< "$steps_arg"
  has_prep=false
  has_nonprep=false
  for _s in "${_steps_list[@]}"; do
    _s="$(echo "$_s" | awk '{$1=$1;print}')"
    [[ -z "$_s" ]] && continue
    if [[ "$_s" == "prep" ]]; then
      has_prep=true
    else
      has_nonprep=true
    fi
  done
  if [[ "$has_prep" == true && "$has_nonprep" == true ]]; then
    >&2 echo "Error: prep must be run on its own. Run:"
    >&2 echo "  --steps prep --sbatch"
    >&2 echo "and then separately:"
    >&2 echo "  --steps sumstat,weights,score,finalize --sbatch -i <sumstat_dir>"
    exit 1
  fi

  # Determine outdir + log directory
  if [[ -z "${cfg_outdir:-}" ]]; then
    >&2 echo "Error: outdir not found in config file and -o not provided"
    exit 1
  fi
  mkdir -p "${cfg_outdir}"
  outdir_host=$(realpath "${cfg_outdir}")

  if [[ "$has_prep" == true ]]; then
    # Prep driver job (runs prep-genotypes/ldref, submits prep-inclusion-list array)
    step_profile="prep"
    step_settings=$(resolve_slurm_step_settings "$step_profile" "$config_file_host")
    slurm_mem=""
    slurm_cpus=""
    slurm_time=""
    if [[ -n "$step_settings" ]]; then
      slurm_mem=$(parse_inline_dict "$step_settings" "mem")
      slurm_cpus=$(parse_inline_dict "$step_settings" "cpus")
      slurm_time=$(parse_inline_dict "$step_settings" "time")
    fi
    slurm_mem="${slurm_mem:-10g}"
    slurm_cpus="${slurm_cpus:-2}"
    slurm_time="${slurm_time:-2:00:00}"

    job_name="pgs_prep_driver"
    log_dir="${outdir_host}/prep/logs/slurm"
    mkdir -p "$log_dir"

    run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps prep --_driver-run"
    [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"
    [[ -n "$force_mode" ]] && run_cmd="${run_cmd} --force"
    [[ "$do_cleanup" == true ]] && run_cmd="${run_cmd} --cleanup"
    run_cmd="${run_cmd} --methods $(methods_to_csv "$CFG_METHODS")"

    sbatch_args=(--parsable)
    sbatch_args+=(--mem="${slurm_mem}")
    sbatch_args+=(--cpus-per-task="${slurm_cpus}")
    sbatch_args+=(--time="${slurm_time}")
    sbatch_args+=(--job-name="${job_name}")
    sbatch_args+=(--output="${log_dir}/${job_name}_%j.out")
    sbatch_args+=(--error="${log_dir}/${job_name}_%j.err")
    [[ -n "$slurm_account" ]] && sbatch_args+=(--account="${slurm_account}")
    [[ -n "$slurm_partition" ]] && sbatch_args+=(--partition="${slurm_partition}")
    sbatch_args+=(--wrap="${run_cmd}")

    echo "Submitting SLURM job..."
    echo "  Job name: ${job_name}"
    echo "  $(format_sbatch_settings "none")"
    echo "  Resources: mem=${slurm_mem}, cpus=${slurm_cpus}, time=${slurm_time}"
    echo "  Command: ${run_cmd}"
    echo ""

    jobid=$(sbatch "${sbatch_args[@]}")
    echo "Submitted: ${jobid}"
    echo "Driver log: ${log_dir}/${job_name}_${jobid}.out"
    echo "Driver err: ${log_dir}/${job_name}_${jobid}.err"
    echo "Watch: tail -f ${log_dir}/${job_name}_${jobid}.out"
    exit 0
  fi

  # Per-sumstat driver job
  if [[ -z "${infold:-}" ]]; then
    >&2 echo "Error: -i (sumstat input) is required for per-sumstat steps when using --sbatch"
    exit 1
  fi
  infold_host=$(realpath "${infold}")
  if [[ ! -d "$infold_host" ]]; then
    >&2 echo "Error: Input directory doesn't exist: $infold_host"
    exit 1
  fi
  sumstat_name=$(basename "$infold_host")

  # Early validation: check that metadata has required N fields before submitting job
  # This provides fast feedback for a common misconfiguration issue
  validate_sumstat_metadata "$infold_host" "$config_file_host"

  # Driver resources: default tiny; configurable via slurm.driver
  step_profile="driver"
  step_settings=$(parse_yaml_nested "slurm" "$step_profile" "$config_file_host")
  slurm_mem=""
  slurm_cpus=""
  slurm_time=""
  if [[ -n "$step_settings" ]]; then
    slurm_mem=$(parse_inline_dict "$step_settings" "mem")
    slurm_cpus=$(parse_inline_dict "$step_settings" "cpus")
    slurm_time=$(parse_inline_dict "$step_settings" "time")
  fi
  slurm_mem="${slurm_mem:-1g}"
  slurm_cpus="${slurm_cpus:-1}"
  slurm_time="${slurm_time:-2:00:00}"

  job_name="pgs_${sumstat_name}_driver"
  log_dir="${outdir_host}/sumstats/${sumstat_name}/logs/slurm"
  mkdir -p "$log_dir"

  run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps ${steps_arg} -i ${infold} --_driver-run"
  [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
  [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"
  [[ -n "$force_mode" ]] && run_cmd="${run_cmd} --force"
  [[ "$do_cleanup" == true ]] && run_cmd="${run_cmd} --cleanup"
  run_cmd="${run_cmd} --methods $(methods_to_csv "$CFG_METHODS")"

  sbatch_args=(--parsable)
  sbatch_args+=(--mem="${slurm_mem}")
  sbatch_args+=(--cpus-per-task="${slurm_cpus}")
  sbatch_args+=(--time="${slurm_time}")
  sbatch_args+=(--job-name="${job_name}")
  sbatch_args+=(--output="${log_dir}/${job_name}_%j.out")
  sbatch_args+=(--error="${log_dir}/${job_name}_%j.err")
  [[ -n "$slurm_account" ]] && sbatch_args+=(--account="${slurm_account}")
  [[ -n "$slurm_partition" ]] && sbatch_args+=(--partition="${slurm_partition}")
  [[ -n "${SLURM_DEPENDENCY:-}" ]] && sbatch_args+=(--dependency="${SLURM_DEPENDENCY}")
  sbatch_args+=(--wrap="${run_cmd}")

  echo "Submitting SLURM driver job..."
  echo "  Job name: ${job_name}"
  echo "  $(format_sbatch_settings "none")"
  echo "  Resources: mem=${slurm_mem}, cpus=${slurm_cpus}, time=${slurm_time}"
  [[ -n "${SLURM_DEPENDENCY:-}" ]] && echo "  Dependency: ${SLURM_DEPENDENCY}"
  echo "  Command: ${run_cmd}"
  echo ""

  jobid=$(sbatch "${sbatch_args[@]}")
  echo "Submitted: ${jobid}"
  echo "Driver log: ${log_dir}/${job_name}_${jobid}.out"
  echo "Driver err: ${log_dir}/${job_name}_${jobid}.err"
  echo "Watch: tail -f ${log_dir}/${job_name}_${jobid}.out"
  exit 0
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

# Liftover reference is required for dual-position mapfile
if [[ -z "$cfg_liftover_reference" ]]; then
  >&2 echo "Error: liftover_reference not found in config file"
  >&2 echo "This file is required for the dual-position variant mapfile."
  >&2 echo "Expected: references/liftover/dbsnp_cleansumstat_reference_GRCh37_GRCh38.txt.gz"
  exit 1
fi

if [[ ! -f "$cfg_liftover_reference" ]]; then
  >&2 echo "Error: Liftover reference file not found: $cfg_liftover_reference"
  exit 1
fi

# Validate genotype_build value
if [[ "$cfg_genotype_build" != "GRCh37" ]] && [[ "$cfg_genotype_build" != "GRCh38" ]]; then
  >&2 echo "Error: genotype_build must be 'GRCh37' or 'GRCh38', got: $cfg_genotype_build"
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
  if [[ -z "$steps_arg" ]] || [[ "$steps_arg" == *"sumstat"* ]] || [[ "$steps_arg" == *"weights"* ]] || [[ "$steps_arg" == *"score"* ]]; then
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

# Resolve liftover reference path
liftover_reference_host=$(realpath "${cfg_liftover_reference}")

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

# Extract sumstat name from input path.
# IMPORTANT: output folder should match the input folder name exactly.
# (e.g. input: /.../sumstat_5759  -> output: sumstats/sumstat_5759/)
sumstat_name=""
if [[ -n "$infold" ]]; then
    sumstat_name=$(basename "$infold_host")
fi

################################################################################
# Prerequisite checking
################################################################################

# Check prerequisites based on requested steps
needs_prep_check=1
if [[ "$steps_arg" == "prep" ]]; then
  needs_prep_check=0
else
  IFS="," read -ra steps_list <<< "$steps_arg"
  needs_prep_check=0
  for step in "${steps_list[@]}"; do
    if [[ ! "$step" =~ ^prep(-genotypes|-ldref|-inclusion-list|-inclusion-combine)?$ ]]; then
      needs_prep_check=1
      break
    fi
  done
fi

if [[ "$needs_prep_check" -eq 1 ]]; then
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

# Check sumstat prerequisite for weights
if [[ "$steps_arg" == *"weights"* ]] && [[ "$steps_arg" != *"sumstat"* ]]; then
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
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps sumstat,weights -i ${infold}"
    exit 1
  fi
fi

# Check weights (posteriors output) prerequisite for score
if [[ "$steps_arg" == *"score"* ]] && [[ "$steps_arg" != *"weights"* ]]; then
  missing_posteriors=$(check_posteriors_exists "$outdir_host" "$sumstat_name")
  if [[ $? -ne 0 ]]; then
    >&2 echo "Error: Weights (posteriors) calculation not completed."
    >&2 echo ""
    >&2 echo "Missing:"
    >&2 echo "$missing_posteriors"
    >&2 echo ""
    >&2 echo "Run weights step first:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps weights -i ${infold}"
    >&2 echo ""
    >&2 echo "Or include weights in your steps:"
    >&2 echo "  ./pgscalculator-v2.sh --config ${config_file} --steps weights,score -i ${infold}"
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
liftover_container="/pgscalculator/references/liftover_reference.gz"

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
# Use a unique filename to avoid collisions between concurrent runs / SLURM array tasks.
# IMPORTANT: array tasks run in parallel and may pass different --_chr values; if they share the same
# config_container.yaml path, they can overwrite each other and end up running the wrong chromosome.
config_base_dir=""
if [[ -n "${sumstat_name:-}" ]]; then
  config_base_dir="${outdir_host}/sumstats/${sumstat_name}/tmp"
else
  config_base_dir="${outdir_host}/prep/tmp"
fi
mkdir -p "$config_base_dir" 2>/dev/null || true
config_yaml_host="$(mktemp "${config_base_dir}/config_container.XXXXXX.yaml")"
# Map host outdir path -> container outdir path
config_yaml_container="${outdir_container}${config_yaml_host#${outdir_host}}"

# Copy original config and update paths for container
cat > "${config_yaml_host}" << EOF
# pgscalculator v2.1.0 - Auto-generated config for container
# Generated from: ${config_file_host}
input: ${indir_container}
outdir: ${outdir_container}
lddir: ${lddir_container}
genodir: ${genodir_container}
genofile: ${genofile_container}
genotype_build: ${cfg_genotype_build}
liftover_reference: ${liftover_container}
EOF

# Copy parameters from user's config (skip path keys and references section - handled separately)
awk '
  BEGIN { in_references = 0 }
  /^references:/ { in_references = 1; next }
  /^[a-zA-Z]/ && in_references { in_references = 0 }
  in_references { next }
  !/^(ld_reference|genotypes|genotype_manifest|outdir|input|lddir|genodir|genofile|genotype_build|liftover_reference):/ {
    print
  }
' "$config_file_host" >> "${config_yaml_host}"

# Ensure container config reflects the resolved active methods (CLI override wins).
sed -i '/^methods:/d' "${config_yaml_host}"
echo "methods: [$(methods_to_csv "$CFG_METHODS")]" >> "${config_yaml_host}"

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

if [[ -n "$force_mode" ]]; then
  cli_cmd="${cli_cmd} ${force_mode}"
fi

cli_cmd="${cli_cmd} --methods $(methods_to_csv "$CFG_METHODS")"
if [[ -n "$method_override" ]]; then
  cli_cmd="${cli_cmd} --method ${method_override}"
fi

################################################################################
# Build mount options
################################################################################
mount_opts=""
# Bind the *host* pgscalculator code into the container so the pipeline logic
# matches the wrapper you are running (and so fixes don't require rebuilding the SIF).
# This intentionally overrides the image's /pgscalculator tree.
mount_opts="${mount_opts} ${mountflag} ${project_dir}:/pgscalculator"
mount_opts="${mount_opts} ${mountflag} ${outdir_host}:${outdir_container}"
mount_opts="${mount_opts} ${mountflag} ${lddir_host}:${lddir_container}"

# Ensure we have a writable temp location with enough space, and bind it as /tmp
# inside the container. Many tools (e.g., sort) spill temporary files to /tmp.
# For sumstat-specific runs, keep tmp inside the sumstat folder for containment.
if [[ -n "${sumstat_name:-}" ]]; then
  tmpdir_host="${outdir_host}/sumstats/${sumstat_name}/tmp"
else
  tmpdir_host="${outdir_host}/prep/tmp"
fi
mkdir -p "${tmpdir_host}"
mount_opts="${mount_opts} ${mountflag} ${tmpdir_host}:/tmp"

if [[ -n "$infold" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${infold_host}:${indir_container}"
fi

if [[ -n "$genodir_host" ]] && [[ -d "$genodir_host" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${genodir_host}:${genodir_container}"
fi

if [[ -n "$genodir2_host" ]] && [[ -d "$genodir2_host" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${genodir2_host}:${genodir2_container}"
fi

# Mount liftover reference file (required for dual-position mapfile)
if [[ -n "$liftover_reference_host" ]] && [[ -f "$liftover_reference_host" ]]; then
  mount_opts="${mount_opts} ${mountflag} ${liftover_reference_host}:${liftover_container}"
fi

# Mount INFO file if provided (or allow "false" to explicitly disable)
if [[ -n "$cfg_info_file" ]] && [[ "${cfg_info_file,,}" == "false" ]]; then
  echo "info_file: false" >> "${config_yaml_host}"
elif [[ -n "$cfg_info_file" ]] && [[ -f "$cfg_info_file" ]]; then
  info_file_host=$(realpath "$cfg_info_file")
  info_file_container="/pgscalculator/references/info_scores.tsv"
  mount_opts="${mount_opts} ${mountflag} ${info_file_host}:${info_file_container}"
  # Add to config
  echo "info_file: ${info_file_container}" >> "${config_yaml_host}"
fi

# Mount MAF file if provided (or allow "false" to explicitly disable)
if [[ -n "$cfg_maf_file" ]] && [[ "${cfg_maf_file,,}" == "false" ]]; then
  echo "maf_file: false" >> "${config_yaml_host}"
elif [[ -n "$cfg_maf_file" ]] && [[ -f "$cfg_maf_file" ]]; then
  maf_file_host=$(realpath "$cfg_maf_file")
  maf_file_container="/pgscalculator/references/maf.tsv"
  mount_opts="${mount_opts} ${mountflag} ${maf_file_host}:${maf_file_container}"
  # Add to config
  echo "maf_file: ${maf_file_container}" >> "${config_yaml_host}"
fi

# Mount user inclusion lists (gt/ss/ld) if provided
if [[ -n "$cfg_inclusion_gt" || -n "$cfg_inclusion_ss" || -n "$cfg_inclusion_ld" ]]; then
  echo "filters:" >> "${config_yaml_host}"
  echo "  inclusion_list:" >> "${config_yaml_host}"
  if [[ -n "$cfg_inclusion_gt" ]] && [[ -f "$cfg_inclusion_gt" ]]; then
    inclusion_gt_host=$(realpath "$cfg_inclusion_gt")
    inclusion_gt_container="/pgscalculator/references/inclusion_gt.tsv"
    mount_opts="${mount_opts} ${mountflag} ${inclusion_gt_host}:${inclusion_gt_container}"
    echo "    gt: ${inclusion_gt_container}" >> "${config_yaml_host}"
  fi
  if [[ -n "$cfg_inclusion_ss" ]] && [[ -f "$cfg_inclusion_ss" ]]; then
    inclusion_ss_host=$(realpath "$cfg_inclusion_ss")
    inclusion_ss_container="/pgscalculator/references/inclusion_ss.tsv"
    mount_opts="${mount_opts} ${mountflag} ${inclusion_ss_host}:${inclusion_ss_container}"
    echo "    ss: ${inclusion_ss_container}" >> "${config_yaml_host}"
  fi
  if [[ -n "$cfg_inclusion_ld" ]] && [[ -f "$cfg_inclusion_ld" ]]; then
    inclusion_ld_host=$(realpath "$cfg_inclusion_ld")
    inclusion_ld_container="/pgscalculator/references/inclusion_ld.tsv"
    mount_opts="${mount_opts} ${mountflag} ${inclusion_ld_host}:${inclusion_ld_container}"
    echo "    ld: ${inclusion_ld_container}" >> "${config_yaml_host}"
  fi
fi

################################################################################
# Execute
################################################################################
echo "Running pgscalculator v2.1.0 in Singularity"
echo "Config: ${config_file_host}"
echo "Output: ${outdir_host}"
  echo "Command: ${cli_cmd}"

  # Force temp usage inside container to /tmp (which we bind to ${outdir_host}/tmp).
  # With --cleanenv, set via APPTAINERENV_*
  APPTAINERENV_TMPDIR=/tmp \
  APPTAINERENV_TMP=/tmp \
  singularity run \
     --contain \
     --cleanenv \
     ${mount_flags} \
  ${mount_opts} \
     "${runimage}" \
     ${cli_cmd}

exit_code=$?

# Optional cleanup for local (non-driver) runs.
# For driver jobs, cleanup happens inside the driver block after arrays/finalization.
if [[ "$exit_code" -eq 0 ]] && [[ "$driver_run" != true ]] && [[ "$do_cleanup" == true ]]; then
  echo "Cleanup enabled: removing work/ and tmp/ for ${sumstat_name:-"(no sumstat)"}"
  cleanup_work_and_tmp "$outdir_host" "${sumstat_name:-}"
fi

exit "$exit_code"
