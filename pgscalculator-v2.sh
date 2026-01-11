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
  echo "                   For per-sumstat steps (sumstat/posteriors/score), this submits a"
  echo "                   lightweight *driver job* that runs sumstat and launches/monitors"
  echo "                   chromosome-parallel arrays for posteriors/score."
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
  echo "  # Optional: SLURM settings for --sbatch / --sbatch-array"
  echo "  slurm:"
  echo "    account: my_account"
  echo "    max_parallel: 22"
  echo "    prep:       { mem: 10g, cpus: 6, time: '1:00:00' }"
  echo "    posteriors: { mem: 20g, cpus: 6, time: '2:00:00', max_parallel: 22 }"
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
driver_run=false

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
    --_driver-run)
      # Internal flag: run inside a driver job; orchestrates arrays, does not submit itself.
      driver_run=true
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

# (legacy) no-op: --sbatch-array removed; keep block absent

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

if [[ "$driver_run" == true ]]; then
  # Support: any combination of:
  #   --steps sumstat
  #   --steps posteriors
  #   --steps score
  #   --steps sumstat,posteriors,score
  # Order is enforced: sumstat -> posteriors -> score
  if [[ -z "${steps_arg:-}" ]]; then
    >&2 echo "Error: --_driver-run requires --steps (sumstat, posteriors, score, or combinations thereof)"
    exit 1
  fi
  has_sumstat=false
  has_posteriors=false
  has_score=false
  IFS=',' read -r -a _sbatch_steps <<< "$steps_arg"
  for _s in "${_sbatch_steps[@]}"; do
    _s="$(echo "$_s" | awk '{$1=$1;print}')"
    [[ -z "$_s" ]] && continue
    if [[ "$_s" == "prep" ]]; then
      >&2 echo "Error: prep must be run on its own (do not include prep in driver jobs)"
      exit 1
    fi
    if [[ "$_s" == "sumstat" ]]; then
      has_sumstat=true
    elif [[ "$_s" == "posteriors" ]]; then
      has_posteriors=true
    elif [[ "$_s" == "score" ]]; then
      has_score=true
    else
      >&2 echo "Error: driver mode only supports --steps sumstat, posteriors, score (or combinations) (got: '${steps_arg}')"
      exit 1
    fi
  done
  if [[ "$has_sumstat" != true && "$has_posteriors" != true && "$has_score" != true ]]; then
    >&2 echo "Error: --sbatch-array requires --steps to include sumstat and/or posteriors and/or score"
    exit 1
  fi

  # Need outdir for logs and chromosome list file
  if [[ -z "$cfg_outdir" ]]; then
    >&2 echo "Error: outdir not found in config file and -o not provided"
    exit 1
  fi
  mkdir -p "${cfg_outdir}"
  outdir_host=$(realpath "${cfg_outdir}")

  # In --sbatch-array mode we exit before the later "Resolve paths" block.
  # Compute sumstat_name here so watcher sanity-checks can run.
  if [[ -z "${infold:-}" ]]; then
    >&2 echo "Error: --sbatch-array requires -i <sumstat_dir>"
    exit 1
  fi
infold_host=$(realpath "${infold}")
  if [[ ! -d "$infold_host" ]]; then
    >&2 echo "Error: Input directory doesn't exist: $infold_host"
    exit 1
  fi
  sumstat_name=$(basename "$infold_host")

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

  # Keep SLURM logs contained within the sumstat folder to avoid collisions across runs.
  # (prep is shared; sumstat/posteriors/score are per-sumstat)
  log_dir="${outdir_host}/sumstats/${sumstat_name}/logs/slurm"
  mkdir -p "$log_dir"

  watch_array() {
    local array_jobid="$1"
    local chr_file="$2"
    local chr_count="$3"

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
    declare -A acct_miss
    declare -A task_state
    declare -A task_failed
    start_ts=$(date +%s)

    while [[ ${#finished[@]} -lt $chr_count ]]; do
      # Snapshot current tasks in queue
      mapfile -t sq_lines < <(squeue -h -j "${array_jobid}" -o "%i|%T" 2>/dev/null || true)

      declare -A in_queue
      for line in "${sq_lines[@]}"; do
        jid="${line%%|*}"
        state="${line##*|}"
        # jid can be like 12345_7
        if [[ "$jid" =~ ^${array_jobid}_[0-9]+$ ]]; then
          idx="${jid#${array_jobid}_}"
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
              finished["$idx"]=1
              task_state["$idx"]="$st"
              if [[ "$st" != COMPLETED* ]]; then
                task_failed["$idx"]=1
              fi
              echo "  Finished: task=${idx} chr=${chr} state=${st} elapsed=${elapsed_fmt} time=$(date)"
            else
              # sacct can lag behind job completion; don't hang forever.
              acct_miss["$idx"]=$(( ${acct_miss["$idx"]:-0} + 1 ))
              if [[ ${acct_miss["$idx"]} -ge 3 ]]; then
                finished["$idx"]=1
                task_state["$idx"]="UNKNOWN"
                echo "  Finished: task=${idx} chr=${chr} (sacct not yet available) time=$(date)"
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
              finished["$idx"]=1
              task_state["$idx"]="$st"
              if [[ "$st" != COMPLETED* ]]; then
                task_failed["$idx"]=1
              fi
              echo "  Finished: task=${idx} chr=${chr} state=${st} elapsed=${elapsed_fmt} time=$(date)"
            else
              finished["$idx"]=1
              task_state["$idx"]="UNKNOWN"
              echo "  Finished: task=${idx} chr=${chr} (array ended; sacct unavailable) time=$(date)"
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

    # Be conservative: if some tasks are still UNKNOWN (no sacct), treat as failure.
    local unknown_count=0
    for ((idx=1; idx<=chr_count; idx++)); do
      if [[ "${task_state[$idx]:-}" == "UNKNOWN" ]]; then
        unknown_count=$((unknown_count + 1))
      fi
    done
    if [[ $unknown_count -gt 0 ]]; then
      echo "Array result: UNKNOWN_TASKS=${unknown_count}/${chr_count}"
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

    step_settings=$(parse_yaml_nested "slurm" "$step_profile" "$config_file_host")
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
    slurm_mem="${slurm_mem:-20g}"
    slurm_cpus="${slurm_cpus:-8}"
    slurm_time="${slurm_time:-2:00:00}"

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
    if [[ "$step_profile" == "score" && -n "${sumstat_name:-}" ]]; then
      mapped_new="${outdir_host}/sumstats/${sumstat_name}/intermediates/posteriors_mapped"
      mapped_old="${outdir_host}/sumstats/${sumstat_name}/posteriors_mapped"
      mapped_dir=""
      if [[ -d "$mapped_new" ]]; then
        mapped_dir="$mapped_new"
      elif [[ -d "$mapped_old" ]]; then
        mapped_dir="$mapped_old"
      fi

      if [[ -n "$mapped_dir" ]]; then
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
    # - posteriors: safe (calc-posteriors + format-posteriors are chr-parallel)
    # - score: only calc-score is chr-parallel; combine-scores/finalize-output must be run once after
    steps_arg_for_task="${step_profile}"
    if [[ "$step_profile" == "score" ]]; then
      steps_arg_for_task="calc-score"
    fi

    run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps ${steps_arg_for_task}"
    [[ -n "$infold" ]] && run_cmd="${run_cmd} -i ${infold}"
    [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"

    task_wrap="CHR=\$(sed -n \"\${SLURM_ARRAY_TASK_ID}p\" \"${chr_file}\"); \
if [[ -z \"\$CHR\" ]]; then echo \"Error: could not resolve chromosome for task \$SLURM_ARRAY_TASK_ID\" >&2; exit 1; fi; \
echo \"[INFO] Starting ${step_profile} chr\${CHR} at \$(date)\"; \
${run_cmd} --chr \"\$CHR\"; \
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
    if [[ "$step_profile" == "score" ]]; then
      echo "  Note: array mode runs 'calc-score' only; combine/finalize will run once after the score array finishes."
    fi
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
      base_sumstat_out="${outdir_host}/sumstats/${sumstat_name}/intermediates"
      if [[ "$step_profile" == "posteriors" ]]; then
        n_post=$(ls "${base_sumstat_out}/posteriors"/chr*.snpRes 2>/dev/null | wc -l | awk '{print $1}')
        n_mapped=$(ls "${base_sumstat_out}/posteriors_mapped"/chr*.snpRes 2>/dev/null | wc -l | awk '{print $1}')
        if [[ "$n_post" -lt "$step_chr_count" || "$n_mapped" -lt "$step_chr_count" ]]; then
          >&2 echo "Warning: posteriors array finished but outputs are missing (continuing)."
          >&2 echo "  Expected >=${step_chr_count} files in:"
          >&2 echo "    - ${base_sumstat_out}/posteriors/chr*.snpRes   (found ${n_post})"
          >&2 echo "    - ${base_sumstat_out}/posteriors_mapped/chr*.snpRes (found ${n_mapped})"
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      elif [[ "$step_profile" == "score" ]]; then
        n_scores=$(ls "${base_sumstat_out}/scores"/chr*.sscore 2>/dev/null | wc -l | awk '{print $1}')
        if [[ "$n_scores" -lt "$step_chr_count" ]]; then
          >&2 echo "Warning: score array finished but outputs are missing (continuing)."
          >&2 echo "  Expected >=${step_chr_count} files in: ${base_sumstat_out}/scores/chr*.sscore (found ${n_scores})"
          >&2 echo "Check logs under: ${log_dir}/"
        fi
      fi
    fi
    echo ""
  }

  # Always run posteriors before score if both requested
  # In driver mode, we run sumstat directly (within this driver job), then launch arrays.
  run_base="${project_dir}/pgscalculator-v2.sh --config ${config_file_host}"
  [[ -n "$infold" ]] && run_base="${run_base} -i ${infold}"
  [[ -n "$outdir" ]] && run_base="${run_base} -o ${outdir}"
  [[ -n "$chromosomes" ]] && run_base="${run_base} --chr ${chromosomes}"
  [[ -n "$devmode" ]] && run_base="${run_base} -d"

  if [[ "$has_sumstat" == true ]]; then
    echo "Running sumstat inside driver job..."
    eval "${run_base} --steps sumstat"

    # Gate downstream arrays on the expected sumstat outputs existing.
    # This prevents submitting 22 tasks that all fail immediately due to missing inputs.
    if ! check_sumstat_exists "$outdir_host" "$sumstat_name" >/dev/null; then
      >&2 echo "Error: sumstat step finished but expected outputs are missing."
      >&2 echo "Missing:"
      check_sumstat_exists "$outdir_host" "$sumstat_name" || true
      exit 1
    fi
  fi
  if [[ "$has_posteriors" == true ]]; then
    submit_array_for_step "posteriors"
  fi
  if [[ "$has_score" == true ]]; then
    submit_array_for_step "score"

    # After score array finishes, run the non-parallel steps once to produce final outputs.
    follow_steps="combine-scores,finalize-output"
    echo "Running post-array finalization: ${follow_steps}"
    follow_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps ${follow_steps}"
    [[ -n "$infold" ]] && follow_cmd="${follow_cmd} -i ${infold}"
    [[ -n "$outdir" ]] && follow_cmd="${follow_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && follow_cmd="${follow_cmd} -d"
    eval "$follow_cmd"
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
    >&2 echo "  --steps sumstat,posteriors,score --sbatch -i <sumstat_dir>"
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
    # Prep job (shared across sumstats)
    step_profile="prep"
    step_settings=$(parse_yaml_nested "slurm" "$step_profile" "$config_file_host")
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

    job_name="pgs_prep"
    log_dir="${outdir_host}/prep/logs/slurm"
    mkdir -p "$log_dir"

    run_cmd="${project_dir}/pgscalculator-v2.sh --config ${config_file_host} --steps prep"
    [[ -n "$outdir" ]] && run_cmd="${run_cmd} -o ${outdir}"
    [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"

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
  [[ -n "$chromosomes" ]] && run_cmd="${run_cmd} --chr ${chromosomes}"
  [[ -n "$devmode" ]] && run_cmd="${run_cmd} -d"

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

  echo "Submitting SLURM driver job..."
  echo "  Job name: ${job_name}"
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
check_prep_exists() {
  local outdir="$1"
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
  
  if [[ -n "$missing" ]]; then
    echo "$missing"
    return 1
  fi
  return 0
}

check_sumstat_exists() {
  local outdir="$1"
  local sumstat="$2"
  
  # Prefer v2.1 intermediates/, but accept legacy locations for backwards compatibility.
  local formatted_new="${outdir}/sumstats/${sumstat}/intermediates/formatted/sumstat_formatted.tsv.gz"
  local formatted_old="${outdir}/sumstats/${sumstat}/formatted/sumstat_formatted.tsv.gz"
  local filtered_new="${outdir}/sumstats/${sumstat}/intermediates/filtered/sumstat_filtered.tsv.gz"
  local filtered_old="${outdir}/sumstats/${sumstat}/filtered/sumstat_filtered.tsv.gz"
  local filtered_chr_new_glob="${outdir}/sumstats/${sumstat}/intermediates/filtered/chr*_filtered.tsv"
  local filtered_chr_old_glob="${outdir}/sumstats/${sumstat}/filtered/chr*_filtered.tsv"
  
  if [[ ! -f "$formatted_new" && ! -f "$formatted_old" ]]; then
    echo "  - Formatted sumstat (expected): ${formatted_new}"
    echo "    (legacy accepted): ${formatted_old}"
    return 1
  fi
  
  if [[ ! -f "$filtered_new" && ! -f "$filtered_old" ]]; then
    # Some workflows may only need per-chromosome filtered files (chrN_filtered.tsv).
    # Accept those as an alternative prereq for posteriors.
    if ! compgen -G "$filtered_chr_new_glob" >/dev/null 2>&1 && ! compgen -G "$filtered_chr_old_glob" >/dev/null 2>&1; then
      echo "  - Filtered sumstat (expected): ${filtered_new}"
      echo "    (legacy accepted): ${filtered_old}"
      echo "    (alt accepted): ${filtered_chr_new_glob}"
      return 1
    fi
  fi
  return 0
}

check_posteriors_exists() {
  local outdir="$1"
  local sumstat="$2"
  
  local post_new="${outdir}/sumstats/${sumstat}/intermediates/posteriors"
  local post_old="${outdir}/sumstats/${sumstat}/posteriors"
  local mapped_new="${outdir}/sumstats/${sumstat}/intermediates/posteriors_mapped"
  local mapped_old="${outdir}/sumstats/${sumstat}/posteriors_mapped"
  
  # calc-posteriors output
  if [[ ( ! -d "$post_new" || -z "$(ls -A "$post_new" 2>/dev/null)" ) && ( ! -d "$post_old" || -z "$(ls -A "$post_old" 2>/dev/null)" ) ]]; then
    echo "  - Posteriors (expected): ${post_new}/"
    echo "    (legacy accepted): ${post_old}/"
    return 1
  fi
  
  # format-posteriors output (required for scoring)
  if [[ ( ! -d "$mapped_new" || -z "$(ls -A "$mapped_new" 2>/dev/null)" ) && ( ! -d "$mapped_old" || -z "$(ls -A "$mapped_old" 2>/dev/null)" ) ]]; then
    echo "  - Posteriors mapped (expected): ${mapped_new}/"
    echo "    (legacy accepted): ${mapped_old}/"
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
  tmpdir_host="${outdir_host}/tmp"
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

################################################################################
# Execute
################################################################################
echo "Running pgscalculator v2.1.0 in Singularity"
echo "Config: ${config_file_host}"
echo "Output: ${outdir_host}"
  echo "Command: ${cli_cmd}"

  # Force temp usage inside container to /tmp (which we bind to ${outdir_host}/tmp).
  # With --cleanenv, set via SINGULARITYENV_*
  SINGULARITYENV_TMPDIR=/tmp \
  SINGULARITYENV_TMP=/tmp \
  singularity run \
     --contain \
     --cleanenv \
     ${mount_flags} \
  ${mount_opts} \
     "${runimage}" \
     ${cli_cmd}
