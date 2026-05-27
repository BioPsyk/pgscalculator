#!/bin/bash
# pgscalculator v2 - Common library functions
# Shared functions for logging, dependency checking, config parsing, and utilities

set -euo pipefail

# =============================================================================
# GLOBAL VARIABLES
# =============================================================================

# Colors for terminal output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[0;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

# Default values
VERBOSE=${VERBOSE:-0}
DRY_RUN=${DRY_RUN:-0}

# =============================================================================
# LOGGING FUNCTIONS
# =============================================================================

log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*" >&2
}

log_debug() {
    if [[ "${VERBOSE}" -eq 1 ]]; then
        echo -e "${BLUE}[DEBUG]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $*"
    fi
}

log_step() {
    echo -e "\n${GREEN}==>${NC} $*"
}

log_substep() {
    echo -e "  ${BLUE}-->${NC} $*"
}

# =============================================================================
# DEPENDENCY CHECKING FUNCTIONS
# =============================================================================

require_file() {
    local filepath="$1"
    local hint="${2:-}"
    
    if [[ ! -f "$filepath" ]]; then
        log_error "Required file not found: $filepath"
        if [[ -n "$hint" ]]; then
            log_error "Hint: $hint"
        fi
        exit 1
    fi
    log_debug "Found required file: $filepath"
}

require_dir() {
    local dirpath="$1"
    local hint="${2:-}"
    
    if [[ ! -d "$dirpath" ]]; then
        log_error "Required directory not found: $dirpath"
        if [[ -n "$hint" ]]; then
            log_error "Hint: $hint"
        fi
        exit 1
    fi
    log_debug "Found required directory: $dirpath"
}

require_command() {
    local cmd="$1"
    local hint="${2:-}"
    
    if ! command -v "$cmd" &> /dev/null; then
        log_error "Required command not found: $cmd"
        if [[ -n "$hint" ]]; then
            log_error "Hint: $hint"
        fi
        exit 1
    fi
    log_debug "Found required command: $cmd"
}

check_step_completed() {
    local step_dir="$1"
    local marker_file="${step_dir}/.completed"
    
    # If FORCE=1 is set, always return false (step not completed)
    if [[ "${FORCE:-0}" -eq 1 ]]; then
        return 1
    fi
    
    if [[ -f "$marker_file" ]]; then
        return 0
    fi
    return 1
}

mark_step_completed() {
    local step_dir="$1"
    local marker_file="${step_dir}/.completed"
    
    date '+%Y-%m-%d %H:%M:%S' > "$marker_file"
    log_debug "Marked step as completed: $step_dir"
}

# =============================================================================
# CONFIG PARSING FUNCTIONS (Simple YAML-like parsing)
# =============================================================================

# Parse a simple YAML config file and export variables
# Supports: key: value and nested key.subkey: value (flattened with _)
parse_config() {
    local config_file="$1"
    local prefix="${2:-CFG}"
    
    require_file "$config_file" "Config file is required"
    
    log_debug "Parsing config file: $config_file"
    
    # Simple YAML parser - handles basic key: value pairs
    # Does not support complex YAML features like arrays, multiline, etc.
    local current_section=""
    local current_subsection=""
    
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Skip comments and empty lines
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line// }" ]] && continue
        
        # Check for section (key with nested values)
        if [[ "$line" =~ ^([a-zA-Z_][a-zA-Z0-9_]*):$ ]]; then
            # methods: uses list syntax handled by load_active_methods(), not nested keys
            [[ "${BASH_REMATCH[1]}" == "methods" ]] && continue
            current_section="${BASH_REMATCH[1]}"
            current_subsection=""
            continue
        fi

        # Check for indented subsection (e.g., "  inclusion_list:")
        if [[ "$line" =~ ^[[:space:]]{2}([a-zA-Z_][a-zA-Z0-9_]*):$ ]]; then
            current_subsection="${BASH_REMATCH[1]}"
            continue
        fi
        
        # Check for indented key: value (nested under section)
        if [[ "$line" =~ ^[[:space:]]+([a-zA-Z_][a-zA-Z0-9_]*):\ *(.*) ]]; then
            local key="${BASH_REMATCH[1]}"
            local value="${BASH_REMATCH[2]}"
            # Remove quotes if present
            value="${value%\"}"
            value="${value#\"}"
            value="${value%\'}"
            value="${value#\'}"
            
            if [[ -n "$current_section" && -n "$current_subsection" ]]; then
                local var_name="${prefix}_${current_section}_${current_subsection}_${key}"
            elif [[ -n "$current_section" ]]; then
                local var_name="${prefix}_${current_section}_${key}"
            else
                local var_name="${prefix}_${key}"
            fi
            var_name="${var_name^^}"  # Uppercase
            export "${var_name}=${value}"
            log_debug "Config: ${var_name}=${value}"
            continue
        fi
        
        # Check for top-level key: value
        if [[ "$line" =~ ^([a-zA-Z_][a-zA-Z0-9_]*):\ *(.*) ]]; then
            current_section=""
            local key="${BASH_REMATCH[1]}"
            # methods: [a, b] or block list — handled by load_active_methods()
            [[ "$key" == "methods" ]] && continue
            local value="${BASH_REMATCH[2]}"
            # Remove quotes if present
            value="${value%\"}"
            value="${value#\"}"
            value="${value%\'}"
            value="${value#\'}"
            
            local var_name="${prefix}_${key}"
            var_name="${var_name^^}"  # Uppercase
            export "${var_name}=${value}"
            log_debug "Config: ${var_name}=${value}"
        fi
    done < "$config_file"
}

# Get a config value with default fallback
get_config() {
    local key="$1"
    local default="${2:-}"
    local prefix="${3:-CFG}"
    
    local var_name="${prefix}_${key}"
    var_name="${var_name^^}"
    
    local value="${!var_name:-$default}"
    echo "$value"
}

# Parse YAML list: block-style (key:\n  - a) or inline (key: [a, b])
parse_yaml_list() {
    local key="$1"
    local file="$2"
    awk -v key="$key" '
        BEGIN { in_list = 0 }
        $0 ~ "^"key":[[:space:]]*\\[" {
            line = $0
            sub("^"key":[[:space:]]*\\[", "", line)
            sub("\\][[:space:]]*$", "", line)
            gsub(/'"'"'/, "", line)
            n = split(line, parts, /,[[:space:]]*/)
            for (i = 1; i <= n; i++) {
                gsub(/^[ \t]+|[ \t]+$/, "", parts[i])
                if (parts[i] != "") print parts[i]
            }
            exit
        }
        $0 ~ "^"key":" { in_list = 1; next }
        in_list && /^  - / { gsub(/^  - */, ""); gsub(/[ \t]+$/, ""); print; next }
        in_list && /^[^ ]/ { exit }
    ' "$file"
}

# Parse nested YAML values (supports one or two levels under a section).
parse_yaml_nested() {
    local section="$1"
    local key="$2"
    local file="$3"
    awk -v section="$section" -v key="$key" '
        BEGIN {
            in_section = 0
            in_subsection = 0
            n = split(key, parts, /\./)
            key1 = parts[1]
            key2 = (n >= 2 ? parts[2] : "")
        }
        $0 ~ "^"section":[[:space:]]*$" {
            in_section = 1
            in_subsection = 0
            next
        }
        in_section && /^[^[:space:]][^:]*:[[:space:]]*$/ {
            in_section = 0
            in_subsection = 0
        }
        !in_section { next }

        n == 1 && $0 ~ "^[[:space:]]{2}"key1":[[:space:]]*" {
            line = $0
            sub("^[[:space:]]{2}"key1":[[:space:]]*", "", line)
            gsub(/[{}]/, "", line)
            gsub(/^[ \t]+|[ \t]+$/, "", line)
            print line
            exit
        }

        n >= 2 && $0 ~ "^[[:space:]]{2}"key1":[[:space:]]*$" {
            in_subsection = 1
            next
        }
        n >= 2 && in_subsection && /^[[:space:]]{2}[a-zA-Z_][a-zA-Z0-9_]*:[[:space:]]*$/ && $0 !~ "^[[:space:]]{2}"key1":[[:space:]]*$" {
            in_subsection = 0
        }
        n >= 2 && in_subsection && $0 ~ "^[[:space:]]{4}"key2":[[:space:]]*" {
            line = $0
            sub("^[[:space:]]{4}"key2":[[:space:]]*", "", line)
            gsub(/[{}]/, "", line)
            gsub(/^[ \t]+|[ \t]+$/, "", line)
            print line
            exit
        }
    ' "$file"
}

# Normalize comma/space-separated method names to a space-separated lowercase list.
normalize_methods_list() {
    local raw="$1"
    local out=""
    local tok

    raw="${raw//,/ }"
    for tok in $raw; do
        tok="${tok,,}"
        tok="${tok// /}"
        [[ -z "$tok" ]] && continue
        case "$tok" in
            sbayesr|ldpred2) out="${out:+${out} }${tok}" ;;
            *)
                log_error "Unknown method: '${tok}' (expected: sbayesr, ldpred2)"
                return 1
                ;;
        esac
    done

    if [[ -z "$out" ]]; then
        out="sbayesr"
    fi
    echo "$out"
}

# Return 0 when method is in the active methods list (space-separated).
has_method() {
    local method="${1,,}"
    local methods_list="${2:-${CFG_METHODS:-sbayesr}}"
    local m

    for m in $methods_list; do
        [[ "${m,,}" == "$method" ]] && return 0
    done
    return 1
}

# Default map basename inside ldpred2.ld_dir from ldpred2.ld_variant_set.
ldpred2_default_map_basename() {
    local set="${1:-${CFG_LDPRED2_LD_VARIANT_SET:-hm3_plus}}"

    case "$set" in
        hm3) echo "map_hm3.rds" ;;
        hm3_plus) echo "map_hm3_plus.rds" ;;
        *)
            log_error "ldpred2.ld_variant_set: '${set}' is not recognised (expected hm3 or hm3_plus)"
            return 1
            ;;
    esac
}

# Resolve ldpred2.ld_meta_file (explicit config value or default under ld_dir).
resolve_ldpred2_ld_meta_file() {
    if [[ -n "${CFG_LDPRED2_LD_META_FILE:-}" ]]; then
        echo "${CFG_LDPRED2_LD_META_FILE}"
        return 0
    fi
    if [[ -z "${CFG_LDPRED2_LD_DIR:-}" ]]; then
        return 1
    fi

    local map_basename
    map_basename=$(ldpred2_default_map_basename) || return 1
    echo "${CFG_LDPRED2_LD_DIR%/}/${map_basename}"
}

# Fail fast when ldpred2 is active but required LDpred2 config is missing.
validate_active_methods_config() {
    local methods_list="${CFG_METHODS:-sbayesr}"

    if ! has_method ldpred2 "$methods_list"; then
        return 0
    fi

    if [[ -z "${CFG_LDPRED2_LD_DIR:-}" ]]; then
        log_error "methods includes 'ldpred2' but ldpred2.ld_dir is not set in config"
        log_error "Set ldpred2.ld_dir (and optionally ldpred2.ld_meta_file) or remove ldpred2 from methods"
        exit 1
    fi

    local meta_file
    meta_file=$(resolve_ldpred2_ld_meta_file) || {
        log_error "methods includes 'ldpred2' but ldpred2.ld_meta_file could not be resolved"
        log_error "Set ldpred2.ld_meta_file explicitly or fix ldpred2.ld_variant_set"
        exit 1
    }
    if [[ -z "$meta_file" ]]; then
        log_error "methods includes 'ldpred2' but ldpred2.ld_meta_file is not set and could not be defaulted"
        exit 1
    fi

    export CFG_LDPRED2_LD_META_FILE_RESOLVED="$meta_file"
    log_debug "Resolved ldpred2.ld_meta_file: ${meta_file}"
}

# Resolve active methods: CLI override > config methods: > default [sbayesr].
# Exports CFG_METHODS (space-separated) and runs validate_active_methods_config().
load_active_methods() {
    local methods_cli="${1:-}"
    local config_file="${2:-}"
    local raw_methods=""

    if [[ -n "$methods_cli" ]]; then
        raw_methods="$methods_cli"
    elif [[ -n "$config_file" && -f "$config_file" ]]; then
        raw_methods=$(parse_yaml_list "methods" "$config_file" | tr '\n' ' ')
    fi

    if [[ -z "${raw_methods// /}" ]]; then
        raw_methods="sbayesr"
    fi

    CFG_METHODS=$(normalize_methods_list "$raw_methods") || exit 1
    export CFG_METHODS
    log_debug "Active methods: ${CFG_METHODS}"
    validate_active_methods_config
}

# Convert space-separated methods to comma-separated (for --methods CLI).
methods_to_csv() {
    echo "${1:-${CFG_METHODS:-sbayesr}}" | tr ' ' ','
}

# SLURM step profile lookup with score_sbayesr/score_ldpred2 inheriting slurm.score.
resolve_slurm_step_settings() {
    local step_profile="$1"
    local config_file="$2"
    local step_settings=""

    step_settings=$(parse_yaml_nested "slurm" "$step_profile" "$config_file")
    if [[ -z "$step_settings" ]]; then
        case "$step_profile" in
            score_sbayesr|score_ldpred2)
                step_settings=$(parse_yaml_nested "slurm" "score" "$config_file")
                ;;
        esac
    fi
    echo "$step_settings"
}

# =============================================================================
# PATH AND DIRECTORY UTILITIES
# =============================================================================

# Get the output directory for a specific step
get_step_dir() {
    local outdir="$1"
    local step_name="$2"
    local sumstat_name="${3:-}"
    
    if [[ -n "$sumstat_name" ]]; then
        # Sumstat step directories live under work/ (per-sumstat containment)
        echo "${outdir}/sumstats/${sumstat_name}/work/${step_name}"
    else
        echo "${outdir}/prep/${step_name}"
    fi
}

# Ensure a directory exists, create if needed
ensure_dir() {
    local dirpath="$1"
    
    if [[ ! -d "$dirpath" ]]; then
        log_debug "Creating directory: $dirpath"
        mkdir -p "$dirpath"
    fi
}

# Get the prep directory
get_prep_dir() {
    local outdir="$1"
    echo "${outdir}/prep"
}

# Get the sumstat directory
get_sumstat_dir() {
    local outdir="$1"
    local sumstat_name="$2"
    echo "${outdir}/sumstats/${sumstat_name}"
}

# Get the work directory within a sumstat run directory (new canonical name)
get_sumstat_work_dir() {
    local sumstat_dir="$1"
    echo "${sumstat_dir}/work"
}

# Migrate step directories into sumstat_dir/work/<step> for v2.1+ output structure.
# Handles:
# - legacy: sumstat_dir/<step>
migrate_sumstat_step_dir() {
    local sumstat_dir="$1"
    local step="$2"

    local work_dir
    work_dir=$(get_sumstat_work_dir "$sumstat_dir")
    ensure_dir "$work_dir"

    local legacy_root="${sumstat_dir}/${step}"
    local new_dir="${work_dir}/${step}"

    if [[ -d "$legacy_root" ]] && [[ ! -d "$new_dir" ]]; then
        mv "$legacy_root" "$new_dir"
    fi
}

get_sumstat_step_dir() {
    local sumstat_dir="$1"
    local step="$2"
    local work_dir
    work_dir=$(get_sumstat_work_dir "$sumstat_dir")
    echo "${work_dir}/${step}"
}

# Rename legacy work/filtered/ -> work/filtered_sbayesr/ (Phase 1 §5.4).
# Idempotent: no-op when legacy dir is absent or the new dir already exists.
migrate_filtered_dirs() {
    local sumstat_dir="$1"

    local work_dir
    work_dir=$(get_sumstat_work_dir "$sumstat_dir")
    ensure_dir "$work_dir"

    local legacy_dir="${work_dir}/filtered"
    local sbayesr_dir="${work_dir}/filtered_sbayesr"

    if [[ -d "$legacy_dir" ]] && [[ ! -d "$sbayesr_dir" ]]; then
        log_info "Migrating legacy work/filtered/ -> work/filtered_sbayesr/"
        mv "$legacy_dir" "$sbayesr_dir"
    fi
}

# Migrate all known per-sumstat step directories into work/.
# Safe to call repeatedly.
migrate_sumstat_all_step_dirs() {
    local sumstat_dir="$1"
    for step in formatted filtered posteriors posteriors_mapped scores scores_combined; do
        migrate_sumstat_step_dir "$sumstat_dir" "$step"
    done
    migrate_filtered_dirs "$sumstat_dir"
}

# =============================================================================
# FILE UTILITIES
# =============================================================================

# Count lines in a file (excluding header)
count_data_lines() {
    local filepath="$1"
    local has_header="${2:-1}"
    
    if [[ ! -f "$filepath" ]]; then
        echo "0"
        return
    fi
    
    local total
    total=$(wc -l < "$filepath")
    
    if [[ "$has_header" -eq 1 ]]; then
        echo $((total - 1))
    else
        echo "$total"
    fi
}

# Check if file is gzipped
is_gzipped() {
    local filepath="$1"
    
    if [[ "$filepath" == *.gz ]]; then
        return 0
    fi
    return 1
}

# Cat file (handles gzipped files)
cat_file() {
    local filepath="$1"
    
    if is_gzipped "$filepath"; then
        zcat "$filepath"
    else
        cat "$filepath"
    fi
}

# =============================================================================
# CHROMOSOME UTILITIES
# =============================================================================

# Get list of chromosomes (default: 1-22, configurable via CFG_CHROMOSOMES)
# CFG_CHROMOSOMES can be:
#   - "21,22" or "21 22" - specific chromosomes
#   - "21-22" - range of chromosomes
#   - not set - defaults to 1-22
get_chromosomes() {
    local chr_spec="${CFG_CHROMOSOMES:-}"
    
    if [[ -z "$chr_spec" ]]; then
        # Default: all chromosomes
        seq 1 22
    elif [[ "$chr_spec" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        # Range format: "21-22"
        seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
    else
        # List format: "21,22" or "21 22"
        echo "$chr_spec" | tr ',' ' ' | tr ' ' '\n' | sort -n | uniq
    fi
}

# =============================================================================
# TEMP DIRECTORY UTILITIES
# =============================================================================

# Create a temp directory.
# Prefer placing temp files under ${CFG_OUTDIR}/tmp (inside output mount) to avoid
# failures on compute nodes with small/full /tmp.
make_tmpdir() {
    local prefix="${1:-pgscalc}"
    local base=""

    if [[ -n "${CFG_OUTDIR:-}" ]]; then
        # If we're running a sumstat-specific pipeline, keep tmp contained within that sumstat folder.
        # This helps avoid collisions and keeps large temporary spill files (sort/plink/etc) co-located.
        if [[ -n "${CFG_SUMSTAT_NAME:-}" ]]; then
            base="${CFG_OUTDIR}/sumstats/${CFG_SUMSTAT_NAME}/tmp"
        else
            # For prep (and other non-sumstat commands), keep tmp contained within prep/.
            base="${CFG_OUTDIR}/prep/tmp"
        fi
        mkdir -p "$base" 2>/dev/null || true
        # If we can write to base, use it; else fall back to system mktemp
        if [[ -d "$base" ]] && [[ -w "$base" ]]; then
            mktemp -d "${base}/${prefix}.XXXXXX"
            return $?
        fi
    fi

    mktemp -d
}

# Validate chromosome number
is_valid_chromosome() {
    local chr="$1"
    
    if [[ "$chr" =~ ^[0-9]+$ ]] && [[ "$chr" -ge 1 ]] && [[ "$chr" -le 22 ]]; then
        return 0
    fi
    return 1
}

# =============================================================================
# STEP EXECUTION UTILITIES
# =============================================================================

# Run a step with logging and error handling
run_step() {
    local step_name="$1"
    local step_func="$2"
    shift 2
    
    log_step "Running step: $step_name"
    
    local start_time
    start_time=$(date +%s)
    
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log_info "DRY RUN: Would execute $step_func $*"
        return 0
    fi
    
    # Execute the step function
    if "$step_func" "$@"; then
        local end_time
        end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_info "Step '$step_name' completed in ${duration}s"
        return 0
    else
        log_error "Step '$step_name' failed"
        return 1
    fi
}

# =============================================================================
# VALIDATION UTILITIES
# =============================================================================

# Validate that required config variables are set
validate_required_config() {
    local prefix="${1:-CFG}"
    shift
    local required_vars=("$@")
    
    local missing=()
    for var in "${required_vars[@]}"; do
        local full_var="${prefix}_${var}"
        full_var="${full_var^^}"
        if [[ -z "${!full_var:-}" ]]; then
            missing+=("$var")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Missing required config variables: ${missing[*]}"
        return 1
    fi
    return 0
}

# =============================================================================
# GENOTYPE FILE UTILITIES
# =============================================================================

# Parse genotype manifest file and return file paths for a chromosome
get_geno_files_for_chr() {
    local genofile="$1"
    local genodir="$2"
    local chr="$3"
    local filetype="$4"  # pgen, pvar, psam, bed, bim, fam
    
    local result
    # Manifest format:
    #   <chr> <filetype> <filename>
    # Supports both TAB and SPACE separated files, and tolerates CRLF line endings.
    # Also accepts "chr1" in the chr column (normalized to "1").
    result=$(awk -v chr="$chr" -v ft="$filetype" -v gdir="$genodir" '
        BEGIN { FS = "[ \t]+" }
        {
            # Drop CR if file has Windows line endings
            sub(/\r$/, "", $0)
            c = $1; t = $2; f = $3
            sub(/^chr/, "", c)
            sub(/\r$/, "", t)
            sub(/\r$/, "", f)
            if (c == chr && t == ft) {
                print gdir "/" f
            }
        }
    ' "$genofile")
    
    echo "$result"
}

# Check if genotype format is plink2 or plink1
detect_geno_format() {
    local genofile="$1"
    
    if grep -q $'\tpgen\t' "$genofile" || grep -q $'\tpvar\t' "$genofile"; then
        echo "plink2"
    elif grep -q $'\tbed\t' "$genofile" || grep -q $'\tbim\t' "$genofile"; then
        echo "plink1"
    else
        echo "unknown"
    fi
}

# =============================================================================
# SBAYESR UTILITIES
# =============================================================================

# Build sbayesR command options from config
build_sbayesr_options() {
    local options=""
    
    # Add options from config
    [[ -n "${CFG_SBAYESR_GAMMA:-}" ]] && options+=" --gamma ${CFG_SBAYESR_GAMMA}"
    [[ -n "${CFG_SBAYESR_PI:-}" ]] && options+=" --pi ${CFG_SBAYESR_PI}"
    [[ -n "${CFG_SBAYESR_BURN_IN:-}" ]] && options+=" --burn-in ${CFG_SBAYESR_BURN_IN}"
    [[ -n "${CFG_SBAYESR_CHAIN_LENGTH:-}" ]] && options+=" --chain-length ${CFG_SBAYESR_CHAIN_LENGTH}"
    [[ -n "${CFG_SBAYESR_OUT_FREQ:-}" ]] && options+=" --out-freq ${CFG_SBAYESR_OUT_FREQ}"
    [[ -n "${CFG_SBAYESR_P_VALUE:-}" ]] && options+=" --p-value ${CFG_SBAYESR_P_VALUE}"
    [[ -n "${CFG_SBAYESR_RSQ:-}" ]] && options+=" --rsq ${CFG_SBAYESR_RSQ}"
    [[ -n "${CFG_SBAYESR_THREADS:-}" ]] && options+=" --thread ${CFG_SBAYESR_THREADS}"
    [[ -n "${CFG_SBAYESR_SEED:-}" ]] && options+=" --seed ${CFG_SBAYESR_SEED}"
    [[ -n "${CFG_SBAYESR_THIN:-}" ]] && options+=" --thin ${CFG_SBAYESR_THIN}"
    
    # Add flags
    [[ "${CFG_SBAYESR_EXCLUDE_MHC:-false}" == "true" ]] && options+=" --exclude-mhc"
    [[ "${CFG_SBAYESR_UNSCALE_GENOTYPE:-false}" == "true" ]] && options+=" --unscale-genotype"
    [[ "${CFG_SBAYESR_NO_MCMC_BIN:-false}" == "true" ]] && options+=" --no-mcmc-bin"
    [[ "${CFG_SBAYESR_IMPUTE_N:-false}" == "true" ]] && options+=" --impute-n"
    
    echo "$options"
}

# =============================================================================
# CLEANUP AND ERROR HANDLING
# =============================================================================

# Cleanup function for traps
cleanup() {
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        log_error "Script exited with code $exit_code"
    fi
    # Add any cleanup tasks here
}

# Set up error trap
setup_error_trap() {
    trap cleanup EXIT
}

# =============================================================================
# HELP UTILITIES
# =============================================================================

# Print a formatted help section
print_help_section() {
    local title="$1"
    shift
    local items=("$@")
    
    echo ""
    echo "$title:"
    for item in "${items[@]}"; do
        echo "  $item"
    done
}




