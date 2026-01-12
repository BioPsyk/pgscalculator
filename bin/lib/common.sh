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
    
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Skip comments and empty lines
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line// }" ]] && continue
        
        # Check for section (key with nested values)
        if [[ "$line" =~ ^([a-zA-Z_][a-zA-Z0-9_]*):$ ]]; then
            current_section="${BASH_REMATCH[1]}"
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
            
            if [[ -n "$current_section" ]]; then
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

# =============================================================================
# PATH AND DIRECTORY UTILITIES
# =============================================================================

# Get the output directory for a specific step
get_step_dir() {
    local outdir="$1"
    local step_name="$2"
    local sumstat_name="${3:-}"
    
    if [[ -n "$sumstat_name" ]]; then
        # Sumstat step directories live under intermediates/ (per-sumstat containment)
        echo "${outdir}/sumstats/${sumstat_name}/intermediates/${step_name}"
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

# Get the intermediates directory within a sumstat run directory
get_sumstat_intermediates_dir() {
    local sumstat_dir="$1"
    echo "${sumstat_dir}/intermediates"
}

# Move legacy step directories (placed directly under sumstat_dir) into
# sumstat_dir/intermediates/<step> for v2.1 output structure.
migrate_sumstat_step_dir() {
    local sumstat_dir="$1"
    local step="$2"

    local inter_dir
    inter_dir=$(get_sumstat_intermediates_dir "$sumstat_dir")
    ensure_dir "$inter_dir"

    local old_dir="${sumstat_dir}/${step}"
    local new_dir="${inter_dir}/${step}"

    if [[ -d "$old_dir" ]] && [[ ! -d "$new_dir" ]]; then
        mv "$old_dir" "$new_dir"
    fi
}

get_sumstat_step_dir() {
    local sumstat_dir="$1"
    local step="$2"
    local inter_dir
    inter_dir=$(get_sumstat_intermediates_dir "$sumstat_dir")
    echo "${inter_dir}/${step}"
}

# Migrate all known per-sumstat step directories into intermediates/.
# Safe to call repeatedly.
migrate_sumstat_all_step_dirs() {
    local sumstat_dir="$1"
    for step in formatted filtered posteriors posteriors_mapped scores scores_combined; do
        migrate_sumstat_step_dir "$sumstat_dir" "$step"
    done
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




