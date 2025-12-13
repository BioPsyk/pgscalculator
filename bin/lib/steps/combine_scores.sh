#!/bin/bash
# pgscalculator v2 - combine-scores step
# Combine per-chromosome scores into merged score file

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_combine_scores_deps() {
    require_command "awk" "awk is required for text processing"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_combine_scores() {
    local sumstat_name="$1"
    
    log_step "Running combine-scores for: $sumstat_name"
    
    # Check dependencies
    check_combine_scores_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local scores_dir="${sumstat_dir}/scores"
    local step_dir="${sumstat_dir}/scores_combined"
    ensure_dir "$step_dir"
    
    # Check that calc-score has been run
    require_dir "$scores_dir" "Run 'pgscalculator calc-score' first"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    # Find all score files
    local score_files=()
    for chr in $(get_chromosomes); do
        local score_file="${scores_dir}/chr${chr}.sscore"
        if [[ -f "$score_file" ]]; then
            score_files+=("$score_file")
        fi
    done
    
    if [[ ${#score_files[@]} -eq 0 ]]; then
        log_error "No score files found in ${scores_dir}"
        exit 1
    fi
    
    log_info "Found ${#score_files[@]} chromosome score files"
    
    # Step 1: Build IID reference from first score file
    log_substep "Building IID reference"
    local ref_file="${step_dir}/iid_ref.txt"
    build_iid_reference "${score_files[0]}" "$ref_file"
    
    # Step 2: Combine scores across chromosomes
    log_substep "Combining chromosome scores"
    local merged_file="${step_dir}/merged.sscore"
    combine_chromosome_scores "$ref_file" "${scores_dir}" "$merged_file"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local sample_count
    sample_count=$(wc -l < "$merged_file")
    sample_count=$((sample_count - 1))
    log_info "Combined scores for ${sample_count} samples"
    log_info "Output: ${merged_file}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

build_iid_reference() {
    local score_file="$1"
    local ref_file="$2"
    
    # Extract IID column and create reference
    # Score file format: FID IID ALLELE_CT NAMED_ALLELE_DOSAGE_SUM SCORE1_SUM
    
    # Find IID column
    local header
    header=$(head -1 "$score_file")
    
    local iid_col
    iid_col=$(echo "$header" | awk -F'\t' '{
        for(i=1; i<=NF; i++) {
            if($i == "IID" || $i == "#IID") {
                print i
                exit
            }
        }
    }')
    
    if [[ -z "$iid_col" ]]; then
        log_error "Could not find IID column in score file"
        exit 1
    fi
    
    # Extract unique IIDs in order
    awk -F'\t' -v iid_col="$iid_col" '
        NR > 1 {
            iid = $iid_col
            if (!(iid in seen)) {
                print iid
                seen[iid] = 1
            }
        }
    ' "$score_file" > "$ref_file"
    
    local count
    count=$(wc -l < "$ref_file")
    log_debug "Reference: ${count} unique IIDs"
}

combine_chromosome_scores() {
    local ref_file="$1"
    local scores_dir="$2"
    local output_file="$3"
    
    local n_samples
    n_samples=$(wc -l < "$ref_file")
    log_debug "Processing ${n_samples} samples"
    
    # Build explicit list of score files (avoid glob issues)
    local score_file_list=""
    for chr in $(get_chromosomes); do
        local score_file="${scores_dir}/chr${chr}.sscore"
        if [[ -f "$score_file" ]]; then
            score_file_list="${score_file_list} ${score_file}"
            log_debug "Adding chr${chr} scores"
        fi
    done
    
    # Combine all chromosome scores
    log_substep "Aggregating scores across chromosomes"
    
    # Simple and robust: use awk to sum scores from all files
    # Score file format: IID<tab>SCORE1_SUM
    awk -F'\t' '
        FNR == 1 { next }  # Skip headers
        {
            iid = $1
            score = $2
            scores[iid] += score
            if (!(iid in order)) {
                order[iid] = ++n
                iids[n] = iid
            }
        }
        END {
            print "IID\tSCORE1_SUM"
            for (i=1; i<=n; i++) {
                iid = iids[i]
                print iid "\t" scores[iid]
            }
        }
    ' $score_file_list > "$output_file"
}




