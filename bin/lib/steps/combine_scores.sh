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
    
    # Create temporary directory for intermediate files
    local tmpdir
    tmpdir=$(mktemp -d)
    
    # Load IID reference into array
    local -A iid_to_idx
    local -a iids
    local idx=0
    while IFS= read -r iid; do
        iids+=("$iid")
        iid_to_idx["$iid"]=$idx
        ((idx++))
    done < "$ref_file"
    
    local n_samples=${#iids[@]}
    log_debug "Processing ${n_samples} samples"
    
    # Initialize score sums for each sample
    declare -a score_sums
    declare -a allele_counts
    for ((i=0; i<n_samples; i++)); do
        score_sums[$i]=0
        allele_counts[$i]=0
    done
    
    # Process each chromosome score file
    for chr in $(get_chromosomes); do
        local score_file="${scores_dir}/chr${chr}.sscore"
        if [[ ! -f "$score_file" ]]; then
            continue
        fi
        
        log_debug "Adding chr${chr} scores"
        
        # Get column indices
        local header
        header=$(head -1 "$score_file")
        
        local iid_col score_col allele_col
        iid_col=$(echo "$header" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="IID"||$i=="#IID") print i}')
        score_col=$(echo "$header" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="SCORE1_SUM") print i}')
        allele_col=$(echo "$header" | awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="ALLELE_CT") print i}')
        
        # Sum scores using awk for efficiency
        awk -F'\t' -v iid_col="$iid_col" -v score_col="$score_col" -v allele_col="$allele_col" '
            BEGIN {
                while ((getline line < "'"$ref_file"'") > 0) {
                    idx++
                    ref_iid[line] = idx
                }
            }
            NR > 1 {
                iid = $iid_col
                if (iid in ref_iid) {
                    idx = ref_iid[iid]
                    scores[idx] += $score_col
                    alleles[idx] += $allele_col
                }
            }
            END {
                for (idx=1; idx<=length(ref_iid); idx++) {
                    print scores[idx], alleles[idx]
                }
            }
        ' "$score_file" > "${tmpdir}/chr${chr}_partial.txt"
        
    done
    
    # Combine all partial sums
    log_substep "Aggregating scores across chromosomes"
    
    # Write output header
    echo -e "IID\tALLELE_CT\tSCORE1_SUM" > "$output_file"
    
    # Use paste to combine all partial files, then sum
    paste "${tmpdir}"/chr*_partial.txt 2>/dev/null | \
    awk -v n_samples="$n_samples" '
        BEGIN {
            # Read IIDs
            idx = 0
            while ((getline line < "'"$ref_file"'") > 0) {
                idx++
                iids[idx] = line
            }
        }
        {
            total_score = 0
            total_alleles = 0
            for (i=1; i<=NF; i+=2) {
                total_score += $i
                total_alleles += $(i+1)
            }
            print iids[NR], total_alleles, total_score
        }
    ' OFS='\t' >> "$output_file" || true
    
    # Alternative simpler approach using awk to process all files
    awk -F'\t' '
        ARGIND == 1 {
            # Load IID reference
            ref_iids[FNR] = $0
            n_iids = FNR
            next
        }
        FNR == 1 {
            # Get column indices from header
            for(i=1; i<=NF; i++) {
                if($i == "IID" || $i == "#IID") iid_col = i
                if($i == "SCORE1_SUM") score_col = i
                if($i == "ALLELE_CT") allele_col = i
            }
            next
        }
        {
            iid = $iid_col
            scores[iid] += $score_col
            alleles[iid] += $allele_col
        }
        END {
            for (i=1; i<=n_iids; i++) {
                iid = ref_iids[i]
                print iid, alleles[iid]+0, scores[iid]+0
            }
        }
    ' OFS='\t' "$ref_file" "${scores_dir}"/chr*.sscore > "${tmpdir}/combined.txt"
    
    # Write final output with header
    echo -e "IID\tALLELE_CT\tSCORE1_SUM" > "$output_file"
    cat "${tmpdir}/combined.txt" >> "$output_file"
    
    # Clean up
    rm -rf "$tmpdir"
}




