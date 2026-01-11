#!/bin/bash
# pgscalculator v2 - finalize-output step
# Generate final output files: sumstat_augmented.tsv.gz and copy config to details/

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_finalize_output_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "gzip" "gzip is required for compression"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_finalize_output() {
    local sumstat_name="$1"
    
    log_step "Running finalize-output for: $sumstat_name"
    
    # Check dependencies
    check_finalize_output_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local step_dir="${sumstat_dir}/details"
    ensure_dir "$step_dir"
    
    # Check prerequisites
    require_file "${sumstat_dir}/scores.tsv.gz" "Run 'pgscalculator combine-scores' first"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    # Step 1: Combine all posteriors into single file
    log_substep "Combining posteriors from all chromosomes"
    local posteriors_combined="${sumstat_dir}/posteriors_combined.tsv"
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors_mapped"
    local posteriors_mapped_dir
    posteriors_mapped_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    combine_posteriors "$posteriors_mapped_dir" "$posteriors_combined"
    
    # Step 2: Generate sumstat_augmented.tsv.gz
    log_substep "Generating augmented sumstat"
    generate_augmented_sumstat "$sumstat_dir" "$prep_dir" "$posteriors_combined"
    
    # Step 3: Copy variant map to sumstat root (v1-compatible artifact)
    log_substep "Writing variant_map.tsv.gz"
    write_variant_map "$prep_dir" "$sumstat_dir"

    # Step 4: Copy config to details/
    log_substep "Copying configuration to details/"
    copy_config_to_details "$outdir" "$step_dir"
    
    # Step 5: Generate run summary
    log_substep "Generating run summary"
    generate_run_summary "$sumstat_dir" "$step_dir"

    # Step 6: Generate stepwise details TSVs
    log_substep "Generating stepwise details"
    generate_stepwise_details "$sumstat_dir" "$step_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    log_info "Finalization complete"
    log_info "Output directory: ${sumstat_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

write_variant_map() {
    local prep_dir="$1"
    local sumstat_dir="$2"

    local variant_map_src="${prep_dir}/variant_map.tsv"
    local variant_map_out="${sumstat_dir}/variant_map.tsv.gz"

    if [[ ! -f "$variant_map_src" ]]; then
        log_warn "variant_map.tsv not found at: ${variant_map_src} (skipping)"
        return 0
    fi

    gzip -c "$variant_map_src" > "$variant_map_out"
    log_debug "Wrote variant map: ${variant_map_out}"
}

combine_posteriors() {
    local posteriors_dir="$1"
    local output_file="$2"
    
    # Header for combined posteriors
    echo -e "RSID\tGENO_ID\tA1\tA2\tFREQ\tEFFECT\tSE\tPIP" > "$output_file"
    
    local total_variants=0
    
    for chr in $(get_chromosomes); do
        local posterior_file="${posteriors_dir}/chr${chr}.snpRes"
        if [[ -f "$posterior_file" ]]; then
            # Append without header
            tail -n +2 "$posterior_file" >> "$output_file"
            local chr_count
            chr_count=$(tail -n +2 "$posterior_file" | wc -l)
            total_variants=$((total_variants + chr_count))
        fi
    done
    
    log_info "Combined ${total_variants} posterior variants"
}

generate_augmented_sumstat() {
    local sumstat_dir="$1"
    local prep_dir="$2"
    local posteriors_file="$3"
    
    migrate_sumstat_step_dir "$sumstat_dir" "formatted"
    local formatted_sumstat
    formatted_sumstat="$(get_sumstat_step_dir "$sumstat_dir" "formatted")/sumstat_formatted.tsv.gz"
    local variant_map="${prep_dir}/variant_map.tsv"
    local output_file="${sumstat_dir}/sumstat_augmented.tsv.gz"
    
    if [[ ! -f "$formatted_sumstat" ]]; then
        log_warn "Formatted sumstat not found, skipping augmented sumstat generation"
        return 0
    fi
    
    # Join formatted sumstat with posteriors and variant map
    # Output: all columns from formatted sumstat + GENO_ID + POST_EFFECT + POST_PIP + IN_ANALYSIS
    
    local tmpdir
    tmpdir=$(make_tmpdir "finalize_output")
    
    # Load posteriors into lookup (RSID -> EFFECT, PIP)
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            rsid = $1
            effect = $6
            pip = $8
            post_effect[rsid] = effect
            post_pip[rsid] = pip
        }
        END {
            for (rsid in post_effect) {
                print rsid, post_effect[rsid], post_pip[rsid]
            }
        }
    ' "$posteriors_file" > "${tmpdir}/posteriors_lookup.tsv"
    
    # Load variant map (ld_rsid -> pvar_snpid)
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            # chrpos, pvar_a1, pvar_a2, pvar_snpid, ld_a1, ld_a2, ld_rsid
            print $7, $4
        }
    ' "$variant_map" > "${tmpdir}/varmap_lookup.tsv"
    
    # Join with formatted sumstat
    zcat "$formatted_sumstat" | awk -F'\t' -v OFS='\t' '
        # Load posteriors lookup
        ARGIND == 1 {
            post_effect[$1] = $2
            post_pip[$1] = $3
            next
        }
        # Load variant map lookup
        ARGIND == 2 {
            geno_id[$1] = $2
            next
        }
        # Process sumstat
        ARGIND == 3 {
            if (FNR == 1) {
                # Find SNP/RSID column
                for (i=1; i<=NF; i++) {
                    if ($i == "SNP" || $i == "RSID" || $i == "rsid" || $i == "ID") snp_col = i
                }
                print $0, "GENO_ID", "POST_EFFECT", "POST_PIP", "IN_ANALYSIS"
                next
            }
            
            rsid = $snp_col
            gid = (rsid in geno_id) ? geno_id[rsid] : "NA"
            pe = (rsid in post_effect) ? post_effect[rsid] : "NA"
            pp = (rsid in post_pip) ? post_pip[rsid] : "NA"
            in_analysis = (pe != "NA") ? "Y" : "N"
            
            print $0, gid, pe, pp, in_analysis
        }
    ' "${tmpdir}/posteriors_lookup.tsv" "${tmpdir}/varmap_lookup.tsv" - | gzip > "$output_file"
    
    # Clean up
    rm -rf "$tmpdir"
    
    local variant_count
    variant_count=$(zcat "$output_file" | wc -l)
    variant_count=$((variant_count - 1))
    log_info "Generated augmented sumstat with ${variant_count} variants"
}

copy_config_to_details() {
    local outdir="$1"
    local details_dir="$2"
    
    # Copy config.yaml to details/
    local config_file="${outdir}/config.yaml"
    if [[ -f "$config_file" ]]; then
        cp "$config_file" "${details_dir}/config.yaml"
        log_debug "Copied config.yaml to details/"
    fi
}

generate_run_summary() {
    local sumstat_dir="$1"
    local details_dir="$2"
    
    local summary_file="${details_dir}/run_summary.txt"
    
    {
        echo "pgscalculator v2 Run Summary"
        echo "============================"
        echo ""
        echo "Run completed: $(date)"
        echo ""
        echo "Output files:"
        
        if [[ -f "${sumstat_dir}/scores.tsv.gz" ]]; then
            local score_count
            score_count=$(zcat "${sumstat_dir}/scores.tsv.gz" | wc -l)
            score_count=$((score_count - 1))
            echo "  - scores.tsv.gz: ${score_count} samples"
        fi
        
        if [[ -f "${sumstat_dir}/sumstat_augmented.tsv.gz" ]]; then
            local var_count
            var_count=$(zcat "${sumstat_dir}/sumstat_augmented.tsv.gz" | wc -l)
            var_count=$((var_count - 1))
            echo "  - sumstat_augmented.tsv.gz: ${var_count} variants"
        fi
        
        if [[ -f "${sumstat_dir}/posteriors_combined.tsv" ]]; then
            local post_count
            post_count=$(wc -l < "${sumstat_dir}/posteriors_combined.tsv")
            post_count=$((post_count - 1))
            echo "  - posteriors_combined.tsv: ${post_count} posteriors"
        fi
        
        echo ""
        echo "Chromosomes processed:"
        for chr in $(get_chromosomes); do
            migrate_sumstat_step_dir "$sumstat_dir" "scores"
            if [[ -f "$(get_sumstat_step_dir "$sumstat_dir" "scores")/chr${chr}.sscore" ]]; then
                echo "  - chr${chr}: OK"
            fi
        done
        
    } > "$summary_file"
    
    log_debug "Generated run summary"
}

generate_stepwise_details() {
    local sumstat_dir="$1"
    local details_dir="$2"

    local steps_file="${details_dir}/steps.tsv"
    local score_file="${details_dir}/score_steps.tsv"

    migrate_sumstat_all_step_dirs "$sumstat_dir"
    local formatted_dir filtered_dir post_dir mapped_dir scores_dir
    formatted_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    filtered_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    post_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors")
    mapped_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    scores_dir=$(get_sumstat_step_dir "$sumstat_dir" "scores")

    # Helpers (avoid hard failure if files missing/broken; report 0)
    local n_formatted n_filtered n_post n_mapped n_scores_files n_samples n_score_variants
    n_formatted=0
    n_filtered=0
    n_post=0
    n_mapped=0

    if [[ -f "${formatted_dir}/sumstat_formatted.tsv.gz" ]]; then
        n_formatted=$(gzip -cd "${formatted_dir}/sumstat_formatted.tsv.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_formatted" -gt 0 ]]; then n_formatted=$((n_formatted - 1)); else n_formatted=0; fi
    fi
    if [[ -f "${filtered_dir}/sumstat_filtered.tsv.gz" ]]; then
        n_filtered=$(gzip -cd "${filtered_dir}/sumstat_filtered.tsv.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_filtered" -gt 0 ]]; then n_filtered=$((n_filtered - 1)); else n_filtered=0; fi
    fi
    if compgen -G "${post_dir}/chr*.snpRes" >/dev/null 2>&1; then
        # sum across chr files: (lines - 1)
        n_post=$(for f in "${post_dir}"/chr*.snpRes; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi
    if compgen -G "${mapped_dir}/chr*.snpRes" >/dev/null 2>&1; then
        n_mapped=$(for f in "${mapped_dir}"/chr*.snpRes; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi

    # Failure markers
    local post_fail mapped_fail score_fail
    post_fail=0
    mapped_fail=0
    score_fail=0
    if compgen -G "${post_dir}/work_chr*/FAILED" >/dev/null 2>&1; then
        post_fail=$(ls -1 "${post_dir}"/work_chr*/FAILED 2>/dev/null | wc -l | awk '{print $1}')
    fi
    if compgen -G "${mapped_dir}/FAILED_chr*" >/dev/null 2>&1; then
        mapped_fail=$(ls -1 "${mapped_dir}"/FAILED_chr* 2>/dev/null | wc -l | awk '{print $1}')
    fi
    if compgen -G "${scores_dir}/FAILED_chr*" >/dev/null 2>&1; then
        score_fail=$(ls -1 "${scores_dir}"/FAILED_chr* 2>/dev/null | wc -l | awk '{print $1}')
    fi

    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "format-sumstat\t${n_formatted}\t${n_formatted}\tformatted sumstat rows"
        echo -e "filter-variants\t${n_formatted}\t${n_filtered}\tfiltered sumstat rows"
        echo -e "calc-posteriors\t${n_filtered}\t${n_post}\tposterior rows (failed_chr=${post_fail})"
        echo -e "format-posteriors\t${n_post}\t${n_mapped}\tmapped posterior rows (failed_chr=${mapped_fail})"
    } > "$steps_file"

    # Score summary
    n_scores_files=0
    n_samples=0
    n_score_variants=0
    if compgen -G "${scores_dir}/chr*.sscore" >/dev/null 2>&1; then
        n_scores_files=$(ls -1 "${scores_dir}"/chr*.sscore 2>/dev/null | wc -l | awk '{print $1}')
        # total variants scored: prefer work_chr*/variants.txt if present
        if compgen -G "${scores_dir}/work_chr*/variants.txt" >/dev/null 2>&1; then
            n_score_variants=$(for f in "${scores_dir}"/work_chr*/variants.txt; do wc -l < "$f"; done | awk '{s+=$1} END{print s+0}')
        fi
    fi
    if [[ -f "${sumstat_dir}/scores.tsv.gz" ]]; then
        n_samples=$(gzip -cd "${sumstat_dir}/scores.tsv.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_samples" -gt 0 ]]; then n_samples=$((n_samples - 1)); else n_samples=0; fi
    fi

    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "calc-score\t${n_mapped}\t${n_score_variants}\tper-chr scoring ran (score_files=${n_scores_files}, failed_chr=${score_fail})"
        echo -e "combine-scores\t${n_score_variants}\t${n_score_variants}\tcombined score table written"
        echo -e "finalize-output\t${n_score_variants}\t${n_score_variants}\tfinal outputs written (samples=${n_samples})"
    } > "$score_file"

    log_debug "Wrote details: ${steps_file}, ${score_file}"
}
