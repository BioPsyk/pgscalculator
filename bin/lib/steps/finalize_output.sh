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
    combine_posteriors "${sumstat_dir}/posteriors_mapped" "$posteriors_combined"
    
    # Step 2: Generate sumstat_augmented.tsv.gz
    log_substep "Generating augmented sumstat"
    generate_augmented_sumstat "$sumstat_dir" "$prep_dir" "$posteriors_combined"
    
    # Step 3: Copy config to details/
    log_substep "Copying configuration to details/"
    copy_config_to_details "$outdir" "$step_dir"
    
    # Step 4: Generate run summary
    log_substep "Generating run summary"
    generate_run_summary "$sumstat_dir" "$step_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    log_info "Finalization complete"
    log_info "Output directory: ${sumstat_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

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
    
    local formatted_sumstat="${sumstat_dir}/formatted/sumstat_formatted.tsv.gz"
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
            if [[ -f "${sumstat_dir}/scores/chr${chr}.sscore" ]]; then
                echo "  - chr${chr}: OK"
            fi
        done
        
    } > "$summary_file"
    
    log_debug "Generated run summary"
}

