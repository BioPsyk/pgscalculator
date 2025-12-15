#!/bin/bash
# pgscalculator v2 - format-sumstat step
# Format summary statistics: add GRCh37 build coordinates only
# (N/EAF/B/SE derivation happens AFTER filtering in filter-variants step)

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_format_sumstat_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "zcat" "zcat is required for reading gzipped files"
    require_command "gzip" "gzip is required for compressing output"
    
    # Check config variables
    validate_required_config "CFG" "INPUT" "OUTDIR"
    
    local input_dir="${CFG_INPUT}"
    
    # Check for cleaned sumstat files from cleansumstats
    require_file "${input_dir}/cleaned_GRCh38.gz" "Cleaned sumstat (GRCh38) not found in input directory"
    require_file "${input_dir}/cleaned_GRCh37.gz" "Cleaned sumstat (GRCh37 map) not found in input directory"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_format_sumstat() {
    local sumstat_name="$1"
    
    log_step "Running format-sumstat for: $sumstat_name"
    
    # Check dependencies
    check_format_sumstat_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "formatted"
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    ensure_dir "$step_dir"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local input_dir="${CFG_INPUT}"
    local input_grch38="${input_dir}/cleaned_GRCh38.gz"
    local input_grch37="${input_dir}/cleaned_GRCh37.gz"
    
    # Step 1: Add GRCh37 coordinates to sumstat
    log_substep "Adding GRCh37 build coordinates"
    add_build_coordinates "$input_grch38" "$input_grch37" "${step_dir}/sumstat_with_b37.tsv.gz"
    
    # Step 2: Filter NA coordinates (remove variants without valid GRCh37 positions)
    log_substep "Filtering NA coordinates"
    filter_na_coordinates_gz "${step_dir}/sumstat_with_b37.tsv.gz" "${step_dir}/sumstat_formatted.tsv.gz"
    
    # Clean up intermediate file
    rm -f "${step_dir}/sumstat_with_b37.tsv.gz"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local final_count
    final_count=$(zcat "${step_dir}/sumstat_formatted.tsv.gz" | wc -l)
    final_count=$((final_count - 1))
    log_info "Formatted sumstat with GRCh37 coordinates: ${final_count} variants"
    log_info "Output directory: ${step_dir}"
    log_info "Note: N/EAF/B/SE derivation will happen after filtering"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

add_build_coordinates() {
    local input_grch38="$1"
    local input_grch37="$2"
    local output_file="$3"
    
    # GRCh37 file has: CHR, POS, RSID (3 columns)
    # GRCh38 file has: CHR, POS, 0, RSID, EffectAllele, ...
    # After paste (b37 first, then b38):
    #   1: CHR_b37, 2: POS_b37, 3: RSID_b37, 4: CHR_b38, 5: POS_b38, 6+: rest
    # We want: CHR_b37, POS_b37, POS_b38, 0, RSID, ... (drop RSID_b37 and CHR_b38)
    # Use cut -f1-2,5- to skip columns 3 and 4
    paste <(zcat "$input_grch37") <(zcat "$input_grch38") | \
        awk -F'\t' -v OFS='\t' '{
            # Replace empty values with NA
            for(i=1; i<=NF; i++) if($i=="") $i="NA"
            print
        }' | cut -f1-2,5- | gzip -c > "$output_file"
}

filter_na_coordinates_gz() {
    local input="$1"
    local output="$2"
    
    # Filter out rows where CHR or POS (for b37) is NA
    zcat "$input" | awk -F'\t' -v OFS='\t' '
        NR == 1 {
            # Find CHR and POS column indices
            for(i=1; i<=NF; i++) {
                if($i == "CHR") chr_col = i
                if($i == "POS") pos_col = i
            }
            print
            next
        }
        {
            if ($chr_col != "NA" && $chr_col != "" && 
                $pos_col != "NA" && $pos_col != "") {
                print
            }
        }
    ' | gzip -c > "$output"
}
