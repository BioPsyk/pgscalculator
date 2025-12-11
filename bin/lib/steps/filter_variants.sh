#!/bin/bash
# pgscalculator v2 - filter-variants step
# Filter sumstat to whitelist variants (early variant reduction for memory optimization)

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_filter_variants_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "sort" "sort is required for sorting"
    require_command "join" "join is required for file merging"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    
    # Check that prep-whitelist has been run
    require_file "${prep_dir}/whitelist/variant_whitelist.tsv" "Run 'pgscalculator prep-whitelist' first"
    require_file "${prep_dir}/whitelist/whitelist_rsids_sorted" "Run 'pgscalculator prep-whitelist' first"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_filter_variants() {
    local sumstat_name="$1"
    
    log_step "Running filter-variants for: $sumstat_name"
    
    # Check dependencies
    check_filter_variants_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local format_dir="${sumstat_dir}/formatted"
    local step_dir="${sumstat_dir}/filtered"
    ensure_dir "$step_dir"
    
    # Check that format-sumstat has been run
    require_file "${format_dir}/sumstat_formatted.tsv.gz" "Run 'pgscalculator format-sumstat' first"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local whitelist_file="${prep_dir}/whitelist/whitelist_rsids_sorted"
    local formatted_sumstat="${format_dir}/sumstat_formatted.tsv.gz"
    
    # Count input variants
    local input_count
    input_count=$(zcat "$formatted_sumstat" | wc -l)
    input_count=$((input_count - 1))
    log_info "Input variants: ${input_count}"
    
    # Step 1: Filter sumstat to whitelist variants
    log_substep "Filtering to whitelist variants"
    filter_to_whitelist "$formatted_sumstat" "$whitelist_file" "${step_dir}/sumstat_filtered.tsv"
    
    # Step 2: Split filtered sumstat by chromosome
    log_substep "Splitting filtered sumstat by chromosome"
    split_filtered_by_chr "${step_dir}/sumstat_filtered.tsv" "$step_dir"
    
    # Step 3: Compress the main filtered file
    gzip -f "${step_dir}/sumstat_filtered.tsv"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local output_count
    output_count=$(zcat "${step_dir}/sumstat_filtered.tsv.gz" | wc -l)
    output_count=$((output_count - 1))
    local reduction_pct
    reduction_pct=$(awk "BEGIN {printf \"%.1f\", (1 - $output_count / $input_count) * 100}")
    
    log_info "Output variants: ${output_count} (${reduction_pct}% reduction)"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

filter_to_whitelist() {
    local input_sumstat="$1"
    local whitelist_file="$2"
    local output_file="$3"
    
    # Get header and find SNP/RSID column
    local header
    header=$(zcat "$input_sumstat" | head -1)
    
    local snp_col
    snp_col=$(echo "$header" | awk -F'\t' '{
        for(i=1; i<=NF; i++) {
            if($i == "SNP" || $i == "RSID" || $i == "rsid" || $i == "ID") {
                print i
                exit
            }
        }
    }')
    
    if [[ -z "$snp_col" ]]; then
        log_error "Could not find SNP/RSID column in sumstat"
        exit 1
    fi
    
    log_debug "SNP column index: $snp_col"
    
    # Create temporary sorted files for join
    local tmpdir
    tmpdir=$(mktemp -d)
    
    # Sort sumstat by SNP column
    zcat "$input_sumstat" | tail -n +2 | \
        sort -t$'\t' -k${snp_col},${snp_col} > "${tmpdir}/sumstat_sorted.tsv"
    
    # Join with whitelist (whitelist is already sorted)
    # Keep all fields from sumstat where SNP matches whitelist
    join -t$'\t' -1 1 -2 ${snp_col} -o 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20 \
        "$whitelist_file" "${tmpdir}/sumstat_sorted.tsv" 2>/dev/null | \
        sed 's/\t$//' | grep -v '^\s*$' > "${tmpdir}/filtered_body.tsv" || true
    
    # Alternative approach using awk for more flexibility
    awk -F'\t' -v OFS='\t' -v snp_col="$snp_col" '
        ARGIND == 1 {
            whitelist[$1] = 1
            next
        }
        {
            if ($snp_col in whitelist) {
                print
            }
        }
    ' "$whitelist_file" <(zcat "$input_sumstat" | tail -n +2) > "${tmpdir}/filtered_body2.tsv"
    
    # Write output with header
    echo "$header" > "$output_file"
    cat "${tmpdir}/filtered_body2.tsv" >> "$output_file"
    
    # Clean up
    rm -rf "$tmpdir"
}

split_filtered_by_chr() {
    local input_file="$1"
    local step_dir="$2"
    
    # Get header
    local header
    header=$(head -1 "$input_file")
    
    # Find CHR column index
    local chr_col
    chr_col=$(echo "$header" | awk -F'\t' '{
        for(i=1; i<=NF; i++) {
            if($i == "CHR" || $i == "chr" || $i == "#CHR") {
                print i
                exit
            }
        }
    }')
    
    if [[ -z "$chr_col" ]]; then
        log_error "Could not find CHR column in filtered sumstat"
        exit 1
    fi
    
    # Split by chromosome
    awk -F'\t' -v OFS='\t' -v chr_col="$chr_col" -v outdir="$step_dir" -v header="$header" '
        NR == 1 { next }
        {
            chr = $chr_col
            if (chr >= 1 && chr <= 22) {
                outfile = outdir "/chr" chr "_filtered.tsv"
                if (!(chr in seen)) {
                    print header > outfile
                    seen[chr] = 1
                }
                print >> outfile
            }
        }
    ' "$input_file"
    
    # Report per-chromosome counts
    for chr in $(get_chromosomes); do
        local chr_file="${step_dir}/chr${chr}_filtered.tsv"
        if [[ -f "$chr_file" ]]; then
            local count
            count=$(wc -l < "$chr_file")
            count=$((count - 1))
            log_debug "chr${chr}: ${count} variants"
        fi
    done
}




