#!/bin/bash
# pgscalculator v2 - format-posteriors step
# Map posteriors to genotype variant IDs for scoring

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_format_posteriors_deps() {
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
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_format_posteriors() {
    local sumstat_name="$1"
    
    log_step "Running format-posteriors for: $sumstat_name"
    
    # Check dependencies
    check_format_posteriors_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local posteriors_dir="${sumstat_dir}/posteriors"
    local step_dir="${sumstat_dir}/posteriors_mapped"
    ensure_dir "$step_dir"
    
    # Check that calc-posteriors has been run
    require_dir "$posteriors_dir" "Run 'pgscalculator calc-posteriors' first"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local whitelist_file="${prep_dir}/whitelist/variant_whitelist.tsv"
    
    # Build RSID to genotype ID mapping
    log_substep "Building RSID to genotype ID mapping"
    local rsid_map="${step_dir}/rsid_to_genoid.tsv"
    create_rsid_mapping "$whitelist_file" "$rsid_map"
    
    # Process each chromosome
    log_substep "Mapping posteriors to genotype IDs"
    local total_mapped=0
    
    for chr in $(get_chromosomes); do
        local posterior_file="${posteriors_dir}/chr${chr}.snpRes"
        
        if [[ ! -f "$posterior_file" ]]; then
            log_debug "No posteriors for chr${chr}, skipping"
            continue
        fi
        
        local output_file="${step_dir}/chr${chr}.snpRes"
        local mapped_count
        mapped_count=$(map_posteriors_for_chr "$chr" "$posterior_file" "$rsid_map" "$output_file")
        
        total_mapped=$((total_mapped + mapped_count))
        log_debug "chr${chr}: ${mapped_count} variants mapped"
    done
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    log_info "Total variants mapped: ${total_mapped}"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

create_rsid_mapping() {
    local whitelist_file="$1"
    local output_file="$2"
    
    # Extract RSID -> pvar_snpid mapping from whitelist
    # whitelist format: ld_rsid, pvar_snpid, chrpos
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            print $1, $2  # rsid, genotype_id
        }
    ' "$whitelist_file" | LC_ALL=C sort -k1,1 > "$output_file"
    
    local count
    count=$(wc -l < "$output_file")
    log_debug "Created mapping with ${count} variants"
}

map_posteriors_for_chr() {
    local chr="$1"
    local posterior_file="$2"
    local rsid_map="$3"
    local output_file="$4"
    
    # sbayesR .snpRes format: SNP A1 A2 b se pval Freq N effect pj
    # We need to:
    # 1. Map SNP (RSID) to genotype variant ID
    # 2. Keep allele information for scoring
    
    local tmpdir
    tmpdir=$(mktemp -d)
    
    # Get header from posterior file
    local header
    header=$(head -1 "$posterior_file")
    
    # Sort posteriors by SNP column (column 1)
    tail -n +2 "$posterior_file" | LC_ALL=C sort -k1,1 > "${tmpdir}/posteriors_sorted.tsv"
    
    # Join with RSID mapping
    # Output: genotype_id + all posterior columns
    LC_ALL=C join -t' ' -1 1 -2 1 \
        -o 2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 \
        "${tmpdir}/posteriors_sorted.tsv" "$rsid_map" 2>/dev/null > "${tmpdir}/mapped.tsv" || true
    
    # Alternative: use awk for more robust joining
    awk -F' ' -v OFS='\t' '
        ARGIND == 1 {
            # Load RSID mapping
            rsid_to_geno[$1] = $2
            next
        }
        NR == 1 { next }  # Skip header in posteriors
        {
            rsid = $1
            if (rsid in rsid_to_geno) {
                geno_id = rsid_to_geno[rsid]
                # Output: genotype_id, A1, A2, b, se, pval, Freq, N, effect, pj
                print geno_id, $2, $3, $4, $5, $6, $7, $8, $9, $10
            }
        }
    ' "$rsid_map" "$posterior_file" > "${tmpdir}/mapped_awk.tsv"
    
    # Write output with modified header
    # Header: ID A1 A2 b se pval Freq N effect pj
    echo -e "ID\tA1\tA2\tb\tse\tpval\tFreq\tN\teffect\tpj" > "$output_file"
    cat "${tmpdir}/mapped_awk.tsv" >> "$output_file"
    
    # Count mapped variants
    local count
    count=$(wc -l < "${tmpdir}/mapped_awk.tsv")
    
    # Clean up
    rm -rf "$tmpdir"
    
    echo "$count"
}




