#!/bin/bash
# pgscalculator v2 - prep-inclusion-list step
# Create variant inclusion list by intersecting genotypes with LD reference and applying filters

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_prep_inclusion_list_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "sort" "sort is required for sorting"
    require_command "join" "join is required for file merging"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    
    # Check that prep-genotypes has been run
    require_dir "${prep_dir}/genotypes" "Run 'pgscalculator prep-genotypes' first"
    require_file "${prep_dir}/genotypes/snplist_sorted" "Run 'pgscalculator prep-genotypes' first"
    
    # Check that prep-ldref has been run  
    require_dir "${prep_dir}/ldref" "Run 'pgscalculator prep-ldref' first"
    require_file "${prep_dir}/ldref/ld_rsids_all" "Run 'pgscalculator prep-ldref' first"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_prep_inclusion_list() {
    log_step "Running prep-inclusion-list"
    
    # Check dependencies
    check_prep_inclusion_list_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local step_dir
    step_dir=$(get_step_dir "$outdir" "inclusion_list")
    ensure_dir "$step_dir"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    # Get filter thresholds from config (with defaults)
    local info_threshold="${CFG_INFO_THRESHOLD:-0.8}"
    local maf_threshold="${CFG_MAF_THRESHOLD:-0.01}"
    
    log_info "INFO threshold: ${info_threshold}"
    log_info "MAF threshold: ${maf_threshold}"
    
    local geno_dir="${prep_dir}/genotypes"
    local ldref_dir="${prep_dir}/ldref"
    
    # Step 1: Create per-chromosome variant maps
    log_substep "Creating per-chromosome variant maps"
    
    for chr in $(get_chromosomes); do
        create_chr_variant_map "$chr" "$geno_dir" "$ldref_dir" "$step_dir"
    done
    
    # Step 2: Combine all chromosome maps (output to prep/ level)
    log_substep "Combining chromosome variant maps"
    combine_variant_maps "$step_dir" "$prep_dir"
    
    # Step 3: Create final inclusion list (variants present in both genotypes and LD ref)
    log_substep "Creating final variant inclusion list"
    create_final_inclusion_list "$step_dir" "$prep_dir" "$info_threshold" "$maf_threshold"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local inclusion_count
    inclusion_count=$(wc -l < "${step_dir}/variant_inclusion_list.tsv")
    inclusion_count=$((inclusion_count - 1))  # Subtract header
    
    log_info "Created inclusion list with ${inclusion_count} variants"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

create_chr_variant_map() {
    local chr="$1"
    local geno_dir="$2"
    local ldref_dir="$3"
    local step_dir="$4"
    
    local pvar_fmt="${geno_dir}/chr${chr}_pvar_fmt"
    local ld_rsids="${ldref_dir}/chr${chr}_ld_rsids"
    local out_map="${step_dir}/chr${chr}_variant_map"
    
    # Check if input files exist
    if [[ ! -f "$pvar_fmt" ]]; then
        log_debug "No genotype data for chr${chr}, skipping"
        return 0
    fi
    
    if [[ ! -f "$ld_rsids" ]]; then
        log_debug "No LD reference data for chr${chr}, skipping"
        return 0
    fi
    
    log_debug "Creating variant map for chr${chr}"
    
    # Sort files for join
    LC_ALL=C sort -k1,1 "$pvar_fmt" > "${step_dir}/chr${chr}_pvar_sorted.tmp"
    LC_ALL=C sort -k1,1 "$ld_rsids" > "${step_dir}/chr${chr}_ld_sorted.tmp"
    
    # Join on chr:pos (column 1)
    # pvar_fmt: chr:pos, pvar_a1, pvar_a2, pvar_snpid
    # ld_rsids: chr:pos, ld_a1, ld_a2, ld_rsid
    # Output: chr:pos, pvar_a1, pvar_a2, pvar_snpid, ld_a1, ld_a2, ld_rsid
    LC_ALL=C join -t$'\t' -1 1 -2 1 \
        -o 1.1,1.2,1.3,1.4,2.2,2.3,2.4 \
        "${step_dir}/chr${chr}_pvar_sorted.tmp" \
        "${step_dir}/chr${chr}_ld_sorted.tmp" \
        > "${step_dir}/chr${chr}_joined.tmp"
    
    # Filter for allele concordance (allow strand flips)
    awk -F'\t' -v OFS='\t' '
    BEGIN {
        # Complement mapping
        c["A"] = "T"; c["T"] = "A"; c["G"] = "C"; c["C"] = "G"
    }
    {
        chrpos = $1
        pvar_a1 = toupper($2); pvar_a2 = toupper($3); pvar_snpid = $4
        ld_a1 = toupper($5); ld_a2 = toupper($6); ld_rsid = $7
        
        # Check direct match
        direct_match = ((pvar_a1 == ld_a1 && pvar_a2 == ld_a2) || 
                       (pvar_a1 == ld_a2 && pvar_a2 == ld_a1))
        
        # Check strand flip match
        flip_match = ((c[pvar_a1] == ld_a1 && c[pvar_a2] == ld_a2) || 
                     (c[pvar_a1] == ld_a2 && c[pvar_a2] == ld_a1))
        
        if (direct_match || flip_match) {
            print chrpos, pvar_a1, pvar_a2, pvar_snpid, ld_a1, ld_a2, ld_rsid
        }
    }
    ' "${step_dir}/chr${chr}_joined.tmp" > "$out_map"
    
    # Clean up temp files
    rm -f "${step_dir}/chr${chr}_pvar_sorted.tmp" \
          "${step_dir}/chr${chr}_ld_sorted.tmp" \
          "${step_dir}/chr${chr}_joined.tmp"
    
    local map_count
    map_count=$(wc -l < "$out_map")
    log_debug "chr${chr}: ${map_count} variants mapped"
}

combine_variant_maps() {
    local step_dir="$1"
    local prep_dir="$2"
    
    # Add header - output to prep/ level (not inclusion_list/)
    echo -e "chrpos\tpvar_a1\tpvar_a2\tpvar_snpid\tld_a1\tld_a2\tld_rsid" > "${prep_dir}/variant_map.tsv"
    
    # Concatenate all chromosome maps
    for chr in $(get_chromosomes); do
        local chr_map="${step_dir}/chr${chr}_variant_map"
        if [[ -f "$chr_map" ]]; then
            cat "$chr_map" >> "${prep_dir}/variant_map.tsv"
        fi
    done
    
    local total_count
    total_count=$(wc -l < "${prep_dir}/variant_map.tsv")
    total_count=$((total_count - 1))  # Subtract header
    log_info "Combined variant map: ${total_count} variants"
    log_info "Variant map saved to: ${prep_dir}/variant_map.tsv"
}

create_final_inclusion_list() {
    local step_dir="$1"
    local prep_dir="$2"
    local info_threshold="$3"
    local maf_threshold="$4"
    
    # For now, the inclusion list is all variants in the map
    # INFO/MAF filtering would require additional data from genotypes
    # This can be enhanced later to include MAF from genotype files
    
    # Create inclusion list with essential columns
    # Format: ld_rsid, pvar_snpid, chrpos
    echo -e "ld_rsid\tpvar_snpid\tchrpos" > "${step_dir}/variant_inclusion_list.tsv"
    
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            # chrpos, pvar_a1, pvar_a2, pvar_snpid, ld_a1, ld_a2, ld_rsid
            print $7, $4, $1
        }
    ' "${prep_dir}/variant_map.tsv" >> "${step_dir}/variant_inclusion_list.tsv"
    
    # Create sorted RSID list for fast lookups (internal file)
    awk -F'\t' 'NR > 1 {print $1}' "${step_dir}/variant_inclusion_list.tsv" | \
        LC_ALL=C sort -u > "${step_dir}/.rsid_index"
    
    log_info "Created inclusion list and sorted RSID lookup file"
}








