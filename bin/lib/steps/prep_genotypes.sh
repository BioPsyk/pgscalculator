#!/bin/bash
# pgscalculator v2 - prep-genotypes step
# Extract variant IDs from genotype files and create formatted pvar files

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_prep_genotypes_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "sort" "sort is required for sorting variant lists"
    
    # Check config variables
    validate_required_config "CFG" "GENODIR" "GENOFILE" "OUTDIR"
    
    # Check genotype files exist
    require_file "${CFG_GENOFILE}" "Genotype manifest file not found"
    require_dir "${CFG_GENODIR}" "Genotype directory not found"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_prep_genotypes() {
    log_step "Running prep-genotypes"
    
    # Check dependencies
    check_prep_genotypes_deps
    
    # Set up output directory
    local outdir="${CFG_OUTDIR}"
    local step_dir
    step_dir=$(get_step_dir "$outdir" "genotypes")
    ensure_dir "$step_dir"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    log_substep "Reading genotype manifest: ${CFG_GENOFILE}"
    
    local genodir="${CFG_GENODIR}"
    local genofile="${CFG_GENOFILE}"
    local geno_format
    geno_format=$(detect_geno_format "$genofile")
    
    log_info "Detected genotype format: $geno_format"
    
    # Process based on format
    if [[ "$geno_format" == "plink2" ]]; then
        process_plink2_genotypes "$genofile" "$genodir" "$step_dir"
    elif [[ "$geno_format" == "plink1" ]]; then
        process_plink1_genotypes "$genofile" "$genodir" "$step_dir"
    else
        log_error "Unknown genotype format. Expected plink1 (bed/bim/fam) or plink2 (pgen/pvar/psam)"
        exit 1
    fi
    
    # Create combined sorted SNP list
    log_substep "Creating combined sorted variant list"
    create_combined_snplist "$step_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local snp_count
    snp_count=$(wc -l < "${step_dir}/snplist_sorted")
    log_info "Extracted ${snp_count} unique variants from genotypes"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

process_plink2_genotypes() {
    local genofile="$1"
    local genodir="$2"
    local step_dir="$3"
    
    log_substep "Processing PLINK2 format genotypes"
    
    for chr in $(get_chromosomes); do
        local pvar_file
        pvar_file=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "pvar")
        
        if [[ -z "$pvar_file" ]] || [[ ! -f "$pvar_file" ]]; then
            log_warn "No pvar file found for chromosome $chr, skipping"
            continue
        fi
        
        log_debug "Processing chr${chr}: $pvar_file"
        
        # Extract and format pvar
        # Output: chr:pos, a1, a2, variant_id
        local out_pvar="${step_dir}/chr${chr}_pvar_fmt"
        
        awk -F'\t' -v OFS='\t' '
            /^#/ { next }
            {
                # PVAR format: #CHROM POS ID REF ALT ...
                chrpos = $1 ":" $2
                print chrpos, $4, $5, $3
            }
        ' "$pvar_file" > "$out_pvar"
        
        # Also extract just the variant IDs for the combined list
        awk -F'\t' '/^#/ {next} {print $3}' "$pvar_file" >> "${step_dir}/snplist_unsorted"
        
        log_debug "Wrote formatted pvar for chr${chr}"
    done
}

process_plink1_genotypes() {
    local genofile="$1"
    local genodir="$2"
    local step_dir="$3"
    
    log_substep "Processing PLINK1 format genotypes"
    
    for chr in $(get_chromosomes); do
        local bim_file
        bim_file=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "bim")
        
        if [[ -z "$bim_file" ]] || [[ ! -f "$bim_file" ]]; then
            log_warn "No bim file found for chromosome $chr, skipping"
            continue
        fi
        
        log_debug "Processing chr${chr}: $bim_file"
        
        # Extract and format bim (similar structure to pvar)
        # BIM format: CHR, VAR_ID, CM, POS, A1, A2
        # Output: chr:pos, a1, a2, variant_id
        local out_pvar="${step_dir}/chr${chr}_pvar_fmt"
        
        awk -F'\t' -v OFS='\t' '{
            chrpos = $1 ":" $4
            print chrpos, $5, $6, $2
        }' "$bim_file" > "$out_pvar"
        
        # Also extract just the variant IDs for the combined list
        awk -F'\t' '{print $2}' "$bim_file" >> "${step_dir}/snplist_unsorted"
        
        log_debug "Wrote formatted pvar for chr${chr}"
    done
}

create_combined_snplist() {
    local step_dir="$1"
    
    # Sort and deduplicate the SNP list
    LC_ALL=C sort -u "${step_dir}/snplist_unsorted" > "${step_dir}/snplist_sorted"
    
    # Clean up temp file
    rm -f "${step_dir}/snplist_unsorted"
}




