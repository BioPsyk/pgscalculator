#!/bin/bash
# pgscalculator v2 - prep-ldref step
# Extract RSIDs from LD reference files

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_prep_ldref_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "sort" "sort is required for sorting"
    require_command "join" "join is required for liftover augmentation"
    require_command "zcat" "zcat is required for reading gzipped liftover reference"
    
    # Check config variables
    validate_required_config "CFG" "LDDIR" "OUTDIR"
    
    # Check LD reference directory exists
    require_dir "${CFG_LDDIR}" "LD reference directory not found"
    
    # Check liftover reference file exists (required for dual-position mapfile)
    local liftover_ref="${CFG_LIFTOVER_REFERENCE:-}"
    if [[ -z "$liftover_ref" ]]; then
        log_error "CFG_LIFTOVER_REFERENCE not set. This is required for the dual-position mapfile."
        exit 1
    fi
    require_file "$liftover_ref" "Liftover reference file is required for dual-position mapfile"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_prep_ldref() {
    log_step "Running prep-ldref"
    
    # Check dependencies
    check_prep_ldref_deps
    
    # Set up output directory
    local outdir="${CFG_OUTDIR}"
    local step_dir
    step_dir=$(get_step_dir "$outdir" "ldref")
    ensure_dir "$step_dir"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local lddir="${CFG_LDDIR}"
    
    log_substep "Processing LD reference files from: $lddir"
    
    local total_rsids=0
    
    for chr in $(get_chromosomes); do
        # Find the .info file for this chromosome
        # Pattern: *chr{N}*.info or *_{N}.info or similar
        local info_file
        info_file=$(find_ld_info_file "$lddir" "$chr")
        
        if [[ -z "$info_file" ]] || [[ ! -f "$info_file" ]]; then
            log_warn "No LD info file found for chromosome $chr, skipping"
            continue
        fi
        
        log_debug "Processing chr${chr}: $info_file"
        
        # Extract RSIDs and position info from LD reference
        # Info file format (sbayesR): Chrom ID GenPos PhysPos A1 A2 A1Freq
        local out_file="${step_dir}/chr${chr}_ld_rsids"
        
        awk -F' ' -v OFS='\t' '
            NR > 1 {
                # Output: chr:pos, a1, a2, rsid
                chrpos = $1 ":" $4
                print chrpos, $5, $6, $2
            }
        ' "$info_file" > "$out_file"
        
        local chr_count
        chr_count=$(wc -l < "$out_file")
        total_rsids=$((total_rsids + chr_count))
        
        log_debug "Extracted $chr_count variants from chr${chr}"
    done
    
    # Create combined LD RSID list
    log_substep "Creating combined LD RSID list"
    cat "${step_dir}"/chr*_ld_rsids | awk -F'\t' '{print $4}' | LC_ALL=C sort -u > "${step_dir}/ld_rsids_all"
    
    # Extract EAF from LD reference for use as fallback in sumstat processing
    log_substep "Extracting allele frequencies from LD reference"
    extract_ldref_eaf "$lddir" "$outdir"
    
    # Augment LD reference with GRCh38 positions (for dual-position mapfile)
    log_substep "Augmenting LD reference with GRCh38 positions"
    local liftover_ref="${CFG_LIFTOVER_REFERENCE}"
    local augmented_dir="${outdir}/prep/ldref_augmented"
    augment_ldref_with_liftover "$step_dir" "$liftover_ref" "$augmented_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    log_info "Extracted ${total_rsids} variants from LD reference"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

extract_ldref_eaf() {
    local lddir="$1"
    local outdir="$2"
    
    # Create references directory
    local ref_dir="${outdir}/prep/references"
    ensure_dir "$ref_dir"
    
    local eaf_file="${ref_dir}/ldref_eaf.tsv"
    
    # Header for EAF file
    echo -e "RSID\tA1\tA2\tA2Freq" > "$eaf_file"
    
    local total_variants=0
    
    for chr in $(get_chromosomes); do
        local info_file
        info_file=$(find_ld_info_file "$lddir" "$chr")
        
        if [[ -z "$info_file" ]] || [[ ! -f "$info_file" ]]; then
            continue
        fi
        
        # Extract: RSID (col 2), A1 (col 5), A2 (col 6), A2Freq (col 7)
        # Info file format (sbayesR): Chrom ID GenPos PhysPos A1 A2 A2Freq ...
        awk -F' ' -v OFS='\t' '
            NR > 1 {
                print $2, $5, $6, $7
            }
        ' "$info_file" >> "$eaf_file"
        
        local chr_count
        chr_count=$(awk 'NR > 1' "$info_file" | wc -l)
        total_variants=$((total_variants + chr_count))
    done
    
    log_info "Extracted EAF for ${total_variants} variants to: ${eaf_file}"
}

augment_ldref_with_liftover() {
    local ldref_dir="$1"
    local liftover_ref="$2"
    local augmented_dir="$3"
    
    ensure_dir "$augmented_dir"
    
    # Check if already augmented (all chromosomes present)
    local all_present=true
    for chr in $(get_chromosomes); do
        if [[ ! -f "${augmented_dir}/chr${chr}_ld_augmented.tsv" ]]; then
            all_present=false
            break
        fi
    done
    
    if [[ "$all_present" == true ]]; then
        log_info "Augmented LD reference already exists, skipping liftover join"
        return 0
    fi
    
    log_info "Joining LD reference with liftover to add GRCh38 positions"
    
    local augmented_count=0
    
    for chr in $(get_chromosomes); do
        local ld_file="${ldref_dir}/chr${chr}_ld_rsids"
        local out_file="${augmented_dir}/chr${chr}_ld_augmented.tsv"
        
        if [[ ! -f "$ld_file" ]]; then
            log_warn "No LD reference file for chr${chr}, skipping"
            continue
        fi
        
        log_debug "Augmenting chr${chr} with GRCh38 positions"
        
        # LD ref format: chr:pos_b37, a1, a2, rsid (tab-separated)
        # Liftover format: chr:pos_b37 chr:pos_b38 rsid a1 a2 (space-separated)
        #
        # Join on chr:pos_b37 (column 1 in both)
        # Output: pos_b37, pos_b38, ldref_a1, ldref_a2, ldref_rsid
        
        LC_ALL=C join -t $'\t' -1 1 -2 1 \
            <(LC_ALL=C sort -t $'\t' -k1,1 "$ld_file") \
            <(zcat "$liftover_ref" | grep "^${chr}:" | tr ' ' '\t' | LC_ALL=C sort -t $'\t' -k1,1) \
            2>/dev/null | \
            awk -F'\t' -v OFS='\t' '{
                # Input after join: pos_b37, ldref_a1, ldref_a2, ldref_rsid, pos_b38, liftover_rsid, liftover_a1, liftover_a2
                # Output: pos_b37, pos_b38, ldref_a1, ldref_a2, ldref_rsid
                print $1, $5, $2, $3, $4
            }' > "$out_file"
        
        local chr_count
        chr_count=$(wc -l < "$out_file")
        augmented_count=$((augmented_count + chr_count))
        log_debug "chr${chr}: ${chr_count} variants augmented with GRCh38 positions"
    done
    
    log_info "Augmented LD reference created: ${augmented_count} total variants"
    log_info "Output directory: ${augmented_dir}"
}

find_ld_info_file() {
    local lddir="$1"
    local chr="$2"
    
    # Try specific file patterns first (exact matches)
    local specific_files=(
        "${lddir}/band_chr${chr}.ldm.sparse.info"
        "${lddir}/chr${chr}.ldm.sparse.info"
        "${lddir}/ukb_chr${chr}.ldm.sparse.info"
    )
    
    for f in "${specific_files[@]}"; do
        if [[ -f "$f" ]]; then
            echo "$f"
            return 0
        fi
    done
    
    # Try glob patterns using compgen (safer than ls with globs)
    local pattern
    for pattern in "*chr${chr}.ldm.sparse.info" "*_chr${chr}.info" "*chr${chr}.info"; do
        local matches
        matches=$(compgen -G "${lddir}/${pattern}" 2>/dev/null | head -1) || true
        if [[ -n "$matches" ]] && [[ -f "$matches" ]]; then
            echo "$matches"
            return 0
        fi
    done
    
    # Broader search with find - be careful with chr number matching
    # Use word boundary to avoid chr2 matching chr22
    local found
    found=$(find "$lddir" -maxdepth 1 -name "*chr${chr}[._]*info" -o -name "*chr${chr}.info" 2>/dev/null | head -1) || true
    if [[ -n "$found" ]] && [[ -f "$found" ]]; then
        echo "$found"
        return 0
    fi
    
    # Return empty string (not an error - just no file found)
    echo ""
}




