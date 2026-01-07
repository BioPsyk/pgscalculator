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
    # Preferred (v2.1): filters.info_threshold / filters.maf_threshold
    # Backwards compatible: info_threshold / maf_threshold
    local info_threshold="${CFG_FILTERS_INFO_THRESHOLD:-${CFG_INFO_THRESHOLD:-0.8}}"
    local maf_threshold="${CFG_FILTERS_MAF_THRESHOLD:-${CFG_MAF_THRESHOLD:-0.01}}"
    
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
    
    local ref_dir="${prep_dir}/references"
    ensure_dir "$ref_dir"
    
    # Get reference file paths from config (optional)
    local info_file="${CFG_INFO_FILE:-}"
    local maf_file="${CFG_MAF_FILE:-}"

    # Allow "false" to explicitly disable these filters via config
    local info_forced_off="no"
    local maf_forced_off="no"
    if [[ -n "$info_file" ]] && [[ "${info_file,,}" == "false" ]]; then
        info_forced_off="yes"
        info_file=""
    fi
    if [[ -n "$maf_file" ]] && [[ "${maf_file,,}" == "false" ]]; then
        maf_forced_off="yes"
        maf_file=""
    fi
    
    # Decide whether we need MAF at all (threshold <= 0 disables MAF filtering)
    local maf_filter_enabled="yes"
    if awk -v t="${maf_threshold}" 'BEGIN{ exit !(t <= 0) }' 2>/dev/null; then
        maf_filter_enabled="no"
    fi
    if [[ "$maf_forced_off" == "yes" ]]; then
        maf_filter_enabled="no"
    fi

    # If MAF filtering enabled and no MAF file provided, compute from genotypes
    if [[ "$maf_filter_enabled" == "yes" ]] && ( [[ -z "$maf_file" ]] || [[ ! -f "$maf_file" ]] ); then
        log_substep "Computing MAF from genotypes"
        compute_maf_from_genotypes "$prep_dir" "$ref_dir"
        maf_file="${ref_dir}/maf_computed.tsv"
    fi
    
    # Start with all variants from variant_map
    local variant_map="${prep_dir}/variant_map.tsv"
    local input_count
    input_count=$(awk 'NR > 1' "$variant_map" | wc -l)
    log_info "Starting with ${input_count} variants from variant map"
    
    # Apply filters
    local tmpdir
    tmpdir=$(make_tmpdir "prep_inclusion_list")
    
    # Copy variant map to temp (add header for filtering output)
    cp "$variant_map" "${tmpdir}/variants.tsv"
    
    # INFO threshold <= 0 disables INFO filtering (even if info_file is provided)
    local info_filter_enabled="yes"
    if awk -v t="${info_threshold}" 'BEGIN{ exit !(t <= 0) }' 2>/dev/null; then
        info_filter_enabled="no"
    fi
    if [[ "$info_forced_off" == "yes" ]]; then
        info_filter_enabled="no"
    fi

    # Apply INFO filter if enabled and info_file provided
    local after_info_count="$input_count"
    if [[ "$info_filter_enabled" == "yes" ]] && [[ -n "$info_file" ]] && [[ -f "$info_file" ]]; then
        log_substep "Applying INFO filter (threshold: ${info_threshold})"
        apply_info_filter "${tmpdir}/variants.tsv" "$info_file" "$info_threshold" "${tmpdir}/after_info.tsv"
        after_info_count=$(awk 'NR > 1' "${tmpdir}/after_info.tsv" | wc -l)
        log_info "After INFO filter: ${after_info_count} variants"
        mv "${tmpdir}/after_info.tsv" "${tmpdir}/variants.tsv"
    elif [[ "$info_filter_enabled" != "yes" ]]; then
        if [[ "$info_forced_off" == "yes" ]]; then
            log_info "INFO filter disabled via config (info_file: false), skipping INFO filter"
        else
            log_info "INFO threshold <= 0, skipping INFO filter"
        fi
    else
        log_info "No INFO file provided, skipping INFO filter"
    fi
    
    # Apply MAF filter if enabled and maf_file exists
    local after_maf_count="$after_info_count"
    if [[ "$maf_filter_enabled" == "yes" ]] && [[ -f "$maf_file" ]]; then
        log_substep "Applying MAF filter (threshold: ${maf_threshold})"
        apply_maf_filter "${tmpdir}/variants.tsv" "$maf_file" "$maf_threshold" "${tmpdir}/after_maf.tsv"
        after_maf_count=$(awk 'NR > 1' "${tmpdir}/after_maf.tsv" | wc -l)
        log_info "After MAF filter: ${after_maf_count} variants"
        mv "${tmpdir}/after_maf.tsv" "${tmpdir}/variants.tsv"
    elif [[ "$maf_filter_enabled" != "yes" ]]; then
        if [[ "$maf_forced_off" == "yes" ]]; then
            log_info "MAF filter disabled via config (maf_file: false), skipping MAF filter"
        else
            log_info "MAF threshold <= 0, skipping MAF filter"
        fi
    else
        log_info "No MAF file available, skipping MAF filter"
    fi
    
    # Create inclusion list with essential columns
    # Format: ld_rsid, pvar_snpid, chrpos
    echo -e "ld_rsid\tpvar_snpid\tchrpos" > "${step_dir}/variant_inclusion_list.tsv"
    
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            # chrpos, pvar_a1, pvar_a2, pvar_snpid, ld_a1, ld_a2, ld_rsid
            print $7, $4, $1
        }
    ' "${tmpdir}/variants.tsv" >> "${step_dir}/variant_inclusion_list.tsv"
    
    # Create sorted RSID list for fast lookups (internal file)
    awk -F'\t' 'NR > 1 {print $1}' "${step_dir}/variant_inclusion_list.tsv" | \
        LC_ALL=C sort -u > "${step_dir}/.rsid_index"
    
    # Clean up
    rm -rf "$tmpdir"
    
    log_info "Created inclusion list with ${after_maf_count} variants"
}

compute_maf_from_genotypes() {
    local prep_dir="$1"
    local ref_dir="$2"
    
    local geno_dir="${CFG_GENODIR}"
    local geno_file="${CFG_GENOFILE}"
    local maf_output="${ref_dir}/maf_computed.tsv"
    
    # Check if plink2 is available
    if ! command -v plink2 &> /dev/null; then
        log_warn "plink2 not available, cannot compute MAF from genotypes"
        return 1
    fi
    
    # Header for MAF file
    echo -e "GENO_ID\tMAF" > "$maf_output"
    
    local total_variants=0
    
    for chr in $(get_chromosomes); do
        # Get genotype file for this chromosome
        local pgen
        pgen=$(get_geno_files_for_chr "$geno_file" "$geno_dir" "$chr" "pgen")
        
        if [[ -n "$pgen" ]] && [[ -f "$pgen" ]]; then
            local geno_prefix="${pgen%.pgen}"
            local tmpdir
            tmpdir=$(make_tmpdir "prep_inclusion_list_plink2_freq")
            
            # Compute allele frequencies
            plink2 --pfile "$geno_prefix" --freq --out "${tmpdir}/freq" \
                --threads 1 > "${tmpdir}/plink2.log" 2>&1 || true
            
            if [[ -f "${tmpdir}/freq.afreq" ]]; then
                # Extract ID and ALT_FREQS, convert to MAF
                awk -F'\t' -v OFS='\t' '
                    NR > 1 {
                        id = $2
                        alt_freq = $5
                        # Convert to MAF (0-0.5)
                        maf = (alt_freq > 0.5) ? (1 - alt_freq) : alt_freq
                        print id, maf
                    }
                ' "${tmpdir}/freq.afreq" >> "$maf_output"
                
                local chr_count
                chr_count=$(awk 'NR > 1' "${tmpdir}/freq.afreq" | wc -l)
                total_variants=$((total_variants + chr_count))
            fi
            
            rm -rf "$tmpdir"
        else
            # Try bed/bim/fam format
            local bed
            bed=$(get_geno_files_for_chr "$geno_file" "$geno_dir" "$chr" "bed")
            
            if [[ -n "$bed" ]] && [[ -f "$bed" ]]; then
                local geno_prefix="${bed%.bed}"
                local tmpdir
                tmpdir=$(make_tmpdir "prep_inclusion_list_plink2_freq")
                
                plink2 --bfile "$geno_prefix" --freq --out "${tmpdir}/freq" \
                    --threads 1 > "${tmpdir}/plink2.log" 2>&1 || true
                
                if [[ -f "${tmpdir}/freq.afreq" ]]; then
                    awk -F'\t' -v OFS='\t' '
                        NR > 1 {
                            id = $2
                            alt_freq = $5
                            maf = (alt_freq > 0.5) ? (1 - alt_freq) : alt_freq
                            print id, maf
                        }
                    ' "${tmpdir}/freq.afreq" >> "$maf_output"
                    
                    local chr_count
                    chr_count=$(awk 'NR > 1' "${tmpdir}/freq.afreq" | wc -l)
                    total_variants=$((total_variants + chr_count))
                fi
                
                rm -rf "$tmpdir"
            fi
        fi
    done
    
    log_info "Computed MAF for ${total_variants} variants"
}

apply_info_filter() {
    local input="$1"
    local info_file="$2"
    local threshold="$3"
    local output="$4"
    
    # INFO file format: GENO_ID, INFO
    # Variant map format: chrpos, pvar_a1, pvar_a2, pvar_snpid, ld_a1, ld_a2, ld_rsid
    # Filter: keep variants where INFO >= threshold
    
    awk -F'\t' -v OFS='\t' -v threshold="$threshold" '
        # Load INFO scores (GENO_ID -> INFO)
        ARGIND == 1 && FNR > 1 {
            info[$1] = $2
            next
        }
        # Process variant map
        ARGIND == 2 && FNR == 1 {
            print  # Header
            next
        }
        ARGIND == 2 {
            geno_id = $4  # pvar_snpid
            # Keep if no INFO data OR INFO >= threshold
            if (!(geno_id in info) || info[geno_id] >= threshold) {
                print
            }
        }
    ' "$info_file" "$input" > "$output"
}

apply_maf_filter() {
    local input="$1"
    local maf_file="$2"
    local threshold="$3"
    local output="$4"
    
    # MAF file format: GENO_ID, MAF
    # Filter: keep variants where MAF >= threshold
    
    awk -F'\t' -v OFS='\t' -v threshold="$threshold" '
        # Load MAF values (GENO_ID -> MAF)
        ARGIND == 1 && FNR > 1 {
            maf[$1] = $2
            next
        }
        # Process variant map
        ARGIND == 2 && FNR == 1 {
            print  # Header
            next
        }
        ARGIND == 2 {
            geno_id = $4  # pvar_snpid
            # Keep if no MAF data OR MAF >= threshold
            if (!(geno_id in maf) || maf[geno_id] >= threshold) {
                print
            }
        }
    ' "$maf_file" "$input" > "$output"
}











