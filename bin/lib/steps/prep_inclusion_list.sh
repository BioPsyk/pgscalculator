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
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    
    # Check that prep-genotypes has been run
    require_dir "${prep_dir}/genotypes" "Run 'pgscalculator prep-genotypes' first"
    require_file "${prep_dir}/genotypes/snplist_sorted" "Run 'pgscalculator prep-genotypes' first"
    
    # Check that prep-ldref has been run (with augmented LD reference)
    require_dir "${prep_dir}/ldref" "Run 'pgscalculator prep-ldref' first"
    require_dir "${prep_dir}/ldref_augmented" "Run 'pgscalculator prep-ldref' first (augmented LD ref missing)"
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
    
    local geno_dir="${prep_dir}/genotypes"
    local ldref_dir="${prep_dir}/ldref"
    
    # If running per-chromosome (driver/array), only build that chromosome map
    local specific_chr=""
    if [[ -n "${CFG_CHROMOSOMES:-}" ]] && [[ "${CFG_CHROMOSOMES}" =~ ^(chr)?[0-9]+$ ]]; then
        specific_chr="${CFG_CHROMOSOMES#chr}"
    fi
    
    if [[ -n "$specific_chr" ]]; then
        local chr_marker="${step_dir}/.completed_chr${specific_chr}"
        if [[ -f "$chr_marker" ]]; then
            log_info "prep-inclusion-list chr${specific_chr} already completed. Use --force to re-run."
            return 0
        fi
        create_chr_variant_map "$specific_chr" "$geno_dir" "$ldref_dir" "$step_dir"
        date '+%Y-%m-%d %H:%M:%S' > "$chr_marker"
        log_info "Completed prep-inclusion-list for chr${specific_chr}"
        return 0
    fi
    
    # Check if already completed (full run)
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    # Step 1: Create per-chromosome variant maps (parallel by chromosome)
    log_substep "Creating per-chromosome variant maps"
    local max_parallel=""
    if [[ -n "${CFG_SLURM_PREP_MAX_PARALLEL:-}" ]]; then
        max_parallel="${CFG_SLURM_PREP_MAX_PARALLEL}"
    elif [[ -n "${CFG_SLURM_PREP:-}" ]]; then
        # Parse inline dict, e.g. "{ mem: 10g, cpus: 6, time: '1:00:00', max_parallel: 22 }"
        max_parallel=$(echo "${CFG_SLURM_PREP}" | sed 's/[{}]//g' | tr ',' '\n' | \
            awk -F': ' '$1 ~ /max_parallel/ {gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}')
    elif [[ -n "${CFG_PREP_MAX_PARALLEL:-}" ]]; then
        # Backward compatibility (legacy top-level key)
        max_parallel="${CFG_PREP_MAX_PARALLEL}"
    else
        max_parallel="4"
    fi
    if ! [[ "$max_parallel" =~ ^[0-9]+$ ]] || [[ "$max_parallel" -lt 1 ]]; then
        log_warn "Invalid prep_max_parallel (${max_parallel}), falling back to 1"
        max_parallel=1
    fi
    log_info "prep-inclusion-list: using up to ${max_parallel} parallel chr jobs"
    
    local -a pids=()
    local job_fail=0
    for chr in $(get_chromosomes); do
        create_chr_variant_map "$chr" "$geno_dir" "$ldref_dir" "$step_dir" &
        pids+=($!)
        
        if [[ "${#pids[@]}" -ge "$max_parallel" ]]; then
            if ! wait -n; then
                job_fail=1
            fi
            # prune finished pids
            local -a still_running=()
            local pid
            for pid in "${pids[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    still_running+=("$pid")
                fi
            done
            pids=("${still_running[@]}")
        fi
    done
    
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            job_fail=1
        fi
    done
    
    if [[ "$job_fail" -ne 0 ]]; then
        log_error "One or more per-chromosome map jobs failed"
        return 1
    fi
    
    # Step 2: Write per-chromosome mapfiles with ldref_a2freq (prep/variant_map/)
    log_substep "Writing per-chromosome mapfiles"
    write_chr_variant_maps "$step_dir" "$prep_dir"
    
    # Step 3: Combine all chromosome maps (output to prep/ level)
    log_substep "Combining chromosome variant maps"
    combine_variant_maps "$prep_dir"
    
    # Step 4: Create final inclusion list (derived from mapfile)
    log_substep "Creating final variant inclusion list"
    create_final_inclusion_list "$step_dir" "$prep_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local inclusion_count
    inclusion_count=$(wc -l < "${step_dir}/variant_inclusion_list.tsv")
    inclusion_count=$((inclusion_count - 1))  # Subtract header
    
    log_info "Created inclusion list with ${inclusion_count} variants"
    log_info "Output directory: ${step_dir}"
    
    # Generate combined prep stepwise details
    generate_prep_stepwise_details "$prep_dir" "$geno_dir"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

create_chr_variant_map() {
    local chr="$1"
    local geno_dir="$2"
    local ldref_dir="$3"
    local step_dir="$4"
    
    # Get genotype build from config (default: GRCh37)
    local genotype_build="${CFG_GENOTYPE_BUILD:-GRCh37}"
    
    local pvar_fmt="${geno_dir}/chr${chr}_pvar_fmt"
    # Use augmented LD reference (has both pos_b37 and pos_b38)
    local augmented_dir="${ldref_dir}/../ldref_augmented"
    local ld_augmented="${augmented_dir}/chr${chr}_ld_augmented.tsv"
    local out_map="${step_dir}/chr${chr}_variant_map"
    
    # Debug: Log the paths being checked
    log_info "prep-inclusion-list chr${chr}: pvar_fmt=${pvar_fmt}"
    log_info "prep-inclusion-list chr${chr}: ld_augmented=${ld_augmented}"
    
    if [[ ! -f "$pvar_fmt" ]]; then
        log_warn "Missing genotype data for chr${chr}: ${pvar_fmt}"
        return 0
    fi
    
    if [[ ! -f "$ld_augmented" ]]; then
        log_warn "Missing augmented LD reference for chr${chr}: ${ld_augmented}"
        return 0
    fi
    
    log_info "Creating variant map for chr${chr} (genotype_build: ${genotype_build})"
    
    # Augmented LD ref format: pos_b37, pos_b38, ldref_a1, ldref_a2, ldref_rsid (tab-separated)
    # Genotype format: chr:pos, a1, a2, snpid (tab-separated)
    #
    # Match strategy:
    #   - If genotype_build == GRCh37: match geno chr:pos to ld pos_b37
    #   - If genotype_build == GRCh38: match geno chr:pos to ld pos_b38
    #
    # Output: chr, pos_b37, pos_b38, geno_snpid, geno_a1, geno_a2, ldref_snpid, ldref_a1, ldref_a2
    
    awk -F'\t' -v OFS='\t' -v build="$genotype_build" '
    BEGIN {
        c["A"] = "T"; c["T"] = "A"; c["G"] = "C"; c["C"] = "G"
    }
    # Load augmented LD reference
    FNR==NR {
        # Format: pos_b37, pos_b38, ldref_a1, ldref_a2, ldref_rsid
        pos_b37 = $1
        pos_b38 = $2
        la1 = toupper($3); la2 = toupper($4); lid = $5
        
        # Determine match key based on which position we will join on
        if (build == "GRCh38") {
            match_key = pos_b38
        } else {
            match_key = pos_b37
        }
        
        # Store with allele variations for matching (all 4 orientations)
        # direct
        k1 = match_key SUBSEP la1 SUBSEP la2
        k2 = match_key SUBSEP la2 SUBSEP la1
        # strand flip
        fa1 = c[la1]; fa2 = c[la2]
        k3 = match_key SUBSEP fa1 SUBSEP fa2
        k4 = match_key SUBSEP fa2 SUBSEP fa1
        
        # Store data for each key variant (all 4 allele orientations)
        ld_pos_b37[k1] = pos_b37; ld_pos_b38[k1] = pos_b38
        ld_a1[k1] = la1; ld_a2[k1] = la2; ld_id[k1] = lid
        
        ld_pos_b37[k2] = pos_b37; ld_pos_b38[k2] = pos_b38
        ld_a1[k2] = la1; ld_a2[k2] = la2; ld_id[k2] = lid
        
        ld_pos_b37[k3] = pos_b37; ld_pos_b38[k3] = pos_b38
        ld_a1[k3] = fa1; ld_a2[k3] = fa2; ld_id[k3] = lid
        
        ld_pos_b37[k4] = pos_b37; ld_pos_b38[k4] = pos_b38
        ld_a1[k4] = fa1; ld_a2[k4] = fa2; ld_id[k4] = lid
        next
    }
    # Process genotypes
    {
        chrpos = $1
        ga1 = toupper($2); ga2 = toupper($3); gid = $4
        
        # Try all allele orientations
        k1 = chrpos SUBSEP ga1 SUBSEP ga2
        k2 = chrpos SUBSEP ga2 SUBSEP ga1
        fa1 = c[ga1]; fa2 = c[ga2]
        k3 = chrpos SUBSEP fa1 SUBSEP fa2
        k4 = chrpos SUBSEP fa2 SUBSEP fa1
        
        matched = ""
        if (k1 in ld_id) matched = k1
        else if (k2 in ld_id) matched = k2
        else if (k3 in ld_id) matched = k3
        else if (k4 in ld_id) matched = k4
        
        if (matched != "" && !(matched in already_matched)) {
            already_matched[matched] = 1
            
            # Extract chr from chrpos
            split(chrpos, cp, ":")
            chr_v = cp[1]
            
            # Extract position numbers from chr:pos strings
            split(ld_pos_b37[matched], pb37, ":")
            split(ld_pos_b38[matched], pb38, ":")
            pos_b37_v = pb37[2]
            pos_b38_v = pb38[2]
            
            # Output: chr, pos_b37, pos_b38, geno_snpid, geno_a1, geno_a2, ldref_snpid, ldref_a1, ldref_a2
            print chr_v, pos_b37_v, pos_b38_v, gid, ga1, ga2, ld_id[matched], ld_a1[matched], ld_a2[matched]
        }
    }
    ' "$ld_augmented" "$pvar_fmt" > "$out_map"
    
    local map_count
    map_count=$(wc -l < "$out_map")
    log_info "chr${chr}: ${map_count} variants mapped to ${out_map}"
}

write_chr_variant_maps() {
    local step_dir="$1"
    local prep_dir="$2"
    local ref_dir="${prep_dir}/references"
    local ldref_eaf="${ref_dir}/ldref_eaf.tsv"
    local out_dir="${prep_dir}/variant_map"
    
    ensure_dir "$out_dir"
    
    for chr in $(get_chromosomes); do
        local in_map="${step_dir}/chr${chr}_variant_map"
        local out_map="${out_dir}/chr${chr}.tsv"
        
        if [[ ! -s "$in_map" ]]; then
            continue
        fi
        
        # New header with dual positions (pos_b37 and pos_b38)
        echo -e "chr\tpos_b37\tpos_b38\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq" > "$out_map"
        if [[ -f "$ldref_eaf" ]]; then
            awk -F'\t' -v OFS='\t' -v eaf_file="$ldref_eaf" '
                BEGIN{
                    while ((getline < eaf_file) > 0) {
                        if (NR==1) continue
                        eaf[$1]=$4
                    }
                    close(eaf_file)
                }
                {
                    # in_map format: chr, pos_b37, pos_b38, geno_snpid, geno_a1, geno_a2, ldref_snpid, ldref_a1, ldref_a2
                    ldid = $7  # ldref_snpid is now column 7
                    freq = (ldid != "NA" && (ldid in eaf)) ? eaf[ldid] : "NA"
                    print $0, freq
                }
            ' "$in_map" >> "$out_map"
        else
            awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "$in_map" >> "$out_map"
        fi
    done
}

combine_variant_maps() {
    local prep_dir="$1"
    local in_dir="${prep_dir}/variant_map"
    
    # New header with dual positions (pos_b37 and pos_b38)
    echo -e "chr\tpos_b37\tpos_b38\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq" > "${prep_dir}/variant_map.tsv"
    
    # Concatenate all chromosome maps (skip header)
    if compgen -G "${in_dir}/chr*.tsv" > /dev/null; then
        awk -F'\t' 'NR==1{next} {print}' "${in_dir}"/chr*.tsv >> "${prep_dir}/variant_map.tsv"
    fi
    
    local total_count
    total_count=$(wc -l < "${prep_dir}/variant_map.tsv")
    total_count=$((total_count - 1))  # Subtract header
    log_info "Combined variant map: ${total_count} variants (with dual positions)"
    log_info "Variant map saved to: ${prep_dir}/variant_map.tsv"
}

create_final_inclusion_list() {
    local step_dir="$1"
    local prep_dir="$2"
    
    local variant_map="${prep_dir}/variant_map.tsv"
    if [[ ! -s "$variant_map" ]]; then
        log_error "Missing or empty variant_map.tsv; cannot derive inclusion list."
        log_error "Expected: ${variant_map}"
        return 1
    fi
    
    local input_count
    input_count=$(awk 'NR > 1' "$variant_map" | wc -l)
    log_info "Starting with ${input_count} variants from variant map"
    
    # Create inclusion list with essential columns
    # New schema: chr(1), pos_b37(2), pos_b38(3), geno_snpid(4), geno_a1(5), geno_a2(6), ldref_snpid(7), ...
    # Format: ldref_snpid, geno_snpid, chr, pos_b37, pos_b38
    echo -e "ldref_snpid\tgeno_snpid\tchr\tpos_b37\tpos_b38" > "${step_dir}/variant_inclusion_list.tsv"
    awk -F'\t' -v OFS='\t' 'NR > 1 { print $7, $4, $1, $2, $3 }' "$variant_map" >> "${step_dir}/variant_inclusion_list.tsv"
    
    # Create sorted LDREF list for fast lookups (internal file)
    awk -F'\t' 'NR > 1 && $1 != "NA" {print $1}' "${step_dir}/variant_inclusion_list.tsv" | \
        LC_ALL=C sort -u > "${step_dir}/.rsid_index"
    
    log_info "Created inclusion list with ${input_count} variants"
}

run_prep_inclusion_list_combine() {
    log_step "Running prep-inclusion-list (combine)"
    
    check_prep_inclusion_list_deps
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local step_dir
    step_dir=$(get_step_dir "$outdir" "inclusion_list")
    ensure_dir "$step_dir"
    
    log_substep "Writing per-chromosome mapfiles"
    write_chr_variant_maps "$step_dir" "$prep_dir"
    
    log_substep "Combining chromosome variant maps"
    combine_variant_maps "$prep_dir"
    
    log_substep "Creating final variant inclusion list"
    create_final_inclusion_list "$step_dir" "$prep_dir"
    
    mark_step_completed "$step_dir"
    log_info "prep-inclusion-list combine completed"
}

compute_maf_from_genotypes() {
    local prep_dir="$1"
    local ref_dir="$2"
    local variant_map="$3"
    
    local geno_dir="${CFG_GENODIR}"
    local geno_file="${CFG_GENOFILE}"
    local maf_output="${ref_dir}/maf_computed.tsv"

    if [[ -z "$variant_map" ]] || [[ ! -s "$variant_map" ]]; then
        log_error "variant_map.tsv is required to compute MAF (genotype ∩ LD-ref)."
        log_error "Got: ${variant_map:-<empty>}"
        return 1
    fi
    
    # Check if plink2 is available
    if ! command -v plink2 &> /dev/null; then
        log_warn "plink2 not available, cannot compute MAF from genotypes"
        return 1
    fi
    
    # Header for MAF file
    echo -e "GENO_ID\tMAF" > "$maf_output"
    
    local total_variants=0
    
    for chr in $(get_chromosomes); do
        # Optional: restrict MAF computation to variants in the combined variant map.
        # variant_map format: chr, pos_b37, pos_b38, geno_snpid, geno_a1, geno_a2, ldref_snpid, ldref_a1, ldref_a2, ldref_a2freq
        local tmpdir
        tmpdir=$(make_tmpdir "prep_inclusion_list_plink2_freq")
        local extract_ids=""
        if [[ -n "$variant_map" ]] && [[ -f "$variant_map" ]]; then
            extract_ids="${tmpdir}/chr${chr}.extract_ids.txt"
            awk -F'\t' -v c="$chr" '
                NR==1{next}
                {
                    # Schema: chr(1), pos_b37(2), pos_b38(3), geno_snpid(4), ...
                    if ($1==c) print $4
                }
            ' "$variant_map" | LC_ALL=C sort -u > "$extract_ids"
            # If no variants for this chromosome, skip plink entirely
            if [[ ! -s "$extract_ids" ]]; then
                rm -rf "$tmpdir"
                continue
            fi
        fi

        # Get genotype file for this chromosome
        local pgen
        pgen=$(get_geno_files_for_chr "$geno_file" "$geno_dir" "$chr" "pgen")
        
        if [[ -n "$pgen" ]] && [[ -f "$pgen" ]]; then
            local geno_prefix="${pgen%.pgen}"
            
            # Compute allele frequencies
            local plink_threads="${CFG_PLINK_THREADS:-1}"
            if ! [[ "$plink_threads" =~ ^[0-9]+$ ]] || [[ "$plink_threads" -lt 1 ]]; then
                plink_threads=1
            fi
            local plink_cmd=(plink2 --pfile "$geno_prefix" --freq --out "${tmpdir}/freq" --threads "${plink_threads}")
            if [[ -n "$extract_ids" ]]; then
                plink_cmd+=(--extract "$extract_ids")
            fi
            "${plink_cmd[@]}" > "${tmpdir}/plink2.log" 2>&1 || true
            
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
                
                local plink_threads="${CFG_PLINK_THREADS:-1}"
                if ! [[ "$plink_threads" =~ ^[0-9]+$ ]] || [[ "$plink_threads" -lt 1 ]]; then
                    plink_threads=1
                fi
                local plink_cmd=(plink2 --bfile "$geno_prefix" --freq --out "${tmpdir}/freq" --threads "${plink_threads}")
                if [[ -n "$extract_ids" ]]; then
                    plink_cmd+=(--extract "$extract_ids")
                fi
                "${plink_cmd[@]}" > "${tmpdir}/plink2.log" 2>&1 || true
                
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
            else
                # No genotype files for this chromosome: clean up tempdir created earlier
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

generate_prep_stepwise_details() {
    local prep_dir="$1"
    local geno_dir="$2"
    
    local details_dir="${prep_dir}/details"
    ensure_dir "$details_dir"
    
    local steps_file="${details_dir}/steps.tsv"
    
    # Count genotype variants
    local n_geno=0
    if [[ -f "${geno_dir}/snplist_sorted" ]]; then
        n_geno=$(wc -l < "${geno_dir}/snplist_sorted")
    fi
    log_info "Genotype variants: ${n_geno}"
    
    # Count augmented LD ref variants
    local n_ldref_augmented=0
    local augmented_dir="${prep_dir}/ldref_augmented"
    for chr in $(get_chromosomes); do
        local aug_file="${augmented_dir}/chr${chr}_ld_augmented.tsv"
        if [[ -f "$aug_file" ]]; then
            local chr_count
            chr_count=$(wc -l < "$aug_file")
            n_ldref_augmented=$((n_ldref_augmented + chr_count))
        fi
    done
    log_info "Augmented LD reference variants: ${n_ldref_augmented}"
    
    # Count final variant map (intersection of genotypes and augmented LD ref)
    local n_variant_map=0
    if [[ -f "${prep_dir}/variant_map.tsv" ]]; then
        n_variant_map=$(wc -l < "${prep_dir}/variant_map.tsv")
        n_variant_map=$((n_variant_map - 1))  # Subtract header
    fi
    log_info "Variant map (geno ∩ ldref): ${n_variant_map}"
    
    # Calculate match rate
    local match_pct="0.00"
    if [[ "$n_ldref_augmented" -gt 0 ]]; then
        match_pct=$(awk "BEGIN {printf \"%.2f\", 100 * ${n_variant_map} / ${n_ldref_augmented}}")
    fi
    
    # Read ldref steps if available
    local n_liftover=0
    local n_ldref=0
    local n_ldref_after_join=0
    local ldref_steps="${details_dir}/prep_ldref_steps.tsv"
    if [[ -f "$ldref_steps" ]]; then
        # Parse the ldref steps file to get counts
        n_liftover=$(awk -F'\t' '$1=="liftover-reference" {print $2}' "$ldref_steps")
        n_ldref=$(awk -F'\t' '$1=="ldref-extract" {print $2}' "$ldref_steps")
        n_ldref_after_join=$(awk -F'\t' '$1=="ldref-augment" {print $3}' "$ldref_steps")
    fi
    
    # Write combined prep steps file
    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "liftover-reference\t${n_liftover}\t${n_liftover}\tliftover reference file rows"
        echo -e "ldref-extract\t${n_ldref}\t${n_ldref}\tLD reference variants extracted"
        echo -e "ldref-augment\t${n_ldref}\t${n_ldref_after_join}\tLD ref ↔ liftover join"
        echo -e "genotypes-extract\t${n_geno}\t${n_geno}\tgenotype variants extracted"
        echo -e "geno-ldref-intersect\t${n_ldref_after_join}\t${n_variant_map}\tgeno ∩ augmented-ldref (match_rate=${match_pct}%)"
    } > "$steps_file"
    
    log_info "Wrote prep stepwise details: ${steps_file}"
}








