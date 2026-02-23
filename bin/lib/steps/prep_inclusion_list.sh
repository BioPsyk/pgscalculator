#!/bin/bash
# pgscalculator v2 - prep-inclusion-list step
# Build variant map from all LD reference (liftover) variants; genotype columns NA where no match (left join from LD ref)

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
    
    # Step 5: Compute MAF from genotypes (runs once during prep, not per sumstat)
    log_substep "Computing MAF from genotypes"
    local ref_dir="${prep_dir}/references"
    ensure_dir "$ref_dir"
    compute_maf_from_genotypes "$prep_dir" "$ref_dir" "${prep_dir}/variant_map.tsv"
    
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
    # Row set: ALL augmented LD ref variants (left join from LD ref). Genotype columns NA where no match.
    # Match strategy: If genotype_build == GRCh37 match on pos_b37; if GRCh38 match on pos_b38.
    #
    # Output: chr, pos_b37, pos_b38, geno_snpid, geno_a1, geno_a2, ldref_snpid, ldref_a1, ldref_a2
    
    awk -F'\t' -v OFS='\t' -v build="$genotype_build" '
    BEGIN {
        c["A"] = "T"; c["T"] = "A"; c["G"] = "C"; c["C"] = "G"
    }
    # Load augmented LD reference (all rows; we output one row per LD ref)
    FNR==NR {
        pos_b37 = $1
        pos_b38 = $2
        la1 = toupper($3); la2 = toupper($4); lid = $5
        n_ld++
        ld_pos_b37[n_ld] = pos_b37
        ld_pos_b38[n_ld] = pos_b38
        ld_a1[n_ld] = la1
        ld_a2[n_ld] = la2
        ld_id[n_ld] = lid
        if (build == "GRCh38") {
            match_key[n_ld] = pos_b38
        } else {
            match_key[n_ld] = pos_b37
        }
        next
    }
    # Load genotypes into lookup (key -> gid, ga1, ga2); only for second file
    FNR != NR {
        chrpos = $1
        ga1 = toupper($2); ga2 = toupper($3); gid = $4
        k1 = chrpos SUBSEP ga1 SUBSEP ga2
        k2 = chrpos SUBSEP ga2 SUBSEP ga1
        fa1 = c[ga1]; fa2 = c[ga2]
        k3 = chrpos SUBSEP fa1 SUBSEP fa2
        k4 = chrpos SUBSEP fa2 SUBSEP fa1
        geno_id[k1] = gid; geno_a1[k1] = ga1; geno_a2[k1] = ga2
        geno_id[k2] = gid; geno_a1[k2] = ga1; geno_a2[k2] = ga2
        geno_id[k3] = gid; geno_a1[k3] = fa1; geno_a2[k3] = fa2
        geno_id[k4] = gid; geno_a1[k4] = fa1; geno_a2[k4] = fa2
        next
    }
    END {
        for (i = 1; i <= n_ld; i++) {
            mk = match_key[i]
            la1 = ld_a1[i]; la2 = ld_a2[i]
            k1 = mk SUBSEP la1 SUBSEP la2
            k2 = mk SUBSEP la2 SUBSEP la1
            fa1 = c[la1]; fa2 = c[la2]
            k3 = mk SUBSEP fa1 SUBSEP fa2
            k4 = mk SUBSEP fa2 SUBSEP fa1
            gid = "NA"; ga1 = "NA"; ga2 = "NA"
            if (k1 in geno_id) { gid = geno_id[k1]; ga1 = geno_a1[k1]; ga2 = geno_a2[k1] }
            else if (k2 in geno_id) { gid = geno_id[k2]; ga1 = geno_a1[k2]; ga2 = geno_a2[k2] }
            else if (k3 in geno_id) { gid = geno_id[k3]; ga1 = geno_a1[k3]; ga2 = geno_a2[k3] }
            else if (k4 in geno_id) { gid = geno_id[k4]; ga1 = geno_a1[k4]; ga2 = geno_a2[k4] }
            split(ld_pos_b37[i], pb37, ":")
            split(ld_pos_b38[i], pb38, ":")
            chr_v = pb37[1]
            pos_b37_v = pb37[2]
            pos_b38_v = pb38[2]
            print chr_v, pos_b37_v, pos_b38_v, gid, ga1, ga2, ld_id[i], ld_a1[i], ld_a2[i]
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
        awk -F'\t' 'FNR==1{next} {print}' "${in_dir}"/chr*.tsv >> "${prep_dir}/variant_map.tsv"
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
    
    # Create inclusion list with essential columns (only rows with genotype match; used for scoring)
    # Schema: chr(1), pos_b37(2), pos_b38(3), geno_snpid(4), geno_a1(5), geno_a2(6), ldref_snpid(7), ...
    # Format: ldref_snpid, geno_snpid, chr, pos_b37, pos_b38
    echo -e "ldref_snpid\tgeno_snpid\tchr\tpos_b37\tpos_b38" > "${step_dir}/variant_inclusion_list.tsv"
    awk -F'\t' -v OFS='\t' 'NR > 1 && $4 != "NA" { print $7, $4, $1, $2, $3 }' "$variant_map" >> "${step_dir}/variant_inclusion_list.tsv"
    
    # Create sorted LDREF list for fast lookups (internal file; only rows with genotype match)
    awk -F'\t' 'NR > 1 && $4 != "NA" && $7 != "NA" {print $7}' "$variant_map" | \
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
    
    # Compute MAF from genotypes (runs once during prep, not per sumstat)
    log_substep "Computing MAF from genotypes"
    local ref_dir="${prep_dir}/references"
    ensure_dir "$ref_dir"
    compute_maf_from_genotypes "$prep_dir" "$ref_dir" "${prep_dir}/variant_map.tsv"
    
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
        log_error "variant_map.tsv is required to compute MAF (rows with genotype match)."
        log_error "Got: ${variant_map:-<empty>}"
        return 1
    fi

    if [[ -f "$maf_output" ]] && [[ $(wc -l < "$maf_output") -gt 1 ]]; then
        log_info "maf_computed.tsv already exists with data, skipping (use --force to recompute)"
        return 0
    fi
    
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
                $4 != "NA" && $1 == c {
                    # Schema: chr(1), pos_b37(2), pos_b38(3), geno_snpid(4), ...
                    print $4
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
                awk -F'\t' -v OFS='\t' '
                    NR == 1 {
                        for (i = 1; i <= NF; i++) {
                            if ($i == "ID") id_col = i
                            if ($i == "ALT_FREQS") af_col = i
                        }
                        next
                    }
                    {
                        af = $af_col + 0
                        maf = (af > 0.5) ? (1 - af) : af
                        print $id_col, maf
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
                        NR == 1 {
                            for (i = 1; i <= NF; i++) {
                                if ($i == "ID") id_col = i
                                if ($i == "ALT_FREQS") af_col = i
                            }
                            next
                        }
                        {
                            af = $af_col + 0
                            maf = (af > 0.5) ? (1 - af) : af
                            print $id_col, maf
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
    
    # Count final variant map (all LD ref liftover variants; genotype NA where no match)
    local n_variant_map=0
    if [[ -f "${prep_dir}/variant_map.tsv" ]]; then
        n_variant_map=$(wc -l < "${prep_dir}/variant_map.tsv")
        n_variant_map=$((n_variant_map - 1))  # Subtract header
    fi
    log_info "Variant map (all LD ref liftover): ${n_variant_map}"

    local n_geno_matched=0
    if [[ -f "${prep_dir}/variant_map.tsv" ]]; then
        n_geno_matched=$(awk -F'\t' 'NR > 1 && $4 != "NA" {count++} END {print count+0}' "${prep_dir}/variant_map.tsv")
    fi
    log_info "Variants with genotype match: ${n_geno_matched}"
    
    # Calculate match rate (genotype match vs augmented LD ref)
    local match_pct="0.00"
    if [[ "$n_ldref_augmented" -gt 0 ]]; then
        match_pct=$(awk "BEGIN {printf \"%.2f\", 100 * ${n_geno_matched} / ${n_ldref_augmented}}")
    fi
    
    # Read ldref steps if available
    local n_liftover=0
    local n_ldref=0
    local n_ldref_after_join=0
    local ldref_steps="${details_dir}/prep_ldref_steps.tsv"
    if [[ -f "$ldref_steps" ]]; then
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
        echo -e "geno-ldref-map\t${n_ldref_after_join}\t${n_variant_map}\tall LD ref in map; geno_matched=${n_geno_matched} (${match_pct}%)"
    } > "$steps_file"
    
    log_info "Wrote prep stepwise details: ${steps_file}"
}








