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
    
    if [[ ! -f "$pvar_fmt" ]] || [[ ! -f "$ld_rsids" ]]; then
        log_debug "Missing genotype or LD reference data for chr${chr}, skipping"
        return 0
    fi
    
    log_debug "Creating variant map for chr${chr}"
    
    awk -F'\t' -v OFS='\t' '
    BEGIN {
        c["A"] = "T"; c["T"] = "A"; c["G"] = "C"; c["C"] = "G"
    }
    FNR==NR {
        chrpos = $1
        a1 = toupper($2); a2 = toupper($3); id = $4
        key = chrpos SUBSEP a1 SUBSEP a2
        pvar_id[key] = id
        pvar_a1[key] = a1
        pvar_a2[key] = a2
        next
    }
    {
        chrpos = $1
        la1 = toupper($2); la2 = toupper($3); lid = $4
        # try to match to a genotype entry (direct or strand flip, either orientation)
        matched = ""
        # direct / swap
        k1 = chrpos SUBSEP la1 SUBSEP la2
        k2 = chrpos SUBSEP la2 SUBSEP la1
        # strand flip
        fa1 = c[la1]; fa2 = c[la2]
        k3 = chrpos SUBSEP fa1 SUBSEP fa2
        k4 = chrpos SUBSEP fa2 SUBSEP fa1
        if (k1 in pvar_id) matched = k1
        else if (k2 in pvar_id) matched = k2
        else if (k3 in pvar_id) matched = k3
        else if (k4 in pvar_id) matched = k4
        
        split(chrpos, cp, ":")
        chr_v = cp[1]; pos_v = cp[2]
        if (matched != "") {
            if (!(matched in pvar_matched)) {
                pvar_matched[matched] = 1
                print chr_v, pos_v, pvar_id[matched], pvar_a1[matched], pvar_a2[matched], lid, la1, la2
            }
        }
    }
    END {
        # intersection only: no unmatched geno or ldref rows
    }
    ' "$pvar_fmt" "$ld_rsids" > "$out_map"
    
    local map_count
    map_count=$(wc -l < "$out_map")
    log_debug "chr${chr}: ${map_count} variants mapped"
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
        
        echo -e "chr\tpos\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq" > "$out_map"
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
                    ldid = $6
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
    
    # Add header - output to prep/ level (not inclusion_list/)
    echo -e "chr\tpos\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq" > "${prep_dir}/variant_map.tsv"
    
    # Concatenate all chromosome maps (skip header)
    if compgen -G "${in_dir}/chr*.tsv" > /dev/null; then
        awk -F'\t' 'NR==1{next} {print}' "${in_dir}"/chr*.tsv >> "${prep_dir}/variant_map.tsv"
    fi
    
    local total_count
    total_count=$(wc -l < "${prep_dir}/variant_map.tsv")
    total_count=$((total_count - 1))  # Subtract header
    log_info "Combined variant map: ${total_count} variants"
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
    # Format: ldref_snpid, geno_snpid, chr, pos
    echo -e "ldref_snpid\tgeno_snpid\tchr\tpos" > "${step_dir}/variant_inclusion_list.tsv"
    awk -F'\t' -v OFS='\t' 'NR > 1 { print $6, $3, $1, $2 }' "$variant_map" >> "${step_dir}/variant_inclusion_list.tsv"
    
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
        # variant_map format: chr, pos, geno_snpid, geno_a1, geno_a2, ldref_snpid, ldref_a1, ldref_a2, ldref_a2freq
        local tmpdir
        tmpdir=$(make_tmpdir "prep_inclusion_list_plink2_freq")
        local extract_ids=""
        if [[ -n "$variant_map" ]] && [[ -f "$variant_map" ]]; then
            extract_ids="${tmpdir}/chr${chr}.extract_ids.txt"
            awk -F'\t' -v c="$chr" '
                NR==1{next}
                {
                    if ($1==c) print $3
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











