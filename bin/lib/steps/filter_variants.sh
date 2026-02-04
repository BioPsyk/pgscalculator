#!/bin/bash
# pgscalculator v2 - filter-variants step
# Filter sumstat to inclusion list variants, then derive N/EAF/B/SE
# (Early filtering reduces data volume before expensive stat derivation)

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_filter_variants_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "sort" "sort is required for sorting"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    
    # Check that prep-inclusion-list has been run
    require_file "${prep_dir}/variant_map.tsv" "Run 'pgscalculator prep-inclusion-list' first"
}

# Validate sumstat metadata has required N field for posterior calculations.
# This is a safety net for non-sbatch runs (sbatch runs check earlier in pgscalculator-v2.sh).
validate_sumstat_n_field() {
    local metadata_file="$1"
    local which_n="$2"
    
    if [[ ! -f "$metadata_file" ]]; then
        log_error "Metadata file not found: ${metadata_file}"
        log_error "The sumstat directory must contain a cleaned_metadata.yaml file (produced by cleansumstats)."
        exit 1
    fi
    
    local n_value=""
    local n_field=""
    local alt_fields=""
    
    if [[ "$which_n" == "effectiveN" ]]; then
        n_field="stats_EffectiveN"
        n_value=$(awk -F': ' '$1=="stats_EffectiveN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
        alt_fields="stats_CaseN + stats_ControlN"
        
        # If effectiveN is missing, check if case/control counts are available (can derive effectiveN)
        if [[ -z "$n_value" ]]; then
            local case_n ctrl_n
            case_n=$(awk -F': ' '$1=="stats_CaseN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
            ctrl_n=$(awk -F': ' '$1=="stats_ControlN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
            if [[ -n "$case_n" ]] && [[ -n "$ctrl_n" ]] && [[ "$case_n" =~ ^[0-9]+$ ]] && [[ "$ctrl_n" =~ ^[0-9]+$ ]]; then
                # Can derive effectiveN from case/control
                n_value="derivable"
            fi
        fi
    else
        # totalN (default)
        n_field="stats_TotalN"
        n_value=$(awk -F': ' '$1=="stats_TotalN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
        alt_fields="stats_CaseN + stats_ControlN (for effectiveN)"
    fi
    
    if [[ -z "$n_value" ]]; then
        log_error "Required sample size field is missing or empty in metadata."
        log_error ""
        log_error "  Config:    whichn: ${which_n}"
        log_error "  Expected:  ${n_field} (or ${alt_fields})"
        log_error "  Metadata:  ${metadata_file}"
        log_error ""
        log_error "The posteriors calculation requires sample size (N) to be present."
        log_error "Please ensure the cleansumstats output includes this field, or"
        log_error "manually add '${n_field}: <value>' to the metadata file."
        exit 1
    fi
    
    log_debug "Metadata validated: ${n_field} = ${n_value}"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_filter_variants() {
    local sumstat_name="$1"
    local specific_chr="${2:-}"  # Optional: run only specific chromosome
    
    log_step "Running filter-variants for: $sumstat_name"
    
    # Check dependencies
    check_filter_variants_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "formatted"
    migrate_sumstat_step_dir "$sumstat_dir" "filtered"
    local format_dir
    format_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    ensure_dir "$step_dir"
    
    # Early validation: check that metadata has required N fields (safety net for non-sbatch runs)
    local input_dir="${CFG_INPUT}"
    local metadata_file="${input_dir}/cleaned_metadata.yaml"
    local which_n="${CFG_WHICHN:-totalN}"
    validate_sumstat_n_field "$metadata_file" "$which_n"
    
    # Auto-detect single-chromosome runs from config (important for --sbatch-array mode)
    if [[ -z "$specific_chr" ]] && [[ -n "${CFG_CHROMOSOMES:-}" ]] && [[ "${CFG_CHROMOSOMES}" =~ ^(chr)?[0-9]+$ ]]; then
        specific_chr="${CFG_CHROMOSOMES#chr}"
    fi
    
    # Detect input format: per-chromosome files or single file
    local use_perchr_input=false
    if [[ -f "${format_dir}/chr1.tsv" ]]; then
        use_perchr_input=true
    elif [[ -f "${format_dir}/sumstat_formatted.tsv.gz" ]]; then
        use_perchr_input=false
    else
        log_error "No formatted sumstat found. Run 'pgscalculator format-sumstat' first"
        log_error "Expected: ${format_dir}/chr*.tsv or ${format_dir}/sumstat_formatted.tsv.gz"
        exit 1
    fi
    
    # Per-chromosome mode
    if [[ -n "$specific_chr" ]]; then
        run_filter_variants_chr "$sumstat_name" "$specific_chr"
        return $?
    fi
    
    # Full mode: check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local input_dir="${CFG_INPUT}"
    local metadata_file="${input_dir}/cleaned_metadata.yaml"
    local which_n="${CFG_WHICHN:-totalN}"
    local prep_mapfile="${prep_dir}/variant_map.tsv"
    local sumstat_mapfile="${sumstat_dir}/variant_map.tsv"
    local sumstat_for_posteriors="${sumstat_dir}/sumstat_for_posteriors.tsv.gz"
    local list_gt="${CFG_FILTERS_INCLUSION_LIST_GT:-}"
    local list_ss="${CFG_FILTERS_INCLUSION_LIST_SS:-}"
    local list_ld="${CFG_FILTERS_INCLUSION_LIST_LD:-}"

    # Normalize user inclusion lists
    for _list_var in list_gt list_ss list_ld; do
        local _val="${!_list_var}"
        if [[ -n "$_val" ]] && [[ "${_val,,}" == "false" ]]; then
            printf -v "$_list_var" ""
        elif [[ -n "$_val" ]] && [[ ! -f "$_val" ]]; then
            log_error "User inclusion list not found: ${_val}"
            exit 1
        fi
    done
    
    if [[ "$use_perchr_input" == true ]]; then
        # New mode: process per-chromosome files
        log_info "Processing per-chromosome formatted files"
        
        local total_input=0
        local total_output=0
        local success_count=0
        local fail_count=0
        
        for chr in $(get_chromosomes); do
            local chr_input="${format_dir}/chr${chr}.tsv"
            if [[ ! -f "$chr_input" ]]; then
                log_debug "chr${chr}: no formatted input file, skipping"
                continue
            fi
            
            if run_filter_variants_chr "$sumstat_name" "$chr"; then
                ((success_count++))
                local chr_out="${step_dir}/chr${chr}_filtered.tsv"
                if [[ -f "$chr_out" ]]; then
                    local n
                    n=$(wc -l < "$chr_out")
                    total_output=$((total_output + n - 1))
                fi
            else
                ((fail_count++))
            fi
            
            local n_in
            n_in=$(wc -l < "$chr_input")
            total_input=$((total_input + n_in - 1))
        done
        
        # Concatenate per-chromosome mapfiles into single sumstat mapfile
        log_substep "Concatenating per-chromosome mapfiles"
        concatenate_chr_mapfiles "$step_dir" "$sumstat_mapfile"
        
        # Create combined sumstat_for_posteriors.tsv.gz from per-chr filtered files
        log_substep "Creating combined filtered sumstat"
        concatenate_chr_filtered "$step_dir" "$sumstat_for_posteriors"
        
        if [[ $fail_count -gt 0 ]]; then
            log_warn "filter-variants had issues for ${fail_count} chromosome(s)"
        fi
        
        log_info "Processed ${success_count} chromosomes: ${total_input} input -> ${total_output} output variants"
        
    else
        # Legacy mode: single formatted file
        local formatted_sumstat="${format_dir}/sumstat_formatted.tsv.gz"
    
    # Count input variants
    local input_count
    input_count=$(zcat "$formatted_sumstat" | wc -l)
    input_count=$((input_count - 1))
    log_info "Input variants: ${input_count}"
    
        # Step 1: Build sumstat-annotated mapfile + reduce sumstat to mapfile intersection
        log_substep "Building sumstat mapfile and reducing to mapfile intersection"
        build_sumstat_map_and_reduce "$prep_mapfile" "$formatted_sumstat" "$sumstat_mapfile" "$sumstat_for_posteriors" "$list_gt" "$list_ss" "$list_ld"
        if [[ ! -s "$sumstat_for_posteriors" ]]; then
            log_error "sumstat_for_posteriors is missing or empty: ${sumstat_for_posteriors}"
            exit 1
        fi
    
    local filtered_count
        filtered_count=$(zcat "$sumstat_for_posteriors" | wc -l)
    filtered_count=$((filtered_count - 1))
    local reduction_pct
    reduction_pct=$(awk "BEGIN {printf \"%.1f\", (1 - $filtered_count / $input_count) * 100}")
    log_info "After inclusion list filter: ${filtered_count} variants (${reduction_pct}% reduction)"
    
        # If the whole-file filtering yields 0 variants, hard-exit.
        if [[ "$filtered_count" -le 0 ]]; then
            log_error "No variants left after mapfile reduction/inclusion-list filtering (0 variants)."
            log_error "Hard exiting: there is nothing to process in downstream steps."
            exit 1
        fi
        
        # Step 2: Derive N/EAF/B/SE on the filtered subset
    log_substep "Deriving N/EAF/B/SE statistics"
        local tmp_filtered
        tmp_filtered="$(mktemp "${step_dir}/sumstat_filtered.tsv.tmp.XXXXXX")"
        zcat "$sumstat_for_posteriors" > "${step_dir}/sumstat_filtered_raw.tsv"
        derive_stats "${step_dir}/sumstat_filtered_raw.tsv" "$tmp_filtered" "$metadata_file" "$which_n" "$prep_dir" "$step_dir"
        mv -f "$tmp_filtered" "${step_dir}/sumstat_filtered.tsv"

        # If derivation yields 0 variants, hard-exit.
        local n_lines
        n_lines=$(wc -l < "${step_dir}/sumstat_filtered.tsv" | awk '{print $1}')
        if [[ "$n_lines" -le 1 ]]; then
            log_error "Filtered sumstat derivation produced 0 variants (lines=${n_lines})."
            exit 1
        fi
        
    rm -f "${step_dir}/sumstat_filtered_raw.tsv"
    
    # Step 3: Split filtered sumstat by chromosome
    log_substep "Splitting filtered sumstat by chromosome"
    split_filtered_by_chr "${step_dir}/sumstat_filtered.tsv" "$step_dir"
    
    # Step 4: Compress the main filtered file
    gzip -f "${step_dir}/sumstat_filtered.tsv"
    fi
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local output_count=0
    if [[ -f "${step_dir}/sumstat_filtered.tsv.gz" ]]; then
    output_count=$(zcat "${step_dir}/sumstat_filtered.tsv.gz" | wc -l)
    output_count=$((output_count - 1))
    elif [[ -f "$sumstat_for_posteriors" ]]; then
        output_count=$(zcat "$sumstat_for_posteriors" | wc -l)
        output_count=$((output_count - 1))
    fi
    
    log_info "Output variants: ${output_count}"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PER-CHROMOSOME PROCESSING
# =============================================================================

run_filter_variants_chr() {
    local sumstat_name="$1"
    local chr="$2"
    
    log_substep "Processing chromosome ${chr}"
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local format_dir
    format_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    ensure_dir "$step_dir"
    
    local chr_input="${format_dir}/chr${chr}.tsv"
    local chr_output="${step_dir}/chr${chr}_filtered.tsv"
    local chr_mapfile="${step_dir}/chr${chr}_map.tsv"
    
    # Check if already processed
    if [[ -f "$chr_output" ]] && [[ $(wc -l < "$chr_output") -gt 1 ]]; then
        log_debug "chr${chr}: already processed, skipping"
        return 0
    fi
    
    # Check input exists
    if [[ ! -f "$chr_input" ]]; then
        log_warn "chr${chr}: no formatted input file: ${chr_input}"
        # Write placeholder
        echo "CHR	POS	RSID	EffectAllele	OtherAllele	EAF	B	SE	P	N	LDREF_SNPID" > "$chr_output"
        echo "formatted_missing" > "${step_dir}/FAILED_chr${chr}"
        return 1
    fi
    
    local input_dir="${CFG_INPUT}"
    local metadata_file="${input_dir}/cleaned_metadata.yaml"
    local which_n="${CFG_WHICHN:-totalN}"
    local prep_mapfile="${prep_dir}/variant_map.tsv"
    local list_gt="${CFG_FILTERS_INCLUSION_LIST_GT:-}"
    local list_ss="${CFG_FILTERS_INCLUSION_LIST_SS:-}"
    local list_ld="${CFG_FILTERS_INCLUSION_LIST_LD:-}"
    
    # Normalize user inclusion lists
    for _list_var in list_gt list_ss list_ld; do
        local _val="${!_list_var}"
        if [[ -n "$_val" ]] && [[ "${_val,,}" == "false" ]]; then
            printf -v "$_list_var" ""
        fi
    done
    
    local tmpdir
    tmpdir=$(make_tmpdir "filter_variants_chr${chr}")
    local reduced_tmp="${tmpdir}/reduced.tsv"
    local map_tmp="${tmpdir}/map.tsv"
    
    # Count input
    local input_count
    input_count=$(wc -l < "$chr_input")
    input_count=$((input_count - 1))
    log_debug "chr${chr}: ${input_count} input variants"
    
    # Build mapfile and reduce for this chromosome only
    # Filter prep_mapfile to just this chromosome for efficiency
    build_sumstat_map_and_reduce_chr "$prep_mapfile" "$chr_input" "$map_tmp" "$reduced_tmp" "$chr" "$list_gt" "$list_ss" "$list_ld"
    
    if [[ ! -s "$reduced_tmp" ]] || [[ $(wc -l < "$reduced_tmp") -le 1 ]]; then
        log_warn "chr${chr}: no variants after mapfile reduction"
        echo "CHR	POS	RSID	EffectAllele	OtherAllele	EAF	B	SE	P	N	LDREF_SNPID" > "$chr_output"
        echo "no_variants_after_reduction" > "${step_dir}/FAILED_chr${chr}"
        rm -rf "$tmpdir"
        return 1
    fi
    
    # Derive stats
    local derived_tmp="${tmpdir}/derived.tsv"
    derive_stats "$reduced_tmp" "$derived_tmp" "$metadata_file" "$which_n" "$prep_dir" ""
    
    if [[ ! -s "$derived_tmp" ]] || [[ $(wc -l < "$derived_tmp") -le 1 ]]; then
        log_warn "chr${chr}: no variants after stat derivation"
        echo "CHR	POS	RSID	EffectAllele	OtherAllele	EAF	B	SE	P	N	LDREF_SNPID" > "$chr_output"
        echo "no_variants_after_derivation" > "${step_dir}/FAILED_chr${chr}"
        rm -rf "$tmpdir"
        return 1
    fi
    
    # Move outputs to final locations
    mv -f "$derived_tmp" "$chr_output"
    mv -f "$map_tmp" "$chr_mapfile"
    rm -rf "$tmpdir"
    
    # Report
    local output_count
    output_count=$(wc -l < "$chr_output")
    output_count=$((output_count - 1))
    log_debug "chr${chr}: ${output_count} output variants"
    
    # Mark chr as completed
    date '+%Y-%m-%d %H:%M:%S' > "${step_dir}/.completed_chr${chr}"
    
    return 0
}

build_sumstat_map_and_reduce_chr() {
    local prep_mapfile="$1"
    local chr_input="$2"
    local out_map="$3"
    local out_sumstat="$4"
    local chr="$5"
    local list_gt="${6:-}"
    local list_ss="${7:-}"
    local list_ld="${8:-}"
    
    # Same logic as build_sumstat_map_and_reduce but for a single chromosome
    # and reading from uncompressed TSV instead of gzipped
    awk -F'\t' -v OFS='\t' \
        -v out_map="$out_map" -v out_sumstat="$out_sumstat" \
        -v target_chr="$chr" \
        -v list_gt="$list_gt" -v list_ss="$list_ss" -v list_ld="$list_ld" '
        BEGIN {
            c["A"]="T"; c["T"]="A"; c["C"]="G"; c["G"]="C"
            if (list_gt != "" && list_gt != "NA") {
                while ((getline < list_gt) > 0) {
                    if ($1 != "") gt[$1]=1
                }
                close(list_gt)
            }
            if (list_ss != "" && list_ss != "NA") {
                while ((getline < list_ss) > 0) {
                    if ($1 != "") ss[$1]=1
                }
                close(list_ss)
            }
            if (list_ld != "" && list_ld != "NA") {
                while ((getline < list_ld) > 0) {
                    if ($1 != "") ld[$1]=1
                }
                close(list_ld)
            }
        }
        NR==FNR {
            if (NR==1) {
                for (i=1;i<=NF;i++) {
                    if ($i=="chr") chr_i=i
                    else if ($i=="pos_b37") pos_b37_i=i
                    else if ($i=="pos_b38") pos_b38_i=i
                    else if ($i=="geno_snpid") geno_id_i=i
                    else if ($i=="geno_a1") geno_a1_i=i
                    else if ($i=="geno_a2") geno_a2_i=i
                    else if ($i=="ldref_snpid") ld_id_i=i
                    else if ($i=="ldref_a1") ld_a1_i=i
                    else if ($i=="ldref_a2") ld_a2_i=i
                    else if ($i=="ldref_a2freq") ld_freq_i=i
                }
                next
            }
            # Only load mapfile entries for target chromosome
            if ($chr_i != target_chr) next
            idx++
            chr_arr[idx]=$chr_i
            pos_b37[idx]=$pos_b37_i
            pos_b38[idx]=$pos_b38_i
            geno_id[idx]=$geno_id_i; geno_a1[idx]=toupper($geno_a1_i); geno_a2[idx]=toupper($geno_a2_i)
            ld_id[idx]=$ld_id_i; ld_a1[idx]=toupper($ld_a1_i); ld_a2[idx]=toupper($ld_a2_i)
            ld_freq[idx]=$ld_freq_i
            a1 = (geno_a1[idx]!="NA" ? geno_a1[idx] : ld_a1[idx])
            a2 = (geno_a2[idx]!="NA" ? geno_a2[idx] : ld_a2[idx])
            # Match on pos_b38 (sumstat is GRCh38)
            if (a1!="NA" && a2!="NA" && chr_arr[idx]!="NA" && pos_b38[idx]!="NA") {
                key = chr_arr[idx] ":" pos_b38[idx] ":" a1 ":" a2
                map_idx[key]=idx
            }
            next
        }
        FNR==1 {
            for (i=1;i<=NF;i++) {
                if (!chr_c && ($i=="CHR" || $i=="chr" || $i=="#CHR")) chr_c=i
                else if (!pos_c && ($i=="POS" || $i=="pos" || $i=="BP" || $i=="Position")) pos_c=i
                else if ($i=="RSID" || $i=="rsid" || $i=="SNP" || $i=="ID") snp_c=i
                else if ($i=="EffectAllele" || $i=="A1" || $i=="effect_allele") a1_c=i
                else if ($i=="OtherAllele" || $i=="A2" || $i=="other_allele") a2_c=i
                else if ($i=="EAF") eaf_c=i
            }
            if (!eaf_c) {
                eaf_c = NF + 1
                header_extra="EAF"
            }
            out = $1
            for (i=2;i<=NF;i++) out = out OFS $i
            if (header_extra!="") out = out OFS header_extra
            out = out OFS "LDREF_SNPID"
            print out > out_sumstat
            next
        }
        {
            chr_v = $chr_c; pos_v = $pos_c
            snp_v = $snp_c
            a1_v = toupper($a1_c); a2_v = toupper($a2_c)
            key1 = chr_v ":" pos_v ":" a1_v ":" a2_v
            key2 = chr_v ":" pos_v ":" a2_v ":" a1_v
            fa1 = c[a1_v]; fa2 = c[a2_v]
            key3 = chr_v ":" pos_v ":" fa1 ":" fa2
            key4 = chr_v ":" pos_v ":" fa2 ":" fa1
            matched_idx = (key1 in map_idx) ? map_idx[key1] : ((key2 in map_idx) ? map_idx[key2] : ((key3 in map_idx) ? map_idx[key3] : ((key4 in map_idx) ? map_idx[key4] : 0)))
            if (matched_idx==0) next
            if (list_gt != "" && list_gt != "NA") {
                if (!(geno_id[matched_idx] in gt)) next
            }
            if (list_ss != "" && list_ss != "NA") {
                if (!(snp_v in ss)) next
            }
            if (list_ld != "" && list_ld != "NA") {
                if (!(ld_id[matched_idx] in ld)) next
            }
            # Fill EAF if missing
            eaf_v = (eaf_c <= NF ? $eaf_c : "NA")
            if (eaf_v=="" || eaf_v=="NA") {
                freq = ld_freq[matched_idx]
                if (freq != "" && freq != "NA" && ld_a1[matched_idx] != "NA" && ld_a2[matched_idx] != "NA") {
                    if (a1_v == ld_a2[matched_idx]) eaf_v = freq
                    else if (a1_v == ld_a1[matched_idx]) eaf_v = 1 - freq
                    else eaf_v = freq
                }
            }
            # write sumstat row
            if (header_extra == "") {
                $eaf_c = eaf_v
            }
            out = $1
            for (i=2;i<=NF;i++) out = out OFS $i
            if (header_extra!="") out = out OFS eaf_v
            out = out OFS ld_id[matched_idx]
            print out > out_sumstat
            # record sumstat columns for mapfile
            if (!(matched_idx in sum_snp)) {
                sum_snp[matched_idx]=snp_v
                sum_a1[matched_idx]=a1_v
                sum_a2[matched_idx]=a2_v
            }
        }
        END {
            # Output mapfile with dual positions (pos_b37 and pos_b38)
            print "chr\tpos_b37\tpos_b38\tsumstat_snpid\tsumstat_effect\tsumstat_other\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq" > out_map
            for (i=1;i<=idx;i++) {
                s_snp = (i in sum_snp) ? sum_snp[i] : "NA"
                s_a1 = (i in sum_a1) ? sum_a1[i] : "NA"
                s_a2 = (i in sum_a2) ? sum_a2[i] : "NA"
                print chr_arr[i], pos_b37[i], pos_b38[i], s_snp, s_a1, s_a2, geno_id[i], geno_a1[i], geno_a2[i], ld_id[i], ld_a1[i], ld_a2[i], ld_freq[i] >> out_map
            }
        }
    ' "$prep_mapfile" "$chr_input"
}

concatenate_chr_mapfiles() {
    local step_dir="$1"
    local output="$2"
    
    # Concatenate per-chromosome mapfiles
    local first=true
    for chr in $(get_chromosomes); do
        local chr_map="${step_dir}/chr${chr}_map.tsv"
        if [[ -f "$chr_map" ]]; then
            if [[ "$first" == true ]]; then
                cat "$chr_map" > "$output"
                first=false
            else
                tail -n +2 "$chr_map" >> "$output"
            fi
        fi
    done
    
    if [[ "$first" == true ]]; then
        # No mapfiles found, write empty header with dual positions
        echo "chr	pos_b37	pos_b38	sumstat_snpid	sumstat_effect	sumstat_other	geno_snpid	geno_a1	geno_a2	ldref_snpid	ldref_a1	ldref_a2	ldref_a2freq" > "$output"
    fi
}

concatenate_chr_filtered() {
    local step_dir="$1"
    local output="$2"
    
    # Concatenate per-chromosome filtered files into gzipped output
    local tmpfile
    tmpfile=$(mktemp)
    local first=true
    
    for chr in $(get_chromosomes); do
        local chr_out="${step_dir}/chr${chr}_filtered.tsv"
        if [[ -f "$chr_out" ]] && [[ $(wc -l < "$chr_out") -gt 1 ]]; then
            if [[ "$first" == true ]]; then
                cat "$chr_out" > "$tmpfile"
                first=false
            else
                tail -n +2 "$chr_out" >> "$tmpfile"
            fi
        fi
    done
    
    if [[ "$first" == true ]]; then
        # No files found, write empty header
        echo "CHR	POS	RSID	EffectAllele	OtherAllele	EAF	B	SE	P	N	LDREF_SNPID" > "$tmpfile"
    fi
    
    gzip -c "$tmpfile" > "$output"
    rm -f "$tmpfile"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

build_sumstat_map_and_reduce() {
    local prep_mapfile="$1"
    local formatted_sumstat="$2"
    local sumstat_mapfile="$3"
    local sumstat_for_posteriors="$4"
    local list_gt="${5:-}"
    local list_ss="${6:-}"
    local list_ld="${7:-}"
    
    local tmpdir
    tmpdir=$(make_tmpdir "sumstat_map_reduce")
    local reduced_tmp="${tmpdir}/sumstat_for_posteriors.tsv"
    local map_tmp="${tmpdir}/variant_map.tsv"
    
    awk -F'\t' -v OFS='\t' \
        -v out_map="$map_tmp" -v out_sumstat="$reduced_tmp" \
        -v list_gt="$list_gt" -v list_ss="$list_ss" -v list_ld="$list_ld" '
        BEGIN {
            c["A"]="T"; c["T"]="A"; c["C"]="G"; c["G"]="C"
            if (list_gt != "" && list_gt != "NA") {
                while ((getline < list_gt) > 0) {
                    if ($1 != "") gt[$1]=1
                }
                close(list_gt)
            }
            if (list_ss != "" && list_ss != "NA") {
                while ((getline < list_ss) > 0) {
                    if ($1 != "") ss[$1]=1
                }
                close(list_ss)
            }
            if (list_ld != "" && list_ld != "NA") {
                while ((getline < list_ld) > 0) {
                    if ($1 != "") ld[$1]=1
                }
                close(list_ld)
            }
        }
        NR==FNR {
            if (NR==1) {
                for (i=1;i<=NF;i++) {
                    if ($i=="chr") chr_i=i
                    else if ($i=="pos_b37") pos_b37_i=i
                    else if ($i=="pos_b38") pos_b38_i=i
                    else if ($i=="geno_snpid") geno_id_i=i
                    else if ($i=="geno_a1") geno_a1_i=i
                    else if ($i=="geno_a2") geno_a2_i=i
                    else if ($i=="ldref_snpid") ld_id_i=i
                    else if ($i=="ldref_a1") ld_a1_i=i
                    else if ($i=="ldref_a2") ld_a2_i=i
                    else if ($i=="ldref_a2freq") ld_freq_i=i
                }
                next
            }
            idx++
            chr[idx]=$chr_i
            pos_b37[idx]=$pos_b37_i
            pos_b38[idx]=$pos_b38_i
            geno_id[idx]=$geno_id_i; geno_a1[idx]=toupper($geno_a1_i); geno_a2[idx]=toupper($geno_a2_i)
            ld_id[idx]=$ld_id_i; ld_a1[idx]=toupper($ld_a1_i); ld_a2[idx]=toupper($ld_a2_i)
            ld_freq[idx]=$ld_freq_i
            a1 = (geno_a1[idx]!="NA" ? geno_a1[idx] : ld_a1[idx])
            a2 = (geno_a2[idx]!="NA" ? geno_a2[idx] : ld_a2[idx])
            # Match on pos_b38 (sumstat is GRCh38)
            if (a1!="NA" && a2!="NA" && chr[idx]!="NA" && pos_b38[idx]!="NA") {
                key = chr[idx] ":" pos_b38[idx] ":" a1 ":" a2
                map_idx[key]=idx
            }
            next
        }
        FNR==1 {
            for (i=1;i<=NF;i++) {
                if ($i=="CHR" || $i=="chr" || $i=="#CHR") chr_c=i
                else if ($i=="POS" || $i=="pos" || $i=="BP" || $i=="Position") pos_c=i
                else if ($i=="RSID" || $i=="rsid" || $i=="SNP" || $i=="ID") snp_c=i
                else if ($i=="EffectAllele" || $i=="A1" || $i=="effect_allele") a1_c=i
                else if ($i=="OtherAllele" || $i=="A2" || $i=="other_allele") a2_c=i
                else if ($i=="EAF") eaf_c=i
            }
            if (!eaf_c) {
                eaf_c = NF + 1
                header_extra="EAF"
            }
            # Write header with appended LDREF_SNPID
            out = $1
            for (i=2;i<=NF;i++) out = out OFS $i
            if (header_extra!="") out = out OFS header_extra
            out = out OFS "LDREF_SNPID"
            print out > out_sumstat
            next
        }
        {
            chr_v = $chr_c; pos_v = $pos_c
            snp_v = $snp_c
            a1_v = toupper($a1_c); a2_v = toupper($a2_c)
            key1 = chr_v ":" pos_v ":" a1_v ":" a2_v
            key2 = chr_v ":" pos_v ":" a2_v ":" a1_v
            fa1 = c[a1_v]; fa2 = c[a2_v]
            key3 = chr_v ":" pos_v ":" fa1 ":" fa2
            key4 = chr_v ":" pos_v ":" fa2 ":" fa1
            idx = (key1 in map_idx) ? map_idx[key1] : ((key2 in map_idx) ? map_idx[key2] : ((key3 in map_idx) ? map_idx[key3] : ((key4 in map_idx) ? map_idx[key4] : 0)))
            if (idx==0) next
            if (list_gt != "" && list_gt != "NA") {
                if (!(geno_id[idx] in gt)) next
            }
            if (list_ss != "" && list_ss != "NA") {
                if (!(snp_v in ss)) next
            }
            if (list_ld != "" && list_ld != "NA") {
                if (!(ld_id[idx] in ld)) next
            }
            # Fill EAF if missing
            eaf_v = (eaf_c <= NF ? $eaf_c : "NA")
            if (eaf_v=="" || eaf_v=="NA") {
                freq = ld_freq[idx]
                if (freq != "" && freq != "NA" && ld_a1[idx] != "NA" && ld_a2[idx] != "NA") {
                    if (a1_v == ld_a2[idx]) eaf_v = freq
                    else if (a1_v == ld_a1[idx]) eaf_v = 1 - freq
                    else eaf_v = freq
                }
            }
            # write sumstat row
            if (header_extra == "") {
                $eaf_c = eaf_v
            }
            out = $1
            for (i=2;i<=NF;i++) out = out OFS $i
            if (header_extra!="") out = out OFS eaf_v
            out = out OFS ld_id[idx]
            print out > out_sumstat
            # record sumstat columns for mapfile
            if (!(idx in sum_snp)) {
                sum_snp[idx]=snp_v
                sum_a1[idx]=a1_v
                sum_a2[idx]=a2_v
            }
        }
        END {
            # Output mapfile with dual positions (pos_b37 and pos_b38)
            print "chr\tpos_b37\tpos_b38\tsumstat_snpid\tsumstat_effect\tsumstat_other\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq" > out_map
            for (i=1;i<=idx;i++) {
                s_snp = (i in sum_snp) ? sum_snp[i] : "NA"
                s_a1 = (i in sum_a1) ? sum_a1[i] : "NA"
                s_a2 = (i in sum_a2) ? sum_a2[i] : "NA"
                print chr[i], pos_b37[i], pos_b38[i], s_snp, s_a1, s_a2, geno_id[i], geno_a1[i], geno_a2[i], ld_id[i], ld_a1[i], ld_a2[i], ld_freq[i] >> out_map
            }
        }
    ' "$prep_mapfile" <(zcat "$formatted_sumstat")
    
    mv -f "$map_tmp" "$sumstat_mapfile"
    gzip -c "$reduced_tmp" > "$sumstat_for_posteriors"
    rm -rf "$tmpdir"
}

derive_stats() {
    local input="$1"
    local output="$2"
    local metadata_file="$3"
    local which_n="$4"
    local prep_dir="$5"
    local step_dir="$6"
    
    local tmpdir
    tmpdir=$(make_tmpdir "filter_variants")

    # We now persist a per-substep audit file so users can see where variants disappear.
    # This is intentionally in the sumstat's filtered/ folder (not tmp), so it survives failures.
    local audit_file=""
    if [[ -n "${step_dir:-}" ]]; then
        audit_file="${step_dir}/filter_variants_steps.tsv"
        {
            echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        } > "$audit_file"
    fi

    # Also persist a compact "modifications" report: how many rows were *modified* by non-filter steps.
    # (e.g., N filled, EAF filled, B/SE derived). This complements removed_lines.tsv.gz.
    local changes_file=""
    if [[ -n "${step_dir:-}" ]]; then
        changes_file="${step_dir}/filter_variants_changes.tsv"
        {
            echo -e "STEP\tFIELD\tN_CHANGED\tDESC"
        } > "$changes_file"
    fi

    count_variants_tsv() {
        local f="$1"
        if [[ ! -s "$f" ]]; then
            echo 0
            return
        fi
        local n
        n=$(wc -l < "$f" | awk '{print $1}')
        if [[ "$n" -le 1 ]]; then
            echo 0
        else
            echo $((n - 1))
        fi
    }

    # Count missing values for one column name (best-effort; returns 0 if col missing)
    count_missing_col() {
        local f="$1"
        local colname="$2"
        awk -F'\t' -v colname="$colname" '
            NR==1{
                for(i=1;i<=NF;i++) if($i==colname){c=i; break}
                next
            }
            c{
                v=$c
                if(v=="" || v=="NA") m++
            }
            END{ print (m+0) }
        ' "$f"
    }

    # Column index in a TSV header; returns 0 if not found.
    get_col_idx() {
        local f="$1"
        local colname="$2"
        awk -F'\t' -v colname="$colname" '
            NR==1{
                for(i=1;i<=NF;i++) if($i==colname){ print i; exit }
                print 0
            }
        ' "$f"
    }

    # Count rows where a field was missing before and non-missing after.
    # Handles "missing column" in before/after by treating it as always-missing / always-missing.
    count_filled_missing() {
        local before="$1"
        local after="$2"
        local colname="$3"
        local ib ia
        ib=$(get_col_idx "$before" "$colname")
        ia=$(get_col_idx "$after" "$colname")
        awk -F'\t' -v ib="$ib" -v ia="$ia" '
            function is_missing(v) { return (v=="" || v=="NA") }
            FNR==1{next} # skip header
            NR==FNR{
                # before
                if (ib==0) miss[FNR]=1
                else miss[FNR]=is_missing($ib)
                next
            }
            {
                # after (second file)
                if (ia==0) { next }
                if (miss[FNR] && !is_missing($ia)) c++
            }
            END{ print (c+0) }
        ' "$before" "$after"
    }

    audit_row() {
        local step="$1"
        local before="$2"
        local after="$3"
        local desc="$4"
        [[ -n "$audit_file" ]] && echo -e "${step}\t${before}\t${after}\t${desc}" >> "$audit_file"
    }

    change_row() {
        local step="$1"
        local field="$2"
        local n_changed="$3"
        local desc="$4"
        [[ -n "$changes_file" ]] && echo -e "${step}\t${field}\t${n_changed}\t${desc}" >> "$changes_file"
    }

    # Track removed lines during filtering so we can quickly see what removed too much.
    # Format: LINE<TAB>CHR:POS:EA:OA<TAB>REASON
    # One file; reason is tagged with pass1/pass2 so we can still attribute which filter removed it.
    local removed_tmp="${tmpdir}/removed_lines.tsv"
    echo -e "LINE\tVARIANT\tREASON" > "$removed_tmp"
    
    # Step 1: Add/fix N (effective or total based on config)
    local n0
    n0=$(count_variants_tsv "$input")
    add_sample_size "$input" "${tmpdir}/step1.tsv" "$metadata_file" "$which_n"
    local n1
    n1=$(count_variants_tsv "${tmpdir}/step1.tsv")
    local n_missing_1
    n_missing_1=$(count_missing_col "${tmpdir}/step1.tsv" "N")
    local n_filled_1
    n_filled_1=$(count_filled_missing "$input" "${tmpdir}/step1.tsv" "N")
    audit_row "add_sample_size" "$n0" "$n1" "ensure N column (missing_N=${n_missing_1})"
    change_row "add_sample_size" "N" "$n_filled_1" "rows where N was filled (before missing -> after present)"
    
    # Step 2: Filter bad values (first pass - remove NA/invalid before derivation)
    filter_bad_values "${tmpdir}/step1.tsv" "${tmpdir}/step2.tsv" "$removed_tmp" "pass1"
    local n2
    n2=$(count_variants_tsv "${tmpdir}/step2.tsv")
    audit_row "filter_bad_values_pass1" "$n1" "$n2" "drop obviously invalid rows before derivation (see removed_lines.tsv.gz; reasons prefixed pass1:)"
    
    # Step 3: Derive B and SE if missing (from Z, N, EAF)
    add_beta_se "${tmpdir}/step2.tsv" "${tmpdir}/step3.tsv"
    local n3
    n3=$(count_variants_tsv "${tmpdir}/step3.tsv")
    local b_missing_4
    local se_missing_4
    b_missing_4=$(count_missing_col "${tmpdir}/step3.tsv" "B")
    se_missing_4=$(count_missing_col "${tmpdir}/step3.tsv" "SE")
    local b_filled_4
    local se_filled_4
    b_filled_4=$(count_filled_missing "${tmpdir}/step2.tsv" "${tmpdir}/step3.tsv" "B")
    se_filled_4=$(count_filled_missing "${tmpdir}/step2.tsv" "${tmpdir}/step3.tsv" "SE")
    audit_row "add_beta_se" "$n2" "$n3" "derive B/SE if missing (missing_B=${b_missing_4}, missing_SE=${se_missing_4})"
    change_row "add_beta_se" "B" "$b_filled_4" "rows where B was derived/filled (before missing -> after present)"
    change_row "add_beta_se" "SE" "$se_filled_4" "rows where SE was derived/filled (before missing -> after present)"
    
    # Step 4: Filter bad values (second pass - ensure derived values are valid)
    filter_bad_values "${tmpdir}/step3.tsv" "$output" "$removed_tmp" "pass2"
    local n4
    n4=$(count_variants_tsv "$output")
    audit_row "filter_bad_values_pass2" "$n3" "$n4" "final validity filter after derivation (see removed_lines.tsv.gz; reasons prefixed pass2:)"
    
    # Persist removal details + counts into the filtered step directory (kept with sumstat outputs).
    if [[ -n "${step_dir:-}" ]]; then
        local removed_out="${step_dir}/removed_lines.tsv.gz"
        local removed_counts="${step_dir}/removed_reason_counts.tsv"
        gzip -c "$removed_tmp" > "$removed_out"
        {
            echo -e "REASON\tN"
            awk -F'\t' '
                NR==1{next}
                { c[$3]++ }
                END{
                    for (r in c) print r "\t" c[r]
                }
            ' "$removed_tmp" | sort -t$'\t' -k2,2nr
        } > "$removed_counts"

        log_debug "Wrote audit: ${audit_file}"
        log_debug "Wrote changes: ${changes_file}"
        log_debug "Wrote removed lines: ${removed_out}"
        log_debug "Wrote removed reason counts: ${removed_counts}"
    fi
    
    # Clean up
    rm -rf "$tmpdir"
    
    local count
    count=$(wc -l < "$output")
    count=$((count - 1))
    log_debug "After stat derivation: ${count} variants"
}

add_sample_size() {
    local input="$1"
    local output="$2"
    local metadata_file="$3"
    local which_n="$4"
    
    # Prefer metadata-derived N when the sumstat file doesn't include N/CaseN/ControlN.
    # cleaned_metadata.yaml is produced by cleansumstats and includes:
    # - stats_TotalN
    # - stats_EffectiveN
    local default_n=""
    if [[ -f "$metadata_file" ]]; then
        case "$which_n" in
            effectiveN)
                default_n=$(awk -F': ' '$1=="stats_EffectiveN"{print $2; exit}' "$metadata_file" | tr -d '[:space:]')
                ;;
            totalN|*)
                default_n=$(awk -F': ' '$1=="stats_TotalN"{print $2; exit}' "$metadata_file" | tr -d '[:space:]')
                ;;
        esac
    fi

    awk -F'\t' -v OFS='\t' -v which_n="$which_n" -v default_n="$default_n" '
        NR == 1 {
            for(i=1; i<=NF; i++) {
                header[i] = $i
                if($i == "N") n_col = i
                if($i == "CaseN") case_col = i
                if($i == "ControlN") ctrl_col = i
            }

            # Ensure an N column exists for downstream steps (e.g., derivations).
            if (!n_col) {
                n_col = NF + 1
                header[n_col] = "N"
                NF = n_col
            }

            # Print (possibly-augmented) header
            out = header[1]
            for (i=2; i<=NF; i++) out = out OFS header[i]
            print out
            next
        }
        {
            if ($n_col != "NA" && $n_col != "") {
                # N already present
                print
            } else if (which_n == "effectiveN" && case_col && ctrl_col) {
                # Calculate effective N: 4 * (cases * controls) / (cases + controls)
                if ($case_col != "NA" && $ctrl_col != "NA" && $case_col > 0 && $ctrl_col > 0) {
                    eff_n = 4 * ($case_col * $ctrl_col) / ($case_col + $ctrl_col)
                        $n_col = eff_n
                }
                print
            } else if (default_n != "" && default_n != "NA" && default_n + 0 > 0) {
                # Fallback: fill N from metadata
                $n_col = default_n
                print
            } else {
                # Keep as is
                print
            }
        }
    ' "$input" > "$output"
}

force_eaf() {
    local input="$1"
    local output="$2"
    local ldref_eaf_file="$3"
    
    # If ldref_eaf.tsv exists, use it as preferred fallback
    # Priority: EAF (from sumstat) > ldref_eaf (from LD reference) > EAF_1KG (last resort)
    if [[ -f "$ldref_eaf_file" ]]; then
        log_debug "Using LD reference EAF as fallback from: $ldref_eaf_file"
        
        # First, find SNP/RSID column in input
        local snp_col
        snp_col=$(head -1 "$input" | awk -F'\t' '{
            for(i=1; i<=NF; i++) {
                if($i == "SNP" || $i == "RSID" || $i == "rsid" || $i == "ID") {
                    print i
                    exit
                }
            }
        }')

        # Fast-path: if EAF exists and has no missing values, don't build an in-memory RSID map.
        # This avoids large awk associative arrays (can OOM on small-memory nodes).
        local has_missing_eaf
        has_missing_eaf=$(
            awk -F'\t' '
                NR==1{
                    for(i=1;i<=NF;i++){
                        if($i=="EAF"){e=i; break}
                    }
                    next
                }
                e{
                    if($e=="" || $e=="NA"){ print 1; exit }
                }
                END{ if(!e) print 1; else if(NR==1) print 0; else print 0 }
            ' "$input"
        )
        if [[ "$has_missing_eaf" -eq 0 ]]; then
            log_debug "EAF present and complete; skipping ldref EAF mapping."
            cp -f "$input" "$output"
            return 0
        fi
        
        awk -F'\t' -v OFS='\t' -v snp_col="$snp_col" '
            # Load LD reference EAF (RSID -> A2Freq), and optionally allele columns for alignment.
            ARGIND == 1 && FNR > 1 {
                # ldref_eaf format: RSID, A1, A2, A2Freq
                ldref_eaf[$1] = $4
                # Store allele columns only if we have an effect allele column to align against.
                if (need_align) {
                    ldref_a1[$1] = $2
                    ldref_a2[$1] = $3
                }
                next
            }
            # Process sumstat
            ARGIND == 2 && FNR == 1 {
                for(i=1; i<=NF; i++) {
                    header[i] = $i
                    if($i == "EAF") eaf_col = i
                    if($i == "EAF_1KG") eaf_1kg_col = i
                    if($i == "A1" || $i == "EffectAllele") a1_col = i
                }
                need_align = (a1_col ? 1 : 0)
                print
                next
            }
            ARGIND == 2 {
                rsid = $snp_col
                
                # If EAF is missing or NA, try to fill from ldref
                if (eaf_col && ($eaf_col == "NA" || $eaf_col == "")) {
                    if (rsid in ldref_eaf) {
                        # Check allele alignment
                        ldref_freq = ldref_eaf[rsid]
                        if (a1_col && $a1_col == ldref_a2[rsid]) {
                            # A1 matches ldref A2, use A2Freq directly
                            $eaf_col = ldref_freq
                        } else if (a1_col && $a1_col == ldref_a1[rsid]) {
                            # A1 matches ldref A1, flip frequency
                            $eaf_col = 1 - ldref_freq
                        } else {
                            # Cannot align, use as-is (assume A2 is effect allele in ldref)
                            $eaf_col = ldref_freq
                        }
                    } else if (eaf_1kg_col && $eaf_1kg_col != "NA" && $eaf_1kg_col != "") {
                        # Last resort: use EAF_1KG
                        $eaf_col = $eaf_1kg_col
                    }
                }
                print
            }
        ' "$ldref_eaf_file" "$input" > "$output"
    else
        # No ldref_eaf file, fall back to old behavior (EAF_1KG)
        log_debug "No LD reference EAF file found, using EAF_1KG as fallback"
    
    awk -F'\t' -v OFS='\t' '
        NR == 1 {
            has_eaf = 0
            has_eaf_1kg = 0
            for(i=1; i<=NF; i++) {
                header[i] = $i
                if($i == "EAF") { eaf_col = i; has_eaf = 1 }
                if($i == "EAF_1KG") { eaf_1kg_col = i; has_eaf_1kg = 1 }
            }
            print
            next
        }
        {
            if (has_eaf && ($eaf_col == "NA" || $eaf_col == "") && has_eaf_1kg) {
                # Use EAF_1KG as fallback
                $eaf_col = $eaf_1kg_col
            }
            print
        }
    ' "$input" > "$output"
    fi
}

filter_bad_values() {
    local input="$1"
    local output="$2"
    local removed_file="${3:-}"
    local pass_tag="${4:-filter}"
    
    awk -F'\t' -v OFS='\t' -v removed_file="$removed_file" -v pass_tag="$pass_tag" '
        NR == 1 {
            for(i=1; i<=NF; i++) {
                header[i] = $i
                if($i == "B" || $i == "BETA") b_col = i
                if($i == "SE") se_col = i
                if($i == "EAF") eaf_col = i
                if($i == "CHR" || $i == "chr" || $i == "#CHR") chr_col = i
                if($i == "POS" || $i == "BP" || $i == "pos" || $i == "Position") pos_col = i
                if($i == "EffectAllele" || $i == "A1") ea_col = i
                if($i == "OtherAllele" || $i == "A2") oa_col = i
            }
            print
            next
        }
        {
            valid = 1
            reason = ""

            # Build identifier string: chr:pos:EA:OA (best-effort)
            chr = (chr_col ? $chr_col : "NA")
            pos = (pos_col ? $pos_col : "NA")
            ea  = (ea_col  ? $ea_col  : "NA")
            oa  = (oa_col  ? $oa_col  : "NA")
            vid = chr ":" pos ":" ea ":" oa
            
            # Check B/BETA
            if (b_col) {
                # Match v1 behavior: treat B==0 as invalid for downstream posterior models.
                if ($b_col == "NA" || $b_col == "") { valid = 0; reason = "beta_missing" }
                else if (($b_col + 0) == 0) { valid = 0; reason = "beta_zero" }
            }
            
            # Check SE
            if (se_col) {
                if (valid && ($se_col == "NA" || $se_col == "")) { valid = 0; reason = "se_missing" }
                else if (valid && (($se_col + 0) == 0)) { valid = 0; reason = "se_zero" }
            }
            
            # Check EAF (pass1 only)
            if (eaf_col && pass_tag != "pass2") {
                if (valid && ($eaf_col == "NA" || $eaf_col == "")) { valid = 0; reason = "eaf_missing" }
                else if (valid && (($eaf_col + 0) == 0 || ($eaf_col + 0) == 1)) { valid = 0; reason = "eaf_boundary" }
            }
            
            if (valid) {
                print
            } else {
                # Best-effort logging of removed rows (appends). Include pass tag in reason for attribution.
                if (removed_file != "") {
                    print NR, vid, (pass_tag ":" reason) >> removed_file
                }
            }
        }
    ' "$input" > "$output"
}

add_beta_se() {
    local input="$1"
    local output="$2"
    
    # Derive B and SE from Z, N, EAF if missing
    # Formula: denom^2 = 2 * EAF * (1 - EAF) * (N + Z^2)
    #          SE = 1 / sqrt(denom^2)
    #          B = Z / sqrt(denom^2)
    
    awk -F'\t' -v OFS='\t' '
        NR == 1 {
            for(i=1; i<=NF; i++) {
                header[i] = $i
                if($i == "B" || $i == "BETA") b_col = i
                if($i == "SE") se_col = i
                if($i == "Z") z_col = i
                if($i == "N") n_col = i
                if($i == "EAF") eaf_col = i
            }
            # If B/SE columns are missing entirely, append them so downstream steps
            # (e.g. sbayesR) always have consistent columns.
            if (!b_col) {
                b_col = NF + 1
                header[b_col] = "B"
                NF = b_col
            }
            if (!se_col) {
                se_col = NF + 1
                header[se_col] = "SE"
                NF = se_col
            }

            # Print (possibly-augmented) header
            out = header[1]
            for (i=2; i<=NF; i++) out = out OFS header[i]
            print out
            next
        }
        {
            # 1) Universal derivation when SE is missing but B and Z exist:
            #    SE = |B / Z| (when Z != 0)
            # This covers common cleansumstats outputs that provide B and Z but not SE.
            if (z_col && b_col && ($se_col == "NA" || $se_col == "") ) {
                z = $z_col
                b = $b_col
                if (z != "NA" && b != "NA" && (z + 0) != 0) {
                    derived_se_bz = (b + 0) / (z + 0)
                    if (derived_se_bz < 0) derived_se_bz = -derived_se_bz
                    if (derived_se_bz > 0) {
                        $se_col = derived_se_bz
                    }
                }
            }

            # 2) If B is missing but Z and SE exist:
            #    B = Z * SE
            if (z_col && se_col && ($b_col == "NA" || $b_col == "") ) {
                z = $z_col
                se = $se_col
                if (z != "NA" && se != "NA" && (se + 0) > 0) {
                    $b_col = (z + 0) * (se + 0)
                }
            }

            # 3) Original derivation for cases where B/SE missing but Z, N, EAF available
            if (z_col && n_col && eaf_col) {
                z = $z_col
                n = $n_col
                eaf = $eaf_col
                
                if (z != "NA" && n != "NA" && eaf != "NA" && 
                    n > 0 && eaf > 0 && eaf < 1) {
                    
                    denom2 = 2 * eaf * (1 - eaf) * (n + z * z)
                    if (denom2 > 0) {
                        sqrt_denom2 = sqrt(denom2)
                        derived_se = 1 / sqrt_denom2
                        derived_b = z / sqrt_denom2
                        
                        # Fill in if missing
                        if (b_col && ($b_col == "NA" || $b_col == "")) {
                            $b_col = derived_b
                        }
                        if (se_col && ($se_col == "NA" || $se_col == "")) {
                            $se_col = derived_se
                        }
                    }
                }
            }
            print
        }
    ' "$input" > "$output"
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
