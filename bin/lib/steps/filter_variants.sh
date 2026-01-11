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
    require_command "join" "join is required for file merging"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    
    # Check that prep-inclusion-list has been run
    require_file "${prep_dir}/inclusion_list/variant_inclusion_list.tsv" "Run 'pgscalculator prep-inclusion-list' first"
    require_file "${prep_dir}/inclusion_list/.rsid_index" "Run 'pgscalculator prep-inclusion-list' first"
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
    migrate_sumstat_step_dir "$sumstat_dir" "formatted"
    migrate_sumstat_step_dir "$sumstat_dir" "filtered"
    local format_dir
    format_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    ensure_dir "$step_dir"
    
    # Check that format-sumstat has been run
    require_file "${format_dir}/sumstat_formatted.tsv.gz" "Run 'pgscalculator format-sumstat' first"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local inclusion_file="${prep_dir}/inclusion_list/.rsid_index"
    local formatted_sumstat="${format_dir}/sumstat_formatted.tsv.gz"
    local input_dir="${CFG_INPUT}"
    local metadata_file="${input_dir}/cleaned_metadata.yaml"
    local which_n="${CFG_WHICHN:-totalN}"
    
    # Count input variants
    local input_count
    input_count=$(zcat "$formatted_sumstat" | wc -l)
    input_count=$((input_count - 1))
    log_info "Input variants: ${input_count}"
    
    # Step 1: Filter sumstat to inclusion list variants
    log_substep "Filtering to inclusion list variants"
    filter_to_inclusion_list "$formatted_sumstat" "$inclusion_file" "${step_dir}/sumstat_filtered_raw.tsv"
    
    local filtered_count
    filtered_count=$(wc -l < "${step_dir}/sumstat_filtered_raw.tsv")
    filtered_count=$((filtered_count - 1))
    local reduction_pct
    reduction_pct=$(awk "BEGIN {printf \"%.1f\", (1 - $filtered_count / $input_count) * 100}")
    log_info "After inclusion list filter: ${filtered_count} variants (${reduction_pct}% reduction)"

    # If the whole-file filtering yields 0 variants, hard-exit.
    # (The "continue-on-failure" model is only for chromosome-parallel steps later on.)
    if [[ "$filtered_count" -le 0 ]]; then
        log_error "No variants left after inclusion-list filtering (0 variants)."
        log_error "Hard exiting: there is nothing to process in downstream steps."
        exit 1
    fi
    
    # Step 2: Derive N/EAF/B/SE on the filtered subset (much faster than on full sumstat)
    log_substep "Deriving N/EAF/B/SE statistics"
    derive_stats "${step_dir}/sumstat_filtered_raw.tsv" "${step_dir}/sumstat_filtered.tsv" "$metadata_file" "$which_n" "$prep_dir"
    # If derivation yields 0 variants, hard-exit (whole-file step).
    if [[ ! -s "${step_dir}/sumstat_filtered.tsv" ]]; then
        log_error "Filtered sumstat derivation produced an empty file: ${step_dir}/sumstat_filtered.tsv"
        log_error "Hard exiting: there is nothing to process in downstream steps."
        exit 1
    fi
    
    # Clean up intermediate file
    rm -f "${step_dir}/sumstat_filtered_raw.tsv"
    
    # Step 3: Split filtered sumstat by chromosome
    log_substep "Splitting filtered sumstat by chromosome"
    split_filtered_by_chr "${step_dir}/sumstat_filtered.tsv" "$step_dir"
    
    # Step 4: Compress the main filtered file
    gzip -f "${step_dir}/sumstat_filtered.tsv"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local output_count
    output_count=$(zcat "${step_dir}/sumstat_filtered.tsv.gz" | wc -l)
    output_count=$((output_count - 1))
    
    log_info "Output variants: ${output_count}"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

filter_to_inclusion_list() {
    local input_sumstat="$1"
    local inclusion_file="$2"
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
    
    # Filter using awk (inclusion list is sorted)
    awk -F'\t' -v OFS='\t' -v snp_col="$snp_col" '
        ARGIND == 1 {
            inclusion[$1] = 1
            next
        }
        FNR == 1 {
            print
            next
        }
        {
            if ($snp_col in inclusion) {
                print
            }
        }
    ' "$inclusion_file" <(zcat "$input_sumstat") > "$output_file"
}

derive_stats() {
    local input="$1"
    local output="$2"
    local metadata_file="$3"
    local which_n="$4"
    local prep_dir="$5"
    
    local tmpdir
    tmpdir=$(make_tmpdir "filter_variants")
    
    # Step 1: Add/fix N (effective or total based on config)
    add_sample_size "$input" "${tmpdir}/step1.tsv" "$metadata_file" "$which_n"
    
    # Step 2: Force EAF column (use ldref_eaf as fallback, preferred over EAF_1KG)
    local ldref_eaf_file="${prep_dir}/references/ldref_eaf.tsv"
    force_eaf "${tmpdir}/step1.tsv" "${tmpdir}/step2.tsv" "$ldref_eaf_file"
    
    # Step 3: Filter bad values (first pass - remove NA/invalid before derivation)
    filter_bad_values "${tmpdir}/step2.tsv" "${tmpdir}/step3.tsv"
    
    # Step 4: Derive B and SE if missing (from Z, N, EAF)
    add_beta_se "${tmpdir}/step3.tsv" "${tmpdir}/step4.tsv"
    
    # Step 5: Filter bad values (second pass - ensure derived values are valid)
    filter_bad_values "${tmpdir}/step4.tsv" "$output"
    
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
    
    awk -F'\t' -v OFS='\t' -v which_n="$which_n" '
        NR == 1 {
            for(i=1; i<=NF; i++) {
                header[i] = $i
                if($i == "N") n_col = i
                if($i == "CaseN") case_col = i
                if($i == "ControlN") ctrl_col = i
            }
            print
            next
        }
        {
            if (n_col && $n_col != "NA" && $n_col != "") {
                # N already exists
                print
            } else if (which_n == "effectiveN" && case_col && ctrl_col) {
                # Calculate effective N: 4 * (cases * controls) / (cases + controls)
                if ($case_col != "NA" && $ctrl_col != "NA" && $case_col > 0 && $ctrl_col > 0) {
                    eff_n = 4 * ($case_col * $ctrl_col) / ($case_col + $ctrl_col)
                    if (n_col) {
                        $n_col = eff_n
                    }
                }
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
        
        awk -F'\t' -v OFS='\t' -v snp_col="$snp_col" '
            # Load LD reference EAF (RSID -> A2Freq)
            ARGIND == 1 && FNR > 1 {
                # ldref_eaf format: RSID, A1, A2, A2Freq
                ldref_eaf[$1] = $4
                ldref_a1[$1] = $2
                ldref_a2[$1] = $3
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
    
    awk -F'\t' -v OFS='\t' '
        NR == 1 {
            for(i=1; i<=NF; i++) {
                header[i] = $i
                if($i == "B" || $i == "BETA") b_col = i
                if($i == "SE") se_col = i
                if($i == "EAF") eaf_col = i
            }
            print
            next
        }
        {
            valid = 1
            
            # Check B/BETA
            if (b_col) {
                # Match v1 behavior: treat B==0 as invalid for downstream posterior models.
                if ($b_col == "NA" || $b_col == "" || $b_col == 0) valid = 0
            }
            
            # Check SE
            if (se_col) {
                if ($se_col == "NA" || $se_col == "" || $se_col == 0) valid = 0
            }
            
            # Check EAF
            if (eaf_col) {
                if ($eaf_col == "NA" || $eaf_col == "" || 
                    $eaf_col == 0 || $eaf_col == 1) valid = 0
            }
            
            if (valid) print
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
            # Try to derive B and SE if missing but Z, N, EAF available
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
