#!/bin/bash
# pgscalculator v2 - format-sumstat step
# Format summary statistics: add build coordinates, derive missing stats

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_format_sumstat_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "zcat" "zcat is required for reading gzipped files"
    require_command "gzip" "gzip is required for compressing output"
    
    # Check config variables
    validate_required_config "CFG" "INPUT" "OUTDIR"
    
    local input_dir="${CFG_INPUT}"
    
    # Check for cleaned sumstat files from cleansumstats
    require_file "${input_dir}/cleaned_GRCh38.gz" "Cleaned sumstat (GRCh38) not found in input directory"
    require_file "${input_dir}/cleaned_GRCh37.gz" "Cleaned sumstat (GRCh37 map) not found in input directory"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_format_sumstat() {
    local sumstat_name="$1"
    
    log_step "Running format-sumstat for: $sumstat_name"
    
    # Check dependencies
    check_format_sumstat_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local step_dir="${sumstat_dir}/formatted"
    ensure_dir "$step_dir"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local input_dir="${CFG_INPUT}"
    local input_grch38="${input_dir}/cleaned_GRCh38.gz"
    local input_grch37="${input_dir}/cleaned_GRCh37.gz"
    local metadata_file="${input_dir}/cleaned_metadata.yaml"
    
    # Step 1: Add GRCh37 coordinates to sumstat
    log_substep "Adding GRCh37 build coordinates"
    add_build_coordinates "$input_grch38" "$input_grch37" "${step_dir}/sumstat_with_b37.tsv.gz"
    
    # Step 2: Split by chromosome
    log_substep "Splitting by chromosome"
    split_sumstat_by_chr "${step_dir}/sumstat_with_b37.tsv.gz" "$step_dir"
    
    # Step 3: Process each chromosome (filter NAs, add N, EAF, B, SE)
    log_substep "Processing per-chromosome files"
    local which_n="${CFG_WHICHN:-totalN}"
    
    for chr in $(get_chromosomes); do
        local chr_file="${step_dir}/chr${chr}_raw.tsv"
        if [[ -f "$chr_file" ]]; then
            process_chr_sumstat "$chr" "$chr_file" "$step_dir" "$metadata_file" "$which_n"
        fi
    done
    
    # Step 4: Concatenate processed files
    log_substep "Concatenating processed files"
    concatenate_processed_sumstats "$step_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local final_count
    final_count=$(zcat "${step_dir}/sumstat_formatted.tsv.gz" | wc -l)
    final_count=$((final_count - 1))
    log_info "Formatted sumstat: ${final_count} variants"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

add_build_coordinates() {
    local input_grch38="$1"
    local input_grch37="$2"
    local output_file="$3"
    
    # GRCh37 file has: CHR, POS, RSID (3 columns)
    # GRCh38 file has: CHR, POS, 0, RSID, EffectAllele, ...
    # After paste (b37 first, then b38):
    #   1: CHR_b37, 2: POS_b37, 3: RSID_b37, 4: CHR_b38, 5: POS_b38, 6+: rest
    # We want: CHR_b37, POS_b37, POS_b38, 0, RSID, ... (drop RSID_b37 and CHR_b38)
    # Use cut -f1-2,5- to skip columns 3 and 4
    paste <(zcat "$input_grch37") <(zcat "$input_grch38") | \
        awk -F'\t' -v OFS='\t' '{
            # Replace empty values with NA
            for(i=1; i<=NF; i++) if($i=="") $i="NA"
            print
        }' | cut -f1-2,5- | gzip -c > "$output_file"
}

split_sumstat_by_chr() {
    local input_file="$1"
    local step_dir="$2"
    
    # Get header
    local header
    header=$(zcat "$input_file" | head -1)
    
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
        log_error "Could not find CHR column in sumstat"
        exit 1
    fi
    
    log_debug "CHR column index: $chr_col"
    
    # Split by chromosome
    zcat "$input_file" | awk -F'\t' -v OFS='\t' -v chr_col="$chr_col" -v outdir="$step_dir" '
        NR == 1 {
            header = $0
            next
        }
        {
            chr = $chr_col
            if (chr >= 1 && chr <= 22) {
                outfile = outdir "/chr" chr "_raw.tsv"
                if (!(chr in seen)) {
                    print header > outfile
                    seen[chr] = 1
                }
                print >> outfile
            }
        }
    '
}

process_chr_sumstat() {
    local chr="$1"
    local input_file="$2"
    local step_dir="$3"
    local metadata_file="$4"
    local which_n="$5"
    
    local output_file="${step_dir}/chr${chr}_processed.tsv"
    
    log_debug "Processing chr${chr}"
    
    # Chain of processing steps
    # 1. Filter NA coordinates
    # 2. Add/fix N
    # 3. Force EAF
    # 4. Filter bad values (round 1)
    # 5. Add B and SE
    # 6. Filter bad values (round 2)
    
    local tmpdir="${step_dir}/tmp_chr${chr}"
    mkdir -p "$tmpdir"
    
    # Step 1: Filter NA coordinates (for b37)
    filter_na_coordinates "$input_file" "${tmpdir}/step1.tsv"
    
    # Step 2: Add N (effective or total based on config)
    add_sample_size "${tmpdir}/step1.tsv" "${tmpdir}/step2.tsv" "$metadata_file" "$which_n"
    
    # Step 3: Force EAF column
    force_eaf "${tmpdir}/step2.tsv" "${tmpdir}/step3.tsv"
    
    # Step 4: Filter bad values (first pass)
    filter_bad_values "${tmpdir}/step3.tsv" "${tmpdir}/step4.tsv"
    
    # Step 5: Add B and SE if missing
    add_beta_se "${tmpdir}/step4.tsv" "${tmpdir}/step5.tsv"
    
    # Step 6: Filter bad values (second pass)
    filter_bad_values "${tmpdir}/step5.tsv" "$output_file"
    
    # Clean up temp files
    rm -rf "$tmpdir"
    
    local count
    count=$(wc -l < "$output_file")
    count=$((count - 1))
    log_debug "chr${chr}: ${count} variants after processing"
}

filter_na_coordinates() {
    local input="$1"
    local output="$2"
    
    # Filter out rows where CHR or POS (for b37) is NA
    awk -F'\t' -v OFS='\t' '
        NR == 1 {
            # Find CHR and POS column indices
            for(i=1; i<=NF; i++) {
                if($i == "CHR") chr_col = i
                if($i == "POS") pos_col = i
            }
            print
            next
        }
        {
            if ($chr_col != "NA" && $chr_col != "" && 
                $pos_col != "NA" && $pos_col != "") {
                print
            }
        }
    ' "$input" > "$output"
}

add_sample_size() {
    local input="$1"
    local output="$2"
    local metadata_file="$3"
    local which_n="$4"
    
    # Check if N column already exists and has values
    # If not, try to get from metadata or use effectiveN calculation
    
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
    
    # Ensure EAF column exists; if not, try to use EAF_1KG
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
}

filter_bad_values() {
    local input="$1"
    local output="$2"
    
    # Filter rows with invalid values
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
            print
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

concatenate_processed_sumstats() {
    local step_dir="$1"
    
    local first=1
    local output="${step_dir}/sumstat_formatted.tsv"
    
    for chr in $(get_chromosomes); do
        local chr_file="${step_dir}/chr${chr}_processed.tsv"
        if [[ -f "$chr_file" ]]; then
            if [[ $first -eq 1 ]]; then
                cat "$chr_file" > "$output"
                first=0
            else
                tail -n +2 "$chr_file" >> "$output"
            fi
        fi
    done
    
    # Compress final output
    gzip -f "$output"
}




