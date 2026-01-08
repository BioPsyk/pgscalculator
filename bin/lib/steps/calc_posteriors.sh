#!/bin/bash
# pgscalculator v2 - calc-posteriors step
# Calculate posterior effects using sbayesR

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_calc_posteriors_deps() {
    require_command "gctb" "gctb (sbayesR) is required for posterior calculation"
    require_command "awk" "awk is required for text processing"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR" "LDDIR"
    
    # Check LD reference directory
    require_dir "${CFG_LDDIR}" "LD reference directory not found"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_calc_posteriors() {
    local sumstat_name="$1"
    local specific_chr="${2:-}"  # Optional: run only specific chromosome
    
    log_step "Running calc-posteriors for: $sumstat_name"
    
    # Check dependencies
    check_calc_posteriors_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "filtered"
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors"
    local filter_dir
    filter_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors")
    ensure_dir "$step_dir"
    
    # Check that filter-variants has been run
    require_dir "$filter_dir" "Run 'pgscalculator filter-variants' first"
    
    # Check if already completed (only if not running specific chr)
    if [[ -z "$specific_chr" ]] && check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local lddir="${CFG_LDDIR}"
    local mapfile="${CFG_MAPFILE:-${SCRIPT_DIR}/../assets/sumstats_column_names_map.tsv}"
    
    # Get sbayesR options from config
    local sbayesr_options
    sbayesr_options=$(build_sbayesr_options)
    
    log_info "sbayesR options: $sbayesr_options"
    
    # Determine which chromosomes to process
    local chromosomes
    if [[ -n "$specific_chr" ]]; then
        chromosomes="$specific_chr"
    else
        chromosomes=$(get_chromosomes)
    fi
    
    # Process each chromosome
    local success_count=0
    local fail_count=0
    
    for chr in $chromosomes; do
        local chr_filtered="${filter_dir}/chr${chr}_filtered.tsv"
        
        if [[ ! -f "$chr_filtered" ]]; then
            log_warn "No filtered sumstat for chr${chr}, skipping"
            continue
        fi
        
        # Check if chr already processed
        if [[ -f "${step_dir}/chr${chr}.snpRes" ]]; then
            log_debug "chr${chr} already processed, skipping"
            ((success_count++))
            continue
        fi
        
        log_substep "Processing chromosome ${chr}"
        
        if process_chr_posteriors "$chr" "$chr_filtered" "$lddir" "$step_dir" "$mapfile" "$sbayesr_options"; then
            ((success_count++))
        else
            ((fail_count++))
            log_warn "chr${chr} posterior calculation failed"
        fi
    done
    
    # Mark step as completed only if all chromosomes succeeded
    if [[ -z "$specific_chr" ]] && [[ $fail_count -eq 0 ]]; then
        mark_step_completed "$step_dir"
    fi
    
    log_info "Completed: ${success_count} chromosomes, Failed: ${fail_count}"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

process_chr_posteriors() {
    local chr="$1"
    local chr_filtered="$2"
    local lddir="$3"
    local step_dir="$4"
    local mapfile="$5"
    local sbayesr_options="$6"
    
    local chr_workdir="${step_dir}/work_chr${chr}"
    mkdir -p "$chr_workdir"
    
    # Step 1: Format sumstat for sbayesR
    local sbayesr_input="${chr_workdir}/chr${chr}_sbayesr.ma"
    if ! format_for_sbayesr "$chr_filtered" "$sbayesr_input" "$mapfile"; then
        log_error "chr${chr}: failed to generate sbayesR input (.ma). Fix the filtered sumstat columns and retry."
        return 1
    fi

    if [[ ! -f "$sbayesr_input" ]]; then
        log_error "chr${chr}: sbayesR input file missing after generation: $sbayesr_input"
        return 1
    fi
    
    # Check we have variants
    local variant_count
    variant_count=$(wc -l < "$sbayesr_input")
    variant_count=$((variant_count - 1))
    
    if [[ $variant_count -lt 10 ]]; then
        log_warn "chr${chr}: Only ${variant_count} variants, skipping"
        # Create empty output
        echo "SNP A1 A2 b se pval Freq N effect pj" > "${step_dir}/chr${chr}.snpRes"
        return 0
    fi
    
    log_debug "chr${chr}: ${variant_count} variants for sbayesR"
    
    # Step 2: Find LD reference files
    local ld_bin
    local ld_info
    ld_bin=$(find_ld_file "$lddir" "$chr" "bin")
    ld_info=$(find_ld_file "$lddir" "$chr" "info")
    
    if [[ -z "$ld_bin" ]] || [[ -z "$ld_info" ]]; then
        log_error "Could not find LD reference files for chr${chr}"
        return 1
    fi
    
    # Get LD prefix (remove .bin extension)
    local ld_prefix="${ld_bin%.bin}"
    
    log_debug "LD reference: $ld_prefix"
    
    # Step 3: Run sbayesR
    local out_prefix="${chr_workdir}/chr${chr}"
    
    # Build command
    local cmd="gctb --sbayes R"
    cmd+=" --gwas-summary ${sbayesr_input}"
    cmd+=" --ldm ${ld_prefix}"
    cmd+=" --out ${out_prefix}"
    cmd+=" ${sbayesr_options}"
    
    log_debug "Running: $cmd"
    
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log_info "DRY RUN: Would execute: $cmd"
        return 0
    fi
    
    # Execute sbayesR
    if eval "$cmd" > "${chr_workdir}/sbayesr.log" 2>&1; then
        # Move output to final location
        if [[ -f "${out_prefix}.snpRes" ]]; then
            mv "${out_prefix}.snpRes" "${step_dir}/chr${chr}.snpRes"
            log_debug "chr${chr}: sbayesR completed successfully"
            return 0
        else
            log_error "chr${chr}: sbayesR did not produce output file"
            return 1
        fi
    else
        log_error "chr${chr}: sbayesR failed. Check ${chr_workdir}/sbayesr.log"
        return 1
    fi
}

format_for_sbayesr() {
    local input_file="$1"
    local output_file="$2"
    local mapfile="$3"
    
    # sbayesR .ma format: SNP A1 A2 freq b se p n
    # Use the mapfile to determine column mapping
    
    # Get header from input
    local header
    header=$(head -1 "$input_file")

    # Detect delimiter (we expect TSV, but be defensive)
    local fs
    if echo "$header" | grep -q $'\t'; then
        fs=$'\t'
    else
        fs='[ \t]+'
    fi

    # Find column indices (pick a single best match, deterministically)
    local snp_col a1_col a2_col freq_col beta_col se_col p_col n_col

    # Prefer RSID over SNP/ID if multiple exist
    snp_col=$(echo "$header" | awk -F"$fs" '
        { for(i=1;i<=NF;i++) if($i=="RSID"||$i=="rsid") {print i; exit} }
        { for(i=1;i<=NF;i++) if($i=="SNP"||$i=="ID") {print i; exit} }
    ')
    a1_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="A1"||$i=="EffectAllele"||$i=="effect_allele") {print i; exit}}')
    a2_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="A2"||$i=="OtherAllele"||$i=="other_allele") {print i; exit}}')
    freq_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="EAF"||$i=="Freq"||$i=="freq"||$i=="EAF_1KG") {print i; exit}}')
    beta_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="B"||$i=="BETA"||$i=="beta") {print i; exit}}')
    se_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="SE"||$i=="se") {print i; exit}}')
    p_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="P"||$i=="p"||$i=="pval") {print i; exit}}')
    n_col=$(echo "$header" | awk -F"$fs" '{for(i=1;i<=NF;i++) if($i=="N"||$i=="n"||$i=="Neff") {print i; exit}}')

    log_debug "Column indices: SNP=$snp_col A1=$a1_col A2=$a2_col FREQ=$freq_col BETA=$beta_col SE=$se_col P=$p_col N=$n_col"

    # Fail fast if required columns are missing (prevents malformed .ma and GCTB segfaults)
    local missing_cols=""
    [[ -z "$snp_col" ]] && missing_cols="${missing_cols} SNP/RSID"
    [[ -z "$a1_col" ]] && missing_cols="${missing_cols} A1"
    [[ -z "$a2_col" ]] && missing_cols="${missing_cols} A2"
    [[ -z "$freq_col" ]] && missing_cols="${missing_cols} EAF/Freq"
    [[ -z "$beta_col" ]] && missing_cols="${missing_cols} B"
    [[ -z "$se_col" ]] && missing_cols="${missing_cols} SE"
    [[ -z "$p_col" ]] && missing_cols="${missing_cols} P"
    [[ -z "$n_col" ]] && missing_cols="${missing_cols} N"
    if [[ -n "$missing_cols" ]]; then
        log_error "format_for_sbayesr: missing required columns:${missing_cols}"
        log_error "Header was: $header"
        return 1
    fi

    # Format for sbayesR (space-separated)
    echo "SNP A1 A2 freq b se p n" > "$output_file"

    # Build .ma and drop obviously bad rows
    awk -F"$fs" -v OFS=' ' \
        -v snp="$snp_col" -v a1="$a1_col" -v a2="$a2_col" \
        -v freq="$freq_col" -v beta="$beta_col" -v se="$se_col" \
        -v p="$p_col" -v n="$n_col" '
        NR == 1 { next }
        {
            # Guard against missing/NA fields
            if ($snp == "" || $a1 == "" || $a2 == "" || $freq == "" || $beta == "" || $se == "" || $p == "" || $n == "") next
            if ($snp == "NA" || $a1 == "NA" || $a2 == "NA" || $freq == "NA" || $beta == "NA" || $se == "NA" || $p == "NA" || $n == "NA") next

            # Numeric sanity (avoid inf/nan strings)
            if ($freq ~ /nan|inf/i || $beta ~ /nan|inf/i || $se ~ /nan|inf/i || $p ~ /nan|inf/i || $n ~ /nan|inf/i) next

            print $snp, $a1, $a2, $freq, $beta, $se, $p, $n
        }
    ' "$input_file" >> "$output_file"

    # Validate: every data row must have exactly 8 fields
    local bad_lines
    bad_lines=$(awk 'NR>1 && NF!=8 {c++} END{print c+0}' "$output_file")
    if [[ "$bad_lines" -gt 0 ]]; then
        log_error "format_for_sbayesr: generated malformed .ma (${bad_lines} data lines did not have 8 columns): $output_file"
        log_error "Example bad lines:"
        awk 'NR<=30 && NF!=8 {print "  " $0}' "$output_file" | head -5 >&2
        return 1
    fi

    return 0
}

find_ld_file() {
    local lddir="$1"
    local chr="$2"
    local ext="$3"  # bin or info
    
    # Try different naming patterns
    local patterns=(
        "${lddir}/band_chr${chr}.ldm.sparse.${ext}"
        "${lddir}/*chr${chr}.ldm.sparse.${ext}"
        "${lddir}/*_chr${chr}.${ext}"
        "${lddir}/*_${chr}.${ext}"
    )
    
    for pattern in "${patterns[@]}"; do
        local matches
        matches=$(ls $pattern 2>/dev/null | head -1)
        if [[ -n "$matches" ]] && [[ -f "$matches" ]]; then
            echo "$matches"
            return 0
        fi
    done
    
    # Broader search
    local found
    found=$(find "$lddir" -name "*chr${chr}*.${ext}" -type f 2>/dev/null | head -1)
    if [[ -n "$found" ]]; then
        echo "$found"
        return 0
    fi
    
    echo ""
    return 1
}




