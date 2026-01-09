#!/bin/bash
# pgscalculator v2 - calc-score step
# Calculate PGS using plink2

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_calc_score_deps() {
    require_command "plink2" "plink2 is required for PGS calculation"
    require_command "awk" "awk is required for text processing"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR" "GENODIR" "GENOFILE"
    
    # Check genotype files
    require_file "${CFG_GENOFILE}" "Genotype manifest file not found"
    require_dir "${CFG_GENODIR}" "Genotype directory not found"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_calc_score() {
    local sumstat_name="$1"
    local specific_chr="${2:-}"  # Optional: run only specific chromosome
    
    log_step "Running calc-score for: $sumstat_name"
    
    # Check dependencies
    check_calc_score_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors_mapped"
    migrate_sumstat_step_dir "$sumstat_dir" "scores"
    local posteriors_mapped_dir
    posteriors_mapped_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "scores")
    ensure_dir "$step_dir"
    
    # Check that format-posteriors has been run
    require_dir "$posteriors_mapped_dir" "Run 'pgscalculator format-posteriors' first"
    
    # Auto-detect single-chromosome runs from config (important for --sbatch-array mode where config is rewritten)
    if [[ -z "$specific_chr" ]] && [[ -n "${CFG_CHROMOSOMES:-}" ]] && [[ "${CFG_CHROMOSOMES}" =~ ^(chr)?[0-9]+$ ]]; then
        specific_chr="${CFG_CHROMOSOMES#chr}"
    fi

    # Check if already completed (only if not running specific chr).
    # Be defensive: if the marker exists but no chr*.sscore outputs exist, re-run.
    if [[ -z "$specific_chr" ]] && check_step_completed "$step_dir"; then
        if ls "${step_dir}"/chr*.sscore >/dev/null 2>&1; then
            log_info "Step already completed. Use --force to re-run."
            return 0
        fi
        log_warn "Found ${step_dir}/.completed but no chr*.sscore outputs; re-running calc-score."
    fi
    
    local genodir="${CFG_GENODIR}"
    local genofile="${CFG_GENOFILE}"
    
    # Get score columns from config.
    # Default must match posteriors_mapped header: ID A1 A2 Freq Effect SE PIP
    # plink2 --score expects: <variant_id_col> <allele_col> <score_col>
    local score_columns="${CFG_SCORE_COLUMNS:-1 2 5}"
    
    log_info "Score columns: $score_columns"
    
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
        local posteriors_file="${posteriors_mapped_dir}/chr${chr}.snpRes"
        
        if [[ ! -f "$posteriors_file" ]]; then
            log_debug "No mapped posteriors for chr${chr}, skipping"
            continue
        fi
        
        # Check if chr already processed
        if [[ -f "${step_dir}/chr${chr}.sscore" ]]; then
            log_debug "chr${chr} already scored, skipping"
            ((success_count++))
            continue
        fi
        
        log_substep "Scoring chromosome ${chr}"
        
        if score_chr "$chr" "$posteriors_file" "$genodir" "$genofile" "$step_dir" "$score_columns"; then
            ((success_count++))
        else
            ((fail_count++))
            log_warn "chr${chr} scoring failed"
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

score_chr() {
    local chr="$1"
    local posteriors_file="$2"
    local genodir="$3"
    local genofile="$4"
    local step_dir="$5"
    local score_columns="$6"
    
    local chr_workdir="${step_dir}/work_chr${chr}"
    mkdir -p "$chr_workdir"
    
    # Check if posteriors file has data (more than just header)
    local variant_count
    variant_count=$(wc -l < "$posteriors_file")
    variant_count=$((variant_count - 1))
    
    if [[ $variant_count -lt 1 ]]; then
        log_warn "chr${chr}: No variants in posteriors file, creating empty score"
        create_empty_score "${step_dir}/chr${chr}.sscore"
        return 0
    fi
    
    log_debug "chr${chr}: ${variant_count} variants for scoring"
    
    # Get genotype files for this chromosome
    local pgen pvar psam
    pgen=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "pgen")
    pvar=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "pvar")
    psam=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "psam")
    
    # Fall back to plink1 format if plink2 not found
    local geno_format="plink2"
    if [[ -z "$pgen" ]] || [[ ! -f "$pgen" ]]; then
        local bed bim fam
        bed=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "bed")
        bim=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "bim")
        fam=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "fam")
        
        if [[ -n "$bed" ]] && [[ -f "$bed" ]]; then
            geno_format="plink1"
        else
            log_error "chr${chr}: No genotype files found"
            return 1
        fi
    fi
    
    # Extract variant IDs from posteriors for --extract
    local extract_file="${chr_workdir}/variants.txt"
    awk -F'\t' 'NR > 1 {print $1}' "$posteriors_file" > "$extract_file"
    
    local extract_count
    extract_count=$(wc -l < "$extract_file")
    log_debug "chr${chr}: Extracting ${extract_count} variants"
    
    # Build plink2 command
    local cmd="plink2"
    
    if [[ "$geno_format" == "plink2" ]]; then
        # Remove extension to get prefix
        local geno_prefix="${pgen%.pgen}"
        cmd+=" --pfile ${geno_prefix}"
    else
        local geno_prefix="${bed%.bed}"
        cmd+=" --bfile ${geno_prefix}"
    fi
    
    cmd+=" --extract ${extract_file}"
    cmd+=" --score ${posteriors_file} ${score_columns} header cols=+scoresums,+denom ignore-dup-ids"
    cmd+=" --out ${chr_workdir}/chr${chr}"
    local plink_threads="${CFG_PLINK_THREADS:-1}"
    if ! [[ "$plink_threads" =~ ^[0-9]+$ ]] || [[ "$plink_threads" -lt 1 ]]; then
        log_warn "Invalid plink.threads='${CFG_PLINK_THREADS:-}', defaulting to 1"
        plink_threads=1
    fi
    cmd+=" --threads ${plink_threads}"
    
    log_debug "Running: $cmd"
    
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log_info "DRY RUN: Would execute: $cmd"
        return 0
    fi
    
    # Execute plink2
    if eval "$cmd" > "${chr_workdir}/plink2.log" 2>&1; then
        local score_file="${chr_workdir}/chr${chr}.sscore"
        
        if [[ -f "$score_file" ]]; then
            # Remove leading # from header if present
            sed -i '1s/^#//' "$score_file"
            
            # Count variants used (from log file)
            local n_variants=0
            if [[ -f "${chr_workdir}/plink2.log" ]]; then
                n_variants=$(grep -oP '\d+(?= variants loaded from --score file)' "${chr_workdir}/plink2.log" || echo "0")
            fi
            
            # Add N_VARIANTS column
            awk -F'\t' -v OFS='\t' -v n_variants="$n_variants" '
                NR == 1 { print $0, "N_VARIANTS"; next }
                { print $0, n_variants }
            ' "$score_file" > "${score_file}.tmp" && mv "${score_file}.tmp" "$score_file"
            
            # Move to final location
            mv "$score_file" "${step_dir}/chr${chr}.sscore"
            
            log_debug "chr${chr}: Scoring completed successfully (${n_variants} variants)"
            return 0
        else
            log_error "chr${chr}: plink2 did not produce score file"
            return 1
        fi
    else
        log_error "chr${chr}: plink2 scoring failed. Check ${chr_workdir}/plink2.log"
        cat "${chr_workdir}/plink2.log" >&2
        return 1
    fi
}

create_empty_score() {
    local output_file="$1"
    
    # Create minimal empty score file with standard header
    echo -e "FID\tIID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_SUM" > "$output_file"
}




