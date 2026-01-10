#!/bin/bash
# pgscalculator v2 - format-posteriors step
# Map posteriors to genotype variant IDs for scoring

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_format_posteriors_deps() {
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
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_format_posteriors() {
    local sumstat_name="$1"
    local specific_chr="${2:-}"  # Optional: run only a specific chromosome (also auto-detected from CFG_CHROMOSOMES)
    
    log_step "Running format-posteriors for: $sumstat_name"
    
    # Check dependencies
    check_format_posteriors_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors"
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors_mapped"
    local posteriors_dir
    posteriors_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    ensure_dir "$step_dir"
    
    # Check that calc-posteriors has been run
    require_dir "$posteriors_dir" "Run 'pgscalculator calc-posteriors' first"
    
    # Auto-detect single-chromosome runs from config (important for --sbatch-array mode where config is rewritten)
    if [[ -z "$specific_chr" ]] && [[ -n "${CFG_CHROMOSOMES:-}" ]] && [[ "${CFG_CHROMOSOMES}" =~ ^(chr)?[0-9]+$ ]]; then
        specific_chr="${CFG_CHROMOSOMES#chr}"
    fi

    # Check if already completed (only for full runs).
    # In chromosome-parallel / single-chr runs we *must not* short-circuit on a global marker.
    if [[ -z "$specific_chr" ]] && check_step_completed "$step_dir"; then
        # Be defensive: only skip if we actually have some mapped outputs
        if ls "${step_dir}"/chr*.snpRes >/dev/null 2>&1; then
            log_info "Step already completed. Use --force to re-run."
            return 0
        fi
        log_warn "Found ${step_dir}/.completed but no mapped chr*.snpRes outputs; re-running format-posteriors."
    fi
    
    local inclusion_file="${prep_dir}/inclusion_list/variant_inclusion_list.tsv"
    
    # Build RSID to genotype ID mapping
    log_substep "Building RSID to genotype ID mapping"
    local rsid_map="${step_dir}/rsid_to_genoid.tsv"
    ensure_rsid_mapping "$inclusion_file" "$rsid_map"
    if [[ ! -s "$rsid_map" ]]; then
        log_error "RSID mapping file is empty: ${rsid_map}"
        log_error "Check inclusion list: ${inclusion_file}"
        exit 1
    fi
    
    # Process each chromosome
    log_substep "Mapping posteriors to genotype IDs"
    local total_mapped=0

    local chromosomes
    if [[ -n "$specific_chr" ]]; then
        chromosomes="$specific_chr"
    else
        chromosomes=$(get_chromosomes)
    fi

    local fail_count=0
    for chr in $chromosomes; do
        local posterior_file="${posteriors_dir}/chr${chr}.snpRes"
        
        if [[ ! -f "$posterior_file" ]]; then
            log_warn "chr${chr}: missing posteriors file: ${posterior_file} (writing empty mapped file and continuing)"
            ((fail_count++))
            # Placeholder mapped file for downstream scoring
            local output_file="${step_dir}/chr${chr}.snpRes"
            echo -e "ID\tA1\tA2\tFreq\tEffect\tSE\tPIP" > "$output_file"
            echo "posteriors_missing" > "${step_dir}/FAILED_chr${chr}"
            continue
        fi
        
        local output_file="${step_dir}/chr${chr}.snpRes"

        # Skip if already mapped for this chr
        if [[ -f "$output_file" ]] && [[ $(wc -l < "$output_file") -gt 1 ]]; then
            log_debug "chr${chr}: already mapped, skipping"
            continue
        fi

        local mapped_count
        mapped_count=$(map_posteriors_for_chr "$chr" "$posterior_file" "$rsid_map" "$output_file")
        
        total_mapped=$((total_mapped + mapped_count))
        log_debug "chr${chr}: ${mapped_count} variants mapped"
    done

    if [[ $fail_count -gt 0 ]]; then
        log_warn "format-posteriors had issues for ${fail_count} chromosome(s) (placeholders written; see ${step_dir}/FAILED_chr*)"
    fi
    
    # Mark completion:
    # - Full runs: mark global .completed
    # - Single-chr runs (e.g. sbatch arrays): mark per-chr marker only (avoid blocking other tasks)
    if [[ -z "$specific_chr" ]]; then
        mark_step_completed "$step_dir"
    else
        date '+%Y-%m-%d %H:%M:%S' > "${step_dir}/.completed_chr${specific_chr}"
        log_debug "Marked chr${specific_chr} as completed: ${step_dir}/.completed_chr${specific_chr}"
    fi
    
    log_info "Total variants mapped: ${total_mapped}"
    log_info "Output directory: ${step_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

ensure_rsid_mapping() {
    local inclusion_file="$1"
    local output_file="$2"

    # If mapping exists and is non-empty, keep it (important for parallel runs)
    if [[ -s "$output_file" ]]; then
        return 0
    fi

    local tmp_out
    tmp_out=$(mktemp "${output_file}.tmp.XXXXXX")
    
    # Extract RSID -> pvar_snpid mapping from inclusion list
    # inclusion list format: ld_rsid, pvar_snpid, chrpos
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            print $1, $2  # rsid, genotype_id
        }
    ' "$inclusion_file" | LC_ALL=C sort -k1,1 > "$tmp_out"
    
    local count
    count=$(wc -l < "$tmp_out")
    log_debug "Created mapping with ${count} variants"

    # Atomically move into place (avoids truncate/read races in parallel runs)
    mv "$tmp_out" "$output_file"
}

map_posteriors_for_chr() {
    local chr="$1"
    local posterior_file="$2"
    local rsid_map="$3"
    local output_file="$4"
    
    # sbayesR .snpRes format: SNP A1 A2 b se pval Freq N effect pj
    # We need to:
    # 1. Map SNP (RSID) to genotype variant ID
    # 2. Keep allele information for scoring
    
    local tmpdir
    tmpdir=$(make_tmpdir "format_posteriors")
    
    # Get header from posterior file
    local header
    header=$(head -1 "$posterior_file")
    
    # Sort posteriors by SNP column (column 1)
    tail -n +2 "$posterior_file" | LC_ALL=C sort -k1,1 > "${tmpdir}/posteriors_sorted.tsv"
    
    # Join with RSID mapping
    # Output: genotype_id + all posterior columns
    LC_ALL=C join -t' ' -1 1 -2 1 \
        -o 2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 \
        "${tmpdir}/posteriors_sorted.tsv" "$rsid_map" 2>/dev/null > "${tmpdir}/mapped.tsv" || true
    
    # Alternative: use awk for more robust joining
    # Note: sbayesR .snpRes format has space-padded columns:
    # Id, Name(RSID), Chrom, Position, A1, A2, A1Frq, A1Effect, SE, PIP, LastSampleEff
    # We need to use column 2 (Name) as RSID, and handle variable whitespace
    awk -v OFS='\t' '
        ARGIND == 1 {
            # Load RSID mapping (tab-separated: rsid, genotype_id)
            rsid_to_geno[$1] = $2
            next
        }
        FNR == 1 { next }  # Skip header in posteriors
        {
            # sbayesR output has whitespace-separated columns
            # Column 2 is Name (RSID)
            rsid = $2
            if (rsid in rsid_to_geno) {
                geno_id = rsid_to_geno[rsid]
                # Output: genotype_id, A1, A2, A1Frq, A1Effect, SE, PIP
                # Fields: $5=A1, $6=A2, $7=A1Frq, $8=A1Effect, $9=SE, $10=PIP
                print geno_id, $5, $6, $7, $8, $9, $10
            }
        }
    ' "$rsid_map" "$posterior_file" > "${tmpdir}/mapped_awk.tsv"
    
    # Write output with modified header (atomically)
    # Header: ID A1 A2 Freq Effect SE PIP
    local tmp_out="${output_file}.tmp.$$"
    echo -e "ID\tA1\tA2\tFreq\tEffect\tSE\tPIP" > "$tmp_out"
    cat "${tmpdir}/mapped_awk.tsv" >> "$tmp_out"
    mv "$tmp_out" "$output_file"
    
    # Count mapped variants
    local count
    count=$(wc -l < "${tmpdir}/mapped_awk.tsv")
    
    # Clean up
    rm -rf "$tmpdir"
    
    echo "$count"
}




