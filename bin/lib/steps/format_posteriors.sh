#!/bin/bash
# pgscalculator v2 - format-posteriors step
# Map posteriors to genotype variant IDs for scoring (method-parameterised)

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_format_posteriors_deps() {
    local method="${1:-sbayesr}"

    require_command "awk" "awk is required for text processing"
    require_command "sort" "sort is required for sorting"
    require_command "join" "join is required for file merging"

    validate_required_config "CFG" "OUTDIR"

    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")

    local mapfile="${prep_dir}/variant_map_${method}.tsv"
    local hint
    case "$method" in
        sbayesr)
            hint="Run 'pgscalculator prep-inclusion-list' first"
            ;;
        ldpred2)
            hint="Run 'pgscalculator prep-inclusion-list-ldpred2' first (LDpred2 is opt-in via ldpred2.ld_dir; see plan §5.3)"
            ;;
        *)
            log_error "Unknown format-posteriors method: '${method}' (expected: sbayesr|ldpred2)"
            exit 1
            ;;
    esac
    require_file "$mapfile" "$hint"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_format_posteriors() {
    local sumstat_name="$1"
    local specific_chr="${2:-}"
    # Method dispatch (Phase 1 §5.6). Default sbayesr for back-compat.
    local method="${3:-${CFG_METHOD:-sbayesr}}"

    local posteriors_step mapped_step calc_hint
    case "$method" in
        sbayesr)
            posteriors_step="posteriors"
            mapped_step="posteriors_mapped"
            calc_hint="calc-posteriors"
            ;;
        ldpred2)
            posteriors_step="posteriors_ldpred2"
            mapped_step="posteriors_mapped_ldpred2"
            calc_hint="calc-ldpred2"
            ;;
        *)
            log_error "Invalid format-posteriors method: '${method}' (expected: sbayesr|ldpred2)"
            exit 1
            ;;
    esac

    log_step "Running format-posteriors for: $sumstat_name (method: ${method})"

    check_format_posteriors_deps "$method"

    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors"
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors_mapped"
    migrate_sumstat_step_dir "$sumstat_dir" "$posteriors_step"
    migrate_sumstat_step_dir "$sumstat_dir" "$mapped_step"

    local posteriors_dir
    posteriors_dir=$(get_sumstat_step_dir "$sumstat_dir" "$posteriors_step")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "$mapped_step")
    ensure_dir "$step_dir"

    require_dir "$posteriors_dir" "Run 'pgscalculator ${calc_hint}' first"

    if [[ -z "$specific_chr" ]] && [[ -n "${CFG_CHROMOSOMES:-}" ]] && [[ "${CFG_CHROMOSOMES}" =~ ^(chr)?[0-9]+$ ]]; then
        specific_chr="${CFG_CHROMOSOMES#chr}"
    fi

    if [[ -z "$specific_chr" ]] && check_step_completed "$step_dir"; then
        if ls "${step_dir}"/chr*.snpRes >/dev/null 2>&1; then
            log_info "Step already completed. Use --force to re-run."
            return 0
        fi
        log_warn "Found ${step_dir}/.completed but no mapped chr*.snpRes outputs; re-running format-posteriors."
    fi

    local mapfile="${prep_dir}/variant_map_${method}.tsv"

    log_substep "Building LDREF to genotype ID mapping"
    local rsid_map="${step_dir}/ldref_to_genoid.tsv"
    ensure_rsid_mapping "$mapfile" "$rsid_map"
    if [[ ! -s "$rsid_map" ]]; then
        log_error "LDREF mapping file is empty: ${rsid_map}"
        log_error "Check mapfile: ${mapfile}"
        exit 1
    fi

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
            local output_file="${step_dir}/chr${chr}.snpRes"
            echo -e "ID\tA1\tA2\tFreq\tEffect\tSE\tPIP" > "$output_file"
            echo "posteriors_missing" > "${step_dir}/FAILED_chr${chr}"
            continue
        fi

        local output_file="${step_dir}/chr${chr}.snpRes"

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

    if [[ -z "$specific_chr" ]]; then
        mark_step_completed "$step_dir"
    else
        date '+%Y-%m-%d %H:%M:%S' > "${step_dir}/.completed_chr${specific_chr}"
        log_debug "Marked chr${specific_chr} as completed: ${step_dir}/.completed_chr${specific_chr}"
    fi

    log_info "Total variants mapped: ${total_mapped}"
    log_info "Output directory: ${step_dir}"
}

run_format_posteriors_sbayesr() {
    run_format_posteriors "$1" "${2:-}" "sbayesr"
}

run_format_posteriors_ldpred2() {
    run_format_posteriors "$1" "${2:-}" "ldpred2"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

ensure_rsid_mapping() {
    local mapfile="$1"
    local output_file="$2"

    if [[ -s "$output_file" ]]; then
        return 0
    fi

    local tmp_out
    tmp_out=$(mktemp "${output_file}.tmp.XXXXXX")

    # variant_map_<method>.tsv: col4=geno_snpid, col7=ldref_snpid (9- or 12-col schema)
    awk -F'\t' -v OFS='\t' '
        NR > 1 {
            if ($7 != "NA" && $4 != "NA") {
                print $7, $4
            }
        }
    ' "$mapfile" | LC_ALL=C sort -k1,1 > "$tmp_out"

    local count
    count=$(wc -l < "$tmp_out")
    log_debug "Created mapping with ${count} variants"

    mv "$tmp_out" "$output_file"
}

map_posteriors_for_chr() {
    local chr="$1"
    local posterior_file="$2"
    local rsid_map="$3"
    local output_file="$4"

    local tmpdir
    tmpdir=$(make_tmpdir "format_posteriors")

    tail -n +2 "$posterior_file" | LC_ALL=C sort -k1,1 > "${tmpdir}/posteriors_sorted.tsv"

    LC_ALL=C join -t' ' -1 1 -2 1 \
        -o 2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 \
        "${tmpdir}/posteriors_sorted.tsv" "$rsid_map" 2>/dev/null > "${tmpdir}/mapped.tsv" || true

    # sbayesR / LDpred2 .snpRes: Name (RSID) in column 2
    awk -v OFS='\t' '
        ARGIND == 1 {
            rsid_to_geno[$1] = $2
            next
        }
        FNR == 1 { next }
        {
            rsid = $2
            if (rsid in rsid_to_geno) {
                geno_id = rsid_to_geno[rsid]
                print geno_id, $5, $6, $7, $8, $9, $10
            }
        }
    ' "$rsid_map" "$posterior_file" > "${tmpdir}/mapped_awk.tsv"

    local tmp_out="${output_file}.tmp.$$"
    echo -e "ID\tA1\tA2\tFreq\tEffect\tSE\tPIP" > "$tmp_out"
    cat "${tmpdir}/mapped_awk.tsv" >> "$tmp_out"
    mv "$tmp_out" "$output_file"

    local count
    count=$(wc -l < "${tmpdir}/mapped_awk.tsv")

    rm -rf "$tmpdir"

    echo "$count"
}
