#!/bin/bash
# pgscalculator v2 - calc-ldpred2 step
# Genome-wide LDpred2 posterior estimation (single job, not chr-parallel).

check_calc_ldpred2_deps() {
    require_command "Rscript" "Rscript is required for LDpred2"
    validate_required_config "CFG" "OUTDIR"

    if [[ -z "${CFG_LDPRED2_LD_DIR:-}" ]]; then
        log_error "ldpred2.ld_dir is not configured."
        log_error "Set ldpred2.ld_dir in config.yaml (see config.template.yaml)."
        exit 1
    fi
}

write_empty_ldpred2_posteriors() {
    local step_dir="$1"
    local hdr="Id Name Chrom Position A1 A2 A1Frq A1Effect SE PIP LastSampleEff"
    local chr
    for chr in $(get_chromosomes); do
        echo "$hdr" > "${step_dir}/chr${chr}.snpRes"
    done
}

run_calc_ldpred2() {
    local sumstat_name="$1"

    log_step "Running calc-ldpred2 for: $sumstat_name"

    check_calc_ldpred2_deps

    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_all_step_dirs "$sumstat_dir"

    local filter_dir
    filter_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered_ldpred2")
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_ldpred2")
    ensure_dir "$step_dir"

    require_dir "$filter_dir" "Run 'pgscalculator filter-variants --method ldpred2' first"

    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi

    local ld_dir="${CFG_LDPRED2_LD_DIR}"
    require_dir "$ld_dir" "LDpred2 LD reference directory not found"

    local ld_meta="${CFG_LDPRED2_LD_META_FILE_RESOLVED:-}"
    if [[ -z "$ld_meta" ]]; then
        ld_meta=$(resolve_ldpred2_ld_meta_file) || return 1
    fi
    require_file "$ld_meta" "LDpred2 LD map file not found"

    local r_script="${SCRIPT_DIR}/lib/scripts/run_ldpred2.R"
    require_file "$r_script" "Bundled helper bin/lib/scripts/run_ldpred2.R not found"

    local mode="${CFG_LDPRED2_MODE:-auto}"
    case "$mode" in
        auto|inf) ;;
        *)
            log_error "ldpred2.mode: '${mode}' is not recognised (expected auto or inf)"
            return 1
            ;;
    esac

    local shrink_corr="${CFG_LDPRED2_SHRINK_CORR:-0.95}"
    local allow_jump_sign="${CFG_LDPRED2_ALLOW_JUMP_SIGN:-false}"
    local hyper_p_max="${CFG_LDPRED2_HYPER_P_MAX:-0.2}"
    local hyper_p_length="${CFG_LDPRED2_HYPER_P_LENGTH:-30}"
    local seed="${CFG_LDPRED2_SEED:-1}"
    local ncores="${CFG_LDPRED2_THREADS:-${CFG_LDPRED2_NCORES:-1}}"
    local genotype_build="${CFG_GENOTYPE_BUILD:-GRCh37}"
    local ld_build="${CFG_LDPRED2_LD_BUILD:-GRCh37}"
    local merge_by_rsid="${CFG_LDPRED2_MERGE_BY_RSID:-false}"

    local logs_dir="${sumstat_dir}/logs"
    ensure_dir "$logs_dir"
    local plot_file="${logs_dir}/ldpred2_chains.png"
    local log_file="${logs_dir}/calc_ldpred2.log"

    local metadata_file="${CFG_INPUT}/cleaned_metadata.yaml"
    local which_n="${CFG_WHICHN:-totalN}"
    local -a meta_args=()
    local eff_n case_n ctrl_n total_n

    if [[ -f "$metadata_file" ]]; then
        if [[ "$which_n" == "effectiveN" ]]; then
            eff_n=$(awk -F': ' '$1=="stats_EffectiveN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
            case_n=$(awk -F': ' '$1=="stats_CaseN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
            ctrl_n=$(awk -F': ' '$1=="stats_ControlN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
            if [[ -n "$eff_n" ]] && [[ "$eff_n" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                meta_args=(--effective-sample-size "$eff_n")
            elif [[ -n "$case_n" ]] && [[ -n "$ctrl_n" ]] \
                && [[ "$case_n" =~ ^[0-9]+$ ]] && [[ "$ctrl_n" =~ ^[0-9]+$ ]]; then
                meta_args=(--n-cases "$case_n" --n-controls "$ctrl_n")
            fi
        else
            total_n=$(awk -F': ' '$1=="stats_TotalN"{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2; exit}' "$metadata_file")
            if [[ -n "$total_n" ]] && [[ "$total_n" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                meta_args=(--effective-sample-size "$total_n")
            fi
        fi
    fi

    local allow_jump_r="FALSE"
    [[ "${allow_jump_sign,,}" == "true" ]] && allow_jump_r="TRUE"
    local merge_r="FALSE"
    [[ "${merge_by_rsid,,}" == "true" ]] && merge_r="TRUE"

    local cmd=(
        Rscript "$r_script"
        --sumstat-dir "$filter_dir"
        --ld-dir "$ld_dir"
        --ld-meta "$ld_meta"
        --out-dir "$step_dir"
        --mode "$mode"
        --shrink-corr "$shrink_corr"
        --allow-jump-sign "$allow_jump_r"
        --hyper-p-max "$hyper_p_max"
        --hyper-p-length "$hyper_p_length"
        --seed "$seed"
        --ncores "$ncores"
        --genotype-build "$genotype_build"
        --ld-build "$ld_build"
        --merge-by-rsid "$merge_r"
        --plot-file "$plot_file"
    )
    cmd+=("${meta_args[@]}")

    log_info "LDpred2 mode: ${mode}"
    log_debug "LD dir: ${ld_dir}"
    log_debug "LD meta: ${ld_meta}"
    log_debug "Filtered dir: ${filter_dir}"
    log_debug "Output dir: ${step_dir}"

    if [[ "${DRY_RUN}" -eq 1 ]]; then
        log_info "DRY RUN: Would execute: ${cmd[*]}"
        return 0
    fi

    local rc=0
    if ! "${cmd[@]}" > "$log_file" 2>&1; then
        rc=1
        log_warn "calc-ldpred2 failed; see ${log_file}"
        if compgen -G "${step_dir}/chr*.snpRes" >/dev/null; then
            :
        else
            write_empty_ldpred2_posteriors "$step_dir"
        fi
        for marker in FAILED_match FAILED_sfbm FAILED_ldsc FAILED_ldpred2; do
            if [[ -f "${step_dir}/${marker}" ]]; then
                log_warn "Failure marker: ${step_dir}/${marker}"
            fi
        done
    else
        log_info "LDpred2 completed successfully"
    fi

    mark_step_completed "$step_dir"

    if [[ $rc -ne 0 ]]; then
        log_warn "calc-ldpred2 finished with errors (empty/header-only posteriors may have been written)"
        return 0
    fi

    log_info "Output directory: ${step_dir}"
    return 0
}
