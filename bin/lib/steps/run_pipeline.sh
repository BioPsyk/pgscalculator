#!/bin/bash
# pgscalculator v2 - run command

declare -A STEP_GROUPS
STEP_GROUPS=([prep]="prep-genotypes prep-ldref prep-ldref-ldpred2 prep-inclusion-list prep-inclusion-list-ldpred2" [sumstat]="format-sumstat filter-variants" [weights]="calc-posteriors format-posteriors calc-benchmark" [score]="calc-score" [finalize]="combine-scores finalize-output")
STEP_GROUP_ORDER=("prep" "sumstat" "weights" "score" "finalize")

format_elapsed() {
    local secs="$1"
    if [[ -z "$secs" ]] || ! [[ "$secs" =~ ^[0-9]+$ ]]; then
        echo "NA"
        return
    fi
    local h=$((secs/3600))
    local m=$(((secs%3600)/60))
    local s=$((secs%60))
    printf "%02d:%02d:%02d" "$h" "$m" "$s"
}

run_pipeline() {
    local sumstat_name="$1" steps_arg="$2" run_all="$3" skip_prep="$4"
    log_step "Running pipeline for: $sumstat_name"
    local pipeline_start_ts
    pipeline_start_ts=$(date +%s)
    # Expose current sumstat name to steps (for tmp/log containment and details reporting)
    if [[ -n "$sumstat_name" ]]; then
        export CFG_SUMSTAT_NAME="$sumstat_name"
    else
        unset CFG_SUMSTAT_NAME || true
    fi
    # Ensure legacy output layouts are migrated into work/
    if [[ -n "$sumstat_name" ]]; then
        local sumstat_dir
        sumstat_dir=$(get_sumstat_dir "${CFG_OUTDIR}" "$sumstat_name")
        migrate_sumstat_all_step_dirs "$sumstat_dir"
    fi
    local -a groups_to_run
    if [[ "$run_all" -eq 1 ]]; then groups_to_run=("${STEP_GROUP_ORDER[@]}"); elif [[ -n "$steps_arg" ]]; then IFS="," read -ra groups_to_run <<< "$steps_arg"; else log_error "Must specify --all or --steps"; exit 1; fi
    [[ "$skip_prep" -eq 1 ]] && check_step_completed "$(get_prep_dir "${CFG_OUTDIR}")/inclusion_list" && groups_to_run=("${groups_to_run[@]/prep/}")
    local failed=0
    for group in "${groups_to_run[@]}"; do [[ -z "$group" ]] && continue; run_step_group "$group" "$sumstat_name" || { failed=1; break; }; done
    local pipeline_end_ts
    pipeline_end_ts=$(date +%s)
    local pipeline_elapsed=$((pipeline_end_ts - pipeline_start_ts))
    if [[ $failed -eq 0 ]]; then
        log_info "Pipeline completed (elapsed=$(format_elapsed "$pipeline_elapsed"))"
    else
        log_error "Pipeline failed (elapsed=$(format_elapsed "$pipeline_elapsed"))"
        exit 1
    fi
}

run_step_group() {
    local group="$1" sumstat_name="$2"; log_substep "Running: $group"
    # `--steps` historically refers to step groups (prep/sumstat/posteriors/score),
    # but the wrapper may also pass concrete step names (e.g. combine-scores,finalize-output).
    # Use a default expansion to avoid "unbound variable" under `set -u`.
    local steps="${STEP_GROUPS[$group]-}"
    if [[ -n "$steps" ]]; then
        for step in $steps; do
            run_single_step "$step" "$sumstat_name" || return 1
        done
        return 0
    fi

    # Fallback: treat the token as a single step name.
    run_single_step "$group" "$sumstat_name" || {
        log_error "Unknown step group/step: '$group' (expected one of: ${!STEP_GROUPS[*]} or a concrete step name)"
        return 1
    }
}

run_single_step() {
    local step="$1" sumstat_name="$2"
    local step_start_ts step_end_ts step_elapsed rc
    log_info "Step: $step"
    step_start_ts=$(date +%s)
    rc=0
    case "$step" in
        prep-genotypes) export CFG_STEP_GROUP="prep"; source "${STEPS_DIR}/prep_genotypes.sh"; run_prep_genotypes; rc=$?;;
        prep-ldref) export CFG_STEP_GROUP="prep"; source "${STEPS_DIR}/prep_ldref.sh"; run_prep_ldref; rc=$?;;
        prep-ldref-ldpred2) export CFG_STEP_GROUP="prep"; source "${STEPS_DIR}/prep_ldref_ldpred2.sh"; run_prep_ldref_ldpred2; rc=$?;;
        prep-inclusion-list) export CFG_STEP_GROUP="prep"; source "${STEPS_DIR}/prep_inclusion_list.sh"; run_prep_inclusion_list; rc=$?;;
        prep-inclusion-list-ldpred2) export CFG_STEP_GROUP="prep"; source "${STEPS_DIR}/prep_inclusion_list_ldpred2.sh"; run_prep_inclusion_list_ldpred2; rc=$?;;
        prep-inclusion-combine) export CFG_STEP_GROUP="prep"; source "${STEPS_DIR}/prep_inclusion_list.sh"; run_prep_inclusion_list_combine; rc=$?;;
        format-sumstat) source "${STEPS_DIR}/format_sumstat.sh"; run_format_sumstat "$sumstat_name"; rc=$?;;
        # filter-variants: method dispatch via CFG_METHOD env var (set from
        # --method by bin/pgscalculator). Default falls back to sbayesr inside
        # run_filter_variants. Per-method driver loops (one filter-variants
        # invocation per requested method) belong to the §6.4 wrapper commit
        # and are intentionally not threaded here yet.
        filter-variants) source "${STEPS_DIR}/filter_variants.sh"; run_filter_variants "$sumstat_name" "" "${CFG_METHOD:-sbayesr}"; rc=$?;;
        calc-posteriors) source "${STEPS_DIR}/calc_posteriors.sh"; run_calc_posteriors "$sumstat_name" ""; rc=$?;;
        calc-ldpred2) source "${STEPS_DIR}/calc_ldpred2.sh"; run_calc_ldpred2 "$sumstat_name"; rc=$?;;
        format-posteriors) source "${STEPS_DIR}/format_posteriors.sh"; run_format_posteriors "$sumstat_name"; rc=$?;;
        calc-score) source "${STEPS_DIR}/calc_score.sh"; run_calc_score "$sumstat_name" ""; rc=$?;;
        combine-scores) source "${STEPS_DIR}/combine_scores.sh"; run_combine_scores "$sumstat_name"; rc=$?;;
        finalize-output) source "${STEPS_DIR}/finalize_output.sh"; run_finalize_output "$sumstat_name"; rc=$?;;
        calc-benchmark) source "${STEPS_DIR}/calc_benchmark.sh"; run_calc_benchmark "$sumstat_name"; rc=$?;;
        *) rc=1;;
    esac
    step_end_ts=$(date +%s)
    step_elapsed=$((step_end_ts - step_start_ts))
    if [[ $rc -eq 0 ]]; then
        log_info "Step completed: ${step} (elapsed=$(format_elapsed "$step_elapsed"))"
    else
        log_error "Step failed: ${step} (elapsed=$(format_elapsed "$step_elapsed"))"
    fi
    return $rc
}
