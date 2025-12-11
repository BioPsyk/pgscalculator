#!/bin/bash
# pgscalculator v2 - run command

declare -A STEP_GROUPS
STEP_GROUPS=([prep]="prep-genotypes prep-ldref prep-whitelist" [sumstat]="format-sumstat filter-variants" [posteriors]="calc-posteriors format-posteriors" [score]="calc-score combine-scores" [benchmark]="calc-benchmark")
STEP_GROUP_ORDER=("prep" "sumstat" "posteriors" "score")

run_pipeline() {
    local sumstat_name="$1" steps_arg="$2" run_all="$3" skip_prep="$4"
    log_step "Running pipeline for: $sumstat_name"
    local -a groups_to_run
    if [[ "$run_all" -eq 1 ]]; then groups_to_run=("${STEP_GROUP_ORDER[@]}"); elif [[ -n "$steps_arg" ]]; then IFS="," read -ra groups_to_run <<< "$steps_arg"; else log_error "Must specify --all or --steps"; exit 1; fi
    [[ "$skip_prep" -eq 1 ]] && check_step_completed "$(get_prep_dir "${CFG_OUTDIR}")/whitelist" && groups_to_run=("${groups_to_run[@]/prep/}")
    local failed=0
    for group in "${groups_to_run[@]}"; do [[ -z "$group" ]] && continue; run_step_group "$group" "$sumstat_name" || { failed=1; break; }; done
    [[ $failed -eq 0 ]] && log_info "Pipeline completed" || { log_error "Pipeline failed"; exit 1; }
}

run_step_group() {
    local group="$1" sumstat_name="$2"; log_substep "Running: $group"
    local steps="${STEP_GROUPS[$group]}"; [[ -z "$steps" ]] && return 1
    for step in $steps; do run_single_step "$step" "$sumstat_name" || return 1; done
}

run_single_step() {
    local step="$1" sumstat_name="$2"; log_info "Step: $step"
    case "$step" in
        prep-genotypes) source "${STEPS_DIR}/prep_genotypes.sh"; run_prep_genotypes;;
        prep-ldref) source "${STEPS_DIR}/prep_ldref.sh"; run_prep_ldref;;
        prep-whitelist) source "${STEPS_DIR}/prep_whitelist.sh"; run_prep_whitelist;;
        format-sumstat) source "${STEPS_DIR}/format_sumstat.sh"; run_format_sumstat "$sumstat_name";;
        filter-variants) source "${STEPS_DIR}/filter_variants.sh"; run_filter_variants "$sumstat_name";;
        calc-posteriors) source "${STEPS_DIR}/calc_posteriors.sh"; run_calc_posteriors "$sumstat_name" "";;
        format-posteriors) source "${STEPS_DIR}/format_posteriors.sh"; run_format_posteriors "$sumstat_name";;
        calc-score) source "${STEPS_DIR}/calc_score.sh"; run_calc_score "$sumstat_name" "";;
        combine-scores) source "${STEPS_DIR}/combine_scores.sh"; run_combine_scores "$sumstat_name";;
        calc-benchmark) source "${STEPS_DIR}/calc_benchmark.sh"; run_calc_benchmark "$sumstat_name" "";;
        *) return 1;;
    esac
}
