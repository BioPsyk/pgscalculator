#!/bin/bash
# pgscalculator v2 - status command

show_status() {
    local sumstat_name="${1:-}"
    local outdir="${CFG_OUTDIR}"
    log_step "Pipeline Status"
    [[ ! -d "$outdir" ]] && { log_info "No runs found."; return 0; }
    echo "Output: $outdir"
    echo ""
    echo "=== Prep Steps ==="
    local prep_dir="${outdir}/prep"
    for step in genotypes ldref whitelist; do
        local d="${prep_dir}/${step}"
        if [[ -d "$d" ]]; then
            check_step_completed "$d" && echo "  [x] prep-${step}" || echo "  [~] prep-${step}"
        else
            echo "  [ ] prep-${step}"
        fi
    done
    echo ""
    if [[ -n "$sumstat_name" ]]; then
        show_one_sumstat "$outdir" "$sumstat_name"
    else
        # Check if any sumstat directories exist
        local found_sumstat=0
        for dir in "${outdir}"/sumstat_*; do
            if [[ -d "$dir" ]]; then
                found_sumstat=1
                show_one_sumstat "$outdir" "$(basename "$dir" | sed 's/^sumstat_//')"
            fi
        done
        if [[ $found_sumstat -eq 0 ]]; then
            echo "=== Sumstats ==="
            echo "  No sumstats processed yet"
            echo ""
        fi
    fi
}

show_one_sumstat() {
    local outdir="$1" name="$2"
    local sd="${outdir}/sumstat_${name}"
    echo "=== Sumstat: ${name} ==="
    [[ ! -d "$sd" ]] && { echo "  Not started"; echo ""; return; }
    for step in formatted filtered posteriors posteriors_mapped scores scores_combined benchmark; do
        local d="${sd}/${step}"
        if [[ -d "$d" ]]; then
            check_step_completed "$d" && echo "  [x] ${step}" || echo "  [~] ${step}"
        else
            echo "  [ ] ${step}"
        fi
    done
    echo ""
}