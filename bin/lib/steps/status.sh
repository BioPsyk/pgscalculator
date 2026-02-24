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
    for step in genotypes ldref inclusion_list; do
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
        for dir in "${outdir}"/sumstats/*; do
            if [[ -d "$dir" ]]; then
                found_sumstat=1
                show_one_sumstat "$outdir" "$(basename "$dir")"
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
    local sd="${outdir}/sumstats/${name}"
    local work="${sd}/work"
    local base="$work"
    echo "=== Sumstat: ${name} ==="
    [[ ! -d "$sd" ]] && { echo "  Not started"; echo ""; return; }
    for step in formatted filtered posteriors posteriors_mapped scores scores_combined; do
        local d="${base}/${step}"
        if [[ -d "$d" ]]; then
            check_step_completed "$d" && echo "  [x] ${step}" || echo "  [~] ${step}"
        else
            echo "  [ ] ${step}"
        fi
    done
    local bd="${sd}/work/benchmark"
    if [[ -d "$bd" ]]; then
        check_step_completed "$bd" && echo "  [x] benchmark" || echo "  [~] benchmark"
    else
        echo "  [ ] benchmark"
    fi
    echo ""
}