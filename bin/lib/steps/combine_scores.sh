#!/bin/bash
# pgscalculator v2 - combine-scores step
# Combine per-chromosome scores into per-method scores_<method>.gz files (§9).

check_combine_scores_deps() {
    require_command "awk" "awk is required for text processing"
    validate_required_config "CFG" "OUTDIR"
}

combine_scores_for_method() {
    local sumstat_name="$1"
    local method="$2"

    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local scores_dir
    scores_dir=$(get_method_scores_dir "$sumstat_dir" "$method") || return 1
    local combine_step="scores_combined_${method}"
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "$combine_step")
    ensure_dir "$step_dir"

    migrate_sumstat_step_dir "$sumstat_dir" "$(method_scores_dir_name "$method")"
    migrate_sumstat_step_dir "$sumstat_dir" "$combine_step"

    if check_step_completed "$step_dir"; then
        log_info "combine-scores (${method}) already completed. Use --force to re-run."
        return 0
    fi

    local out_gz
    out_gz=$(method_scores_gz_name "$method") || return 1
    local out_path="${sumstat_dir}/${out_gz}"

    local score_files=()
    local chr
    for chr in $(get_chromosomes); do
        local score_file="${scores_dir}/chr${chr}.sscore"
        if [[ -f "$score_file" ]]; then
            score_files+=("$score_file")
        fi
    done

    if [[ ${#score_files[@]} -eq 0 ]]; then
        log_warn "No score files for ${method} in ${scores_dir}; writing header-only ${out_gz}."
        local empty_merged="${step_dir}/merged.sscore"
        echo -e "IID\tSCORE_SUM\tALLELE_CT\tN_VARIANTS" > "$empty_merged"
        create_final_scores "$empty_merged" "$out_path"
        mark_step_completed "$step_dir"
        return 0
    fi

    log_info "combine-scores (${method}): ${#score_files[@]} chromosome files"

    local total_nvar
    total_nvar=$(compute_total_scored_variants "$scores_dir")

    local ref_file="${step_dir}/iid_ref.txt"
    build_iid_reference "${score_files[0]}" "$ref_file"

    local merged_file="${step_dir}/merged.sscore"
    combine_chromosome_scores "$ref_file" "$scores_dir" "$merged_file" "$total_nvar"
    create_final_scores "$merged_file" "$out_path"

    if [[ "$method" == "sbayesr" ]]; then
        ln -sf "$out_gz" "${sumstat_dir}/scores.gz"
        log_debug "Symlink: scores.gz -> ${out_gz}"
    fi

    mark_step_completed "$step_dir"

    local sample_count
    sample_count=$(zcat "$out_path" | wc -l)
    sample_count=$((sample_count - 1))
    log_info "Wrote ${out_path} (${sample_count} samples)"
}

run_combine_scores() {
    local sumstat_name="$1"
    local method="${2:-}"

    log_step "Running combine-scores for: $sumstat_name"

    check_combine_scores_deps

    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_all_step_dirs "$sumstat_dir"

    if [[ -n "$method" ]]; then
        case "$method" in
            sbayesr|ldpred2) ;;
            *)
                log_error "Invalid combine-scores method: '${method}'"
                return 1
                ;;
        esac
        combine_scores_for_method "$sumstat_name" "$method"
        return $?
    fi

    local discovered ordered m
    discovered=$(discover_score_methods "$sumstat_dir")
    if [[ -z "$discovered" ]]; then
        log_error "No per-chromosome score files found. Run calc-score first."
        return 1
    fi
    ordered=$(order_discovered_methods "$discovered")
    log_info "Combining scores for method(s): ${ordered}"

    for m in $ordered; do
        combine_scores_for_method "$sumstat_name" "$m" || return 1
    done
}

build_iid_reference() {
    local score_file="$1"
    local ref_file="$2"

    local header
    header=$(head -1 "$score_file")

    local iid_col
    iid_col=$(echo "$header" | awk -F'\t' '{
        for (i=1; i<=NF; i++) {
            if ($i == "IID" || $i == "#IID") {
                print i
                exit
            }
        }
    }')

    if [[ -z "$iid_col" ]]; then
        log_error "Could not find IID column in score file"
        exit 1
    fi

    awk -F'\t' -v iid_col="$iid_col" '
        NR > 1 {
            iid = $iid_col
            if (!(iid in seen)) {
                print iid
                seen[iid] = 1
            }
        }
    ' "$score_file" > "$ref_file"

    local count
    count=$(wc -l < "$ref_file")
    log_debug "Reference: ${count} unique IIDs"
}

combine_chromosome_scores() {
    local ref_file="$1"
    local scores_dir="$2"
    local output_file="$3"
    local total_nvar="$4"

    local score_file_list=""
    for chr in $(get_chromosomes); do
        local score_file="${scores_dir}/chr${chr}.sscore"
        if [[ -f "$score_file" ]]; then
            score_file_list="${score_file_list} ${score_file}"
        fi
    done

    awk -F'\t' -v total_nvar="${total_nvar:-0}" '
        FNR == 1 {
            for (i=1; i<=NF; i++) {
                if ($i == "IID" || $i == "#IID") iid_col = i
                if ($i == "SCORE1_SUM") score_col = i
                if ($i == "ALLELE_CT") allele_col = i
                if ($i == "N_VARIANTS") nvar_col = i
            }
            next
        }
        {
            iid = $iid_col
            if (score_col) scores[iid] += $score_col
            if (allele_col) alleles[iid] += $allele_col
            if (nvar_col) {
                v = $nvar_col + 0
                nvars[iid] += v
                if (v > 0) any_nvar = 1
            }
            if (!(iid in order)) {
                order[iid] = ++n
                iids[n] = iid
            }
        }
        END {
            print "IID\tSCORE_SUM\tALLELE_CT\tN_VARIANTS"
            for (i=1; i<=n; i++) {
                iid = iids[i]
                out_nvar = (any_nvar ? (nvars[iid] + 0) : (total_nvar + 0))
                printf "%s\t%.6g\t%d\t%d\n", iid, scores[iid], alleles[iid]+0, out_nvar
            }
        }
    ' $score_file_list > "$output_file"
}

compute_total_scored_variants() {
    local scores_dir="$1"
    local total=0

    for chr in $(get_chromosomes); do
        local vfile="${scores_dir}/work_chr${chr}/variants.txt"
        if [[ -f "$vfile" ]]; then
            local n
            n=$(wc -l < "$vfile")
            total=$((total + n))
        fi
    done

    echo "$total"
}

create_final_scores() {
    local merged_file="$1"
    local output_file="$2"
    gzip -c "$merged_file" > "$output_file"
    log_debug "Created: $output_file"
}
