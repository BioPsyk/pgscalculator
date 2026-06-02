#!/bin/bash
# pgscalculator v2 - finalize-output step
# Generate final output files: scores.gz, augmented_sumstat.gz, variant_map.gz, bench_score.gz, and details/

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_finalize_output_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "gzip" "gzip is required for compression"
    
    # Check config variables
    validate_required_config "CFG" "OUTDIR"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_finalize_output() {
    local sumstat_name="$1"
    
    log_step "Running finalize-output for: $sumstat_name"
    
    # Check dependencies
    check_finalize_output_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local step_dir="${sumstat_dir}/details"
    ensure_dir "$step_dir"
    
    if ! sumstat_has_any_scores_gz "$sumstat_dir"; then
        log_error "No scores_*.gz found. Run 'pgscalculator combine-scores' first."
        exit 1
    fi
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    # Steps 1-2: Per-method, self-contained outputs (§11).
    # Each method that produced mapped posteriors gets its own augmented_<method>.gz
    # (restricted to that method's LD-reference variant set, with its own
    # benchEffect column kept inside the file) and a per-sample
    # bench_score_<method>.gz. sBayesR additionally gets back-compat aliases
    # (augmented_sumstat.gz, bench_score.gz).
    migrate_sumstat_all_step_dirs "$sumstat_dir"
    local discovered ordered
    discovered=$(discover_posterior_methods "$sumstat_dir")
    ordered=$(order_discovered_methods "$discovered")
    if [[ -z "$ordered" ]]; then
        log_warn "No mapped posteriors on disk; skipping augmented/benchmark outputs"
    else
        log_info "finalize: producing per-method outputs for: ${ordered}"
        local m
        for m in $ordered; do
            log_substep "Per-sample benchmark score (${m})"
            create_bench_score "$sumstat_dir" "$m"
            log_substep "Augmented sumstat (${m})"
            write_augmented_sumstat_for_method "$sumstat_dir" "$prep_dir" "$m"
        done
        link_legacy_finalize_aliases "$sumstat_dir" "$ordered"
    fi
    
    # Step 3: Copy variant map to sumstat root (with rsid as col1)
    log_substep "Writing variant_map.gz"
    write_variant_map "$prep_dir" "$sumstat_dir"

    # Step 4: Copy config to details/
    log_substep "Copying configuration to details/"
    copy_config_to_details "$outdir" "$step_dir"
    
    # Step 5: Generate run summary
    log_substep "Generating run summary"
    generate_run_summary "$sumstat_dir" "$step_dir"

    # Step 6: Generate stepwise details TSVs
    log_substep "Generating stepwise details"
    generate_stepwise_details "$sumstat_dir" "$step_dir"

    # Step 7: Surface LDpred2 run diagnostics under details/ldpred2/
    if has_method ldpred2 "${ordered:-}"; then
        log_substep "Collecting LDpred2 diagnostics"
        generate_ldpred2_details "$sumstat_dir" "$step_dir"
    fi

    # Mark step as completed
    mark_step_completed "$step_dir"
    
    log_info "Finalization complete"
    log_info "Output directory: ${sumstat_dir}"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

write_variant_map() {
    local prep_dir="$1"
    local sumstat_dir="$2"

    local variant_map_src="${sumstat_dir}/variant_map.tsv"
    if [[ ! -f "$variant_map_src" ]]; then
        variant_map_src="${prep_dir}/variant_map_sbayesr.tsv"
    fi
    local variant_map_out="${sumstat_dir}/variant_map.gz"

    if [[ ! -f "$variant_map_src" ]]; then
        log_warn "variant map not found at: ${variant_map_src} (skipping)"
        return 0
    fi

    # Output variant_map with rsid (ldref_snpid) as column 1 for user-facing joins
    # Sumstat map: chr, pos_b37, pos_b38, sumstat_snpid, sumstat_effect, sumstat_other, geno_snpid, geno_a1, geno_a2, ldref_snpid(10), ldref_a1, ldref_a2, ldref_a2freq
    # Prep map: chr, pos_b37, pos_b38, geno_snpid, geno_a1, geno_a2, ldref_snpid(7), ldref_a1, ldref_a2, ldref_a2freq
    awk -F'\t' -v OFS='\t' '
        NR==1 { print "rsid", $0; next }
        NF >= 13 { print $10, $0; next }
        { print $7, $0 }
    ' "$variant_map_src" | gzip -c > "$variant_map_out"
    log_debug "Wrote variant map: ${variant_map_out}"
}

create_bench_score() {
    local sumstat_dir="$1"
    local method="${2:-sbayesr}"
    local bench_dir
    bench_dir=$(get_method_benchmark_dir "$sumstat_dir" "$method")
    local bench_combined="${bench_dir}/benchmark.sscore"
    local out_name
    out_name=$(method_bench_score_gz_name "$method") || return 0
    local output_file="${sumstat_dir}/${out_name}"

    if [[ ! -f "$bench_combined" ]]; then
        log_debug "[${method}] No benchmark.sscore found, skipping ${out_name}"
        return 0
    fi

    gzip -c "$bench_combined" > "$output_file"
    local sample_count
    sample_count=$(wc -l < "$bench_combined")
    sample_count=$((sample_count - 1))
    log_info "Created ${out_name} with ${sample_count} samples"
}

# Back-compat aliases: legacy consumers expect augmented_sumstat.gz / bench_score.gz.
# Point them at the sBayesR method outputs (mirrors scores.gz -> scores_sbayesr.gz).
link_legacy_finalize_aliases() {
    local sumstat_dir="$1" ordered="$2"
    has_method sbayesr "$ordered" || return 0
    local aug bench
    aug=$(method_augmented_gz_name sbayesr)
    bench=$(method_bench_score_gz_name sbayesr)
    if [[ -f "${sumstat_dir}/${aug}" ]]; then
        ln -sf "$aug" "${sumstat_dir}/augmented_sumstat.gz"
        log_debug "Symlink: augmented_sumstat.gz -> ${aug}"
    fi
    if [[ -f "${sumstat_dir}/${bench}" ]]; then
        ln -sf "$bench" "${sumstat_dir}/bench_score.gz"
        log_debug "Symlink: bench_score.gz -> ${bench}"
    fi
}

# augmented_<method>.gz: a fully self-contained, per-method augmented sumstat (§11).
# Restricted to the method's own LD-reference variant set (its variant_map), with
# B/SE/Z/P from the method's matched sumstat, postEffect (+ postp for LDpred2) from
# its mapped posteriors, EAF/MAF, and benchEffect (the method's P+T benchmark weight)
# kept inside the same file. All rows are keyed by rsid (ldref_snpid) for cross-file
# joins via the variant map.
write_augmented_sumstat_for_method() {
    local sumstat_dir="$1" prep_dir="$2" method="$3"

    local variant_map
    variant_map=$(get_method_variant_map "$prep_dir" "$sumstat_dir" "$method") || return 0
    local out_name
    out_name=$(method_augmented_gz_name "$method") || return 0
    local output_file="${sumstat_dir}/${out_name}"

    if [[ ! -f "$variant_map" ]]; then
        log_warn "[${method}] variant map not found (${variant_map}); skipping ${out_name}"
        return 0
    fi

    local mapped_dir
    mapped_dir=$(get_method_posteriors_mapped_dir "$sumstat_dir" "$method")
    if ! method_dir_has_snpres_data "$mapped_dir"; then
        log_warn "[${method}] no mapped posteriors; skipping ${out_name}"
        return 0
    fi

    local filtered_dir bench_dir
    filtered_dir=$(get_method_filtered_dir "$sumstat_dir" "$method")
    bench_dir=$(get_method_benchmark_dir "$sumstat_dir" "$method")
    local eaf_file="${prep_dir}/references/ldref_eaf.tsv"
    local maf_file="${prep_dir}/references/maf_computed.tsv"

    local tmpdir
    tmpdir=$(make_tmpdir "finalize_aug_${method}")

    # Base (keyed by ldref_snpid/rsid): rsid, geno_snpid, EffectAllele, OtherAllele.
    # Prep map:        ldref_snpid=$7,  geno_snpid=$4, ldref_a1=$8,  ldref_a2=$9
    # Sumstat map (>=13): ldref_snpid=$10, geno_snpid=$7, ldref_a1=$11, ldref_a2=$12
    awk -F'\t' -v OFS='\t' '
        NR==1 || $1=="chr" {next}
        NF >= 13 { if ($10 != "NA" && $10 != "") print $10, $7, $11, $12; next }
        $7 != "NA" && $7 != "" { print $7, $4, $8, $9 }
    ' "$variant_map" | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/vm_base.tsv"

    # geno_snpid -> rsid crosswalk (to re-key genotype-keyed inputs to rsid).
    awk -F'\t' -v OFS='\t' '
        NR==1 || $1=="chr" {next}
        NF >= 13 { if ($7 != "NA" && $7 != "" && $10 != "NA" && $10 != "") print $7, $10; next }
        $4 != "NA" && $4 != "" && $4 != "." && $7 != "NA" && $7 != "" { print $4, $7 }
    ' "$variant_map" | LC_ALL=C sort -t $'\t' -k1,1 -u > "${tmpdir}/geno_to_ldref.tsv"

    # sumstat B/SE/Z/P keyed by ldref_snpid (rsid), from chr*_matched.tsv.
    local matched_header=""
    for chr in $(get_chromosomes); do
        local mf="${filtered_dir}/chr${chr}_matched.tsv"
        if [[ -f "$mf" ]]; then
            head -1 "$mf" > "${tmpdir}/matched_header.txt"
            matched_header="${tmpdir}/matched_header.txt"
            break
        fi
    done
    if [[ -z "$matched_header" ]]; then
        log_warn "[${method}] no per-chr matched files (chr*_matched.tsv) in ${filtered_dir}; skipping ${out_name}"
        rm -rf "$tmpdir"
        return 0
    fi

    for chr in $(get_chromosomes); do
        local mf="${filtered_dir}/chr${chr}_matched.tsv"
        [[ -f "$mf" ]] && tail -n +2 "$mf"
    done | awk -F'\t' -v OFS='\t' -v hdr="$matched_header" '
        BEGIN {
            getline line < hdr
            n = split(line, cols, "\t")
            for (i = 1; i <= n; i++) {
                if (cols[i] == "LDREF_SNPID") lc = i
                if (cols[i] == "B" || cols[i] == "BETA") bc = i
                if (cols[i] == "SE") sc = i
                if (cols[i] == "Z") zc = i
                if (cols[i] == "P" || cols[i] == "PVAL") pc = i
            }
        }
        { print $lc, (bc ? $bc : "NA"), (sc ? $sc : "NA"), (zc ? $zc : "NA"), (pc ? $pc : "NA") }
    ' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/ss_sorted.tsv"

    # j1: rsid(1) geno(2) EA(3) OA(4) B(5) SE(6) Z(7) P(8)
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/vm_base.tsv" "${tmpdir}/ss_sorted.tsv" > "${tmpdir}/j1.tsv"

    # Posteriors: .snpRes is keyed by geno_snpid (col1): ID A1 A2 Freq Effect SE PIP.
    # Re-key to rsid via geno_to_ldref -> EAF(Freq), postEffect(Effect) [, postp(PIP)].
    local is_ldpred2=0
    [[ "$method" == "ldpred2" ]] && is_ldpred2=1
    {
        for chr in $(get_chromosomes); do
            local pf="${mapped_dir}/chr${chr}.snpRes"
            [[ -f "$pf" ]] && tail -n +2 "$pf"
        done
    } | awk -F'\t' -v OFS='\t' -v ld="$is_ldpred2" '
        { if (ld) print $1, $4, $5, $7; else print $1, $4, $5 }
    ' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/pp_geno.tsv"

    LC_ALL=C join -t $'\t' "${tmpdir}/geno_to_ldref.tsv" "${tmpdir}/pp_geno.tsv" | \
        awk -F'\t' -v OFS='\t' -v ld="$is_ldpred2" '
            { if (ld) print $2, $3, $4, $5; else print $2, $3, $4 }
        ' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/pp_rsid.tsv"

    # j2 = j1 + EAF + postEffect [+ postp]
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j1.tsv" "${tmpdir}/pp_rsid.tsv" > "${tmpdir}/j2.tsv"

    # Column bookkeeping: after j2 we have 8 base cols + EAF + postEffect [+ postp].
    local eaf_col=9 pe_col=10 postp_col=0 ncol=10
    if [[ $is_ldpred2 -eq 1 ]]; then
        postp_col=11
        ncol=11
    fi

    # benchEffect: benchmark score_input is keyed by geno_snpid (ID A1 BETA);
    # re-key to rsid and left-join.
    local current="${tmpdir}/j2.tsv"
    if [[ -d "$bench_dir" ]]; then
        {
            for chr in $(get_chromosomes); do
                local si="${bench_dir}/work_chr${chr}/score_input.tsv"
                [[ -f "$si" ]] && tail -n +2 "$si"
            done
        } | awk -F'\t' -v OFS='\t' '{print $1, $3}' | \
            LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/bench_geno.tsv"

        LC_ALL=C join -t $'\t' -o 1.2,2.2 "${tmpdir}/geno_to_ldref.tsv" "${tmpdir}/bench_geno.tsv" | \
            LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/bench_sorted.tsv"

        if [[ -s "${tmpdir}/bench_sorted.tsv" ]]; then
            LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "$current" "${tmpdir}/bench_sorted.tsv" > "${tmpdir}/j_bench.tsv"
        else
            awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "$current" > "${tmpdir}/j_bench.tsv"
        fi
    else
        awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "$current" > "${tmpdir}/j_bench.tsv"
    fi
    current="${tmpdir}/j_bench.tsv"
    local bench_col=$((ncol + 1))
    ncol=$bench_col

    # MAF: prefer computed MAF (keyed by geno_snpid), else derive from LD-ref EAF (keyed by rsid).
    if [[ -f "$maf_file" ]] && [[ $(wc -l < "$maf_file") -gt 1 ]]; then
        tail -n +2 "$maf_file" | awk -F'\t' -v OFS='\t' '{print $1, $2}' | \
            LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/geno_maf.tsv"
        LC_ALL=C join -t $'\t' -o 1.2,2.2 "${tmpdir}/geno_to_ldref.tsv" "${tmpdir}/geno_maf.tsv" | \
            LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/maf_sorted.tsv"
        LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "$current" "${tmpdir}/maf_sorted.tsv" > "${tmpdir}/j_final.tsv"
    elif [[ -f "$eaf_file" ]]; then
        tail -n +2 "$eaf_file" | awk -F'\t' -v OFS='\t' '{
            f = $4 + 0; maf = (f > 0.5 ? 1 - f : f); print $1, maf
        }' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/maf_sorted.tsv"
        LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "$current" "${tmpdir}/maf_sorted.tsv" > "${tmpdir}/j_final.tsv"
    else
        awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "$current" > "${tmpdir}/j_final.tsv"
    fi
    local maf_col=$((ncol + 1))

    # Emit: RSID, EffectAllele, OtherAllele, B, SE, Z, P, EAF, MAF, postEffect [, postp_ldpred2], benchEffect
    local header="RSID\tEffectAllele\tOtherAllele\tB\tSE\tZ\tP\tEAF\tMAF\tpostEffect"
    [[ $postp_col -gt 0 ]] && header="${header}\tpostp_ldpred2"
    header="${header}\tbenchEffect"

    local awk_script="${tmpdir}/emit_augmented.awk"
    {
        echo 'BEGIN { OFS = "\t" }'
        printf '%s' '{ line = $1 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8'
        printf ' OFS $%s' "$eaf_col"
        printf ' OFS $%s' "$maf_col"
        printf ' OFS $%s' "$pe_col"
        [[ $postp_col -gt 0 ]] && printf ' OFS $%s' "$postp_col"
        printf ' OFS $%s; print line }\n' "$bench_col"
    } > "$awk_script"

    {
        echo -e "$header"
        awk -F'\t' -f "$awk_script" "${tmpdir}/j_final.tsv"
    } | gzip -c > "$output_file"

    rm -rf "$tmpdir"

    local variant_count
    variant_count=$(zcat "$output_file" | wc -l)
    variant_count=$((variant_count - 1))
    log_info "Generated ${out_name} with ${variant_count} variants"
}

# Copy the LDpred2 run diagnostics (summary.tsv, chains.png) produced by
# calc-ldpred2 into details/ldpred2/ for user-facing inspection.
generate_ldpred2_details() {
    local sumstat_dir="$1"
    local details_dir="$2"

    local src_dir
    src_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_ldpred2")
    local out_dir="${details_dir}/ldpred2"

    local copied=0 f
    for f in summary.tsv chains.png; do
        if [[ -f "${src_dir}/${f}" ]]; then
            ensure_dir "$out_dir"
            cp -f "${src_dir}/${f}" "${out_dir}/${f}"
            copied=1
        fi
    done

    if [[ $copied -eq 1 ]]; then
        log_info "Wrote LDpred2 diagnostics to ${out_dir}/"
    else
        log_debug "No LDpred2 diagnostics (summary.tsv/chains.png) found in ${src_dir}"
    fi
}

copy_config_to_details() {
    local outdir="$1"
    local details_dir="$2"
    
    # Copy config.yaml to details/
    local config_file="${outdir}/config.yaml"
    if [[ -f "$config_file" ]]; then
        cp "$config_file" "${details_dir}/config.yaml"
        log_debug "Copied config.yaml to details/"
    fi
}

generate_run_summary() {
    local sumstat_dir="$1"
    local details_dir="$2"
    
    local summary_file="${details_dir}/run_summary.txt"
    
    {
        echo "pgscalculator v2 Run Summary"
        echo "============================"
        echo ""
        echo "Run completed: $(date)"
        echo ""
        echo "Output files:"
        
        local gz
        for gz in scores_sbayesr.gz scores_ldpred2.gz; do
            if [[ -f "${sumstat_dir}/${gz}" ]]; then
                local score_count
                score_count=$(zcat "${sumstat_dir}/${gz}" | wc -l)
                score_count=$((score_count - 1))
                echo "  - ${gz}: ${score_count} samples"
            fi
        done
        for gz in bench_score_sbayesr.gz bench_score_ldpred2.gz; do
            if [[ -f "${sumstat_dir}/${gz}" ]]; then
                local bench_count
                bench_count=$(zcat "${sumstat_dir}/${gz}" | wc -l)
                bench_count=$((bench_count - 1))
                echo "  - ${gz}: ${bench_count} samples"
            fi
        done
        for gz in augmented_sbayesr.gz augmented_ldpred2.gz; do
            if [[ -f "${sumstat_dir}/${gz}" ]]; then
                local aug_count
                aug_count=$(zcat "${sumstat_dir}/${gz}" | wc -l)
                aug_count=$((aug_count - 1))
                echo "  - ${gz}: ${aug_count} variants"
            fi
        done
        if [[ -f "${sumstat_dir}/variant_map.gz" ]]; then
            local vm_count
            vm_count=$(zcat "${sumstat_dir}/variant_map.gz" | wc -l)
            vm_count=$((vm_count - 1))
            echo "  - variant_map.gz: ${vm_count} variants"
        fi
        
        echo ""
        echo "Chromosomes processed:"
        for chr in $(get_chromosomes); do
            migrate_sumstat_step_dir "$sumstat_dir" "scores"
            if [[ -f "$(get_sumstat_step_dir "$sumstat_dir" "scores")/chr${chr}.sscore" ]]; then
                echo "  - chr${chr}: OK"
            fi
        done
        
    } > "$summary_file"
    
    log_debug "Generated run summary"
}

generate_stepwise_details() {
    local sumstat_dir="$1"
    local details_dir="$2"

    local steps_file="${details_dir}/steps.tsv"

    migrate_sumstat_all_step_dirs "$sumstat_dir"
    local formatted_dir filtered_dir post_dir mapped_dir scores_dir
    formatted_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    filtered_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    post_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors")
    mapped_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    scores_dir=$(get_sumstat_step_dir "$sumstat_dir" "scores")

    # Helpers (avoid hard failure if files missing/broken; report 0)
    local n_formatted n_filtered n_post n_mapped
    n_formatted=0
    n_filtered=0
    n_post=0
    n_mapped=0

    if compgen -G "${formatted_dir}/chr*.tsv" >/dev/null 2>&1; then
        n_formatted=$(for f in "${formatted_dir}"/chr*.tsv; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi
    if compgen -G "${filtered_dir}/chr*_filtered.tsv" >/dev/null 2>&1; then
        n_filtered=$(for f in "${filtered_dir}"/chr*_filtered.tsv; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi
    if compgen -G "${post_dir}/chr*.snpRes" >/dev/null 2>&1; then
        # sum across chr files: (lines - 1)
        n_post=$(for f in "${post_dir}"/chr*.snpRes; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi
    if compgen -G "${mapped_dir}/chr*.snpRes" >/dev/null 2>&1; then
        n_mapped=$(for f in "${mapped_dir}"/chr*.snpRes; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi

    # Failure markers
    local post_fail mapped_fail score_fail
    post_fail=0
    mapped_fail=0
    score_fail=0
    if compgen -G "${post_dir}/work_chr*/FAILED" >/dev/null 2>&1; then
        post_fail=$(ls -1 "${post_dir}"/work_chr*/FAILED 2>/dev/null | wc -l | awk '{print $1}')
    fi
    if compgen -G "${mapped_dir}/FAILED_chr*" >/dev/null 2>&1; then
        mapped_fail=$(ls -1 "${mapped_dir}"/FAILED_chr* 2>/dev/null | wc -l | awk '{print $1}')
    fi
    if compgen -G "${scores_dir}/FAILED_chr*" >/dev/null 2>&1; then
        score_fail=$(ls -1 "${scores_dir}"/FAILED_chr* 2>/dev/null | wc -l | awk '{print $1}')
    fi

    # Matched variants (pre-QC, used for augmented_sumstat)
    local n_matched=0
    if compgen -G "${filtered_dir}/chr*_matched.tsv" >/dev/null 2>&1; then
        n_matched=$(for f in "${filtered_dir}"/chr*_matched.tsv; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi

    # Score counts
    local n_scores_files=0 n_samples=0 n_score_variants=0
    if compgen -G "${scores_dir}/chr*.sscore" >/dev/null 2>&1; then
        n_scores_files=$(ls -1 "${scores_dir}"/chr*.sscore 2>/dev/null | wc -l | awk '{print $1}')
        if compgen -G "${scores_dir}/work_chr*/variants.txt" >/dev/null 2>&1; then
            n_score_variants=$(for f in "${scores_dir}"/work_chr*/variants.txt; do wc -l < "$f"; done | awk '{s+=$1} END{print s+0}')
        fi
    fi
    if [[ -f "${sumstat_dir}/scores.gz" ]]; then
        n_samples=$(gzip -cd "${sumstat_dir}/scores.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_samples" -gt 0 ]]; then n_samples=$((n_samples - 1)); else n_samples=0; fi
    fi

    # Benchmark counts
    local n_bench=0
    local bench_dir
    bench_dir=$(get_method_benchmark_dir "$sumstat_dir" "sbayesr")
    if compgen -G "${bench_dir}/work_chr*/score_input.tsv" >/dev/null 2>&1; then
        n_bench=$(for f in "${bench_dir}"/work_chr*/score_input.tsv; do c=$(wc -l < "$f"); echo $((c-1)); done | awk '{s+=$1} END{print s+0}')
    fi

    # Augmented sumstat count (sBayesR primary; alias augmented_sumstat.gz also points here)
    local n_augmented=0
    local aug_primary="${sumstat_dir}/$(method_augmented_gz_name sbayesr)"
    [[ -f "$aug_primary" ]] || aug_primary="${sumstat_dir}/$(method_augmented_gz_name ldpred2)"
    if [[ -f "$aug_primary" ]]; then
        n_augmented=$(gzip -cd "$aug_primary" 2>/dev/null | wc -l || true)
        if [[ "$n_augmented" -gt 0 ]]; then n_augmented=$((n_augmented - 1)); else n_augmented=0; fi
    fi

    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "format-sumstat\t-\t${n_formatted}\tformatted sumstat rows"
        echo -e "filter-variants\t${n_formatted}\t${n_filtered}\tfiltered sumstat rows (matched=${n_matched})"
        echo -e "calc-posteriors\t${n_filtered}\t${n_post}\tposterior rows (failed_chr=${post_fail})"
        echo -e "format-posteriors\t${n_post}\t${n_mapped}\tmapped posterior rows (failed_chr=${mapped_fail})"
        echo -e "calc-benchmark\t${n_filtered}\t${n_bench}\tbenchmark variants (MAF+LD pruned)"
        echo -e "calc-score\t${n_mapped}\t${n_score_variants}\tper-chr scoring (score_files=${n_scores_files}, failed_chr=${score_fail})"
        echo -e "combine-scores\t${n_score_variants}\t${n_samples}\tsamples scored"
        echo -e "finalize-output\t${n_matched}\t${n_augmented}\taugmented sumstat variants"
    } > "$steps_file"

    log_debug "Wrote details: ${steps_file}"
}
