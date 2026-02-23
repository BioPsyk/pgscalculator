#!/bin/bash
# pgscalculator v2 - finalize-output step
# Generate final output files: augmented_sumstat.gz, variant_map.tsv.gz, and details/

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
    
    # Check prerequisites
    require_file "${sumstat_dir}/scores.tsv.gz" "Run 'pgscalculator combine-scores' first"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    # Step 1: Combine all posteriors into single file
    log_substep "Combining posteriors from all chromosomes"
    local posteriors_combined="${sumstat_dir}/posteriors_combined.tsv"
    migrate_sumstat_step_dir "$sumstat_dir" "posteriors_mapped"
    local posteriors_mapped_dir
    posteriors_mapped_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    combine_posteriors "$posteriors_mapped_dir" "$posteriors_combined"
    
    # Step 2: Create bench_score.gz from benchmark results (if available)
    log_substep "Creating benchmark score output"
    create_bench_score "$sumstat_dir"

    # Step 3: Generate augmented_sumstat.gz (v2 output schema)
    log_substep "Generating augmented sumstat"
    write_augmented_sumstat_v2 "$sumstat_dir" "$prep_dir" "$posteriors_combined"
    
    # Step 4: Copy variant map to sumstat root (with rsid as col1)
    log_substep "Writing variant_map.tsv.gz"
    write_variant_map "$prep_dir" "$sumstat_dir"

    # Step 5: Copy config to details/
    log_substep "Copying configuration to details/"
    copy_config_to_details "$outdir" "$step_dir"
    
    # Step 6: Generate run summary
    log_substep "Generating run summary"
    generate_run_summary "$sumstat_dir" "$step_dir"

    # Step 7: Generate stepwise details TSVs
    log_substep "Generating stepwise details"
    generate_stepwise_details "$sumstat_dir" "$step_dir"
    
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
        variant_map_src="${prep_dir}/variant_map.tsv"
    fi
    local variant_map_out="${sumstat_dir}/variant_map.tsv.gz"

    if [[ ! -f "$variant_map_src" ]]; then
        log_warn "variant_map.tsv not found at: ${variant_map_src} (skipping)"
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
    local bench_dir="${sumstat_dir}/benchmark"
    local bench_combined="${bench_dir}/benchmark.sscore"
    local output_file="${sumstat_dir}/bench_score.gz"

    if [[ ! -f "$bench_combined" ]]; then
        log_debug "No benchmark.sscore found, skipping bench_score.gz"
        return 0
    fi

    gzip -c "$bench_combined" > "$output_file"
    local sample_count
    sample_count=$(wc -l < "$bench_combined")
    sample_count=$((sample_count - 1))
    log_info "Created bench_score.gz with ${sample_count} samples"
}

combine_posteriors() {
    local posteriors_dir="$1"
    local output_file="$2"
    # Use LC_ALL=C sort + join (no in-memory lookup). Posteriors: ID, A1, A2, Freq, Effect, SE, PIP → RSID, A1, A2, FREQ, EFFECT, SE, PIP. ldref_to_genoid: RSID, GENO_ID.
    local tmpdir
    tmpdir=$(make_tmpdir "combine_posteriors")
    local post_body="${tmpdir}/post_body.tsv"
    for chr in $(get_chromosomes); do
        local f="${posteriors_dir}/chr${chr}.snpRes"
        [[ -f "$f" ]] && tail -n +2 "$f" >> "$post_body"
    done
    if [[ ! -s "$post_body" ]]; then
        echo -e "RSID\tGENO_ID\tA1\tA2\tFREQ\tEFFECT\tSE\tPIP" > "$output_file"
        rm -rf "$tmpdir"
        log_info "Combined 0 posterior variants (no chr files)"
        return 0
    fi
    LC_ALL=C sort -t $'\t' -k1,1 "$post_body" > "${tmpdir}/post_sorted.tsv"
    local rsid_to_genoid="${posteriors_dir}/ldref_to_genoid.tsv"
    if [[ -f "$rsid_to_genoid" ]]; then
        tail -n +2 "$rsid_to_genoid" | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/genoid_sorted.tsv"
        {
            echo -e "RSID\tGENO_ID\tA1\tA2\tFREQ\tEFFECT\tSE\tPIP"
            LC_ALL=C join -t $'\t' -a 1 -e NA -o 1.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7 "${tmpdir}/post_sorted.tsv" "${tmpdir}/genoid_sorted.tsv"
        } > "$output_file"
    else
        log_warn "ldref_to_genoid.tsv not found; GENO_ID will be NA in posteriors_combined.tsv"
        {
            echo -e "RSID\tGENO_ID\tA1\tA2\tFREQ\tEFFECT\tSE\tPIP"
            awk -F'\t' -v OFS='\t' '{ print $1, "NA", $2, $3, $4, $5, $6, $7 }' "${tmpdir}/post_sorted.tsv"
        } > "$output_file"
    fi
    local total_variants
    total_variants=$(wc -l < "$output_file")
    total_variants=$((total_variants - 1))
    rm -rf "$tmpdir"
    log_info "Combined ${total_variants} posterior variants"
}

# augmented_sumstat.gz: user-facing output with same row set as variant map.
# Schema: RSID, EffectAllele, OtherAllele, B, SE, Z, P, EAF, MAF, postEffect, benchEffect
# Reads B/SE/Z/P from per-chr matched files, EAF+postEffect from posteriors, MAF from genotypes.
write_augmented_sumstat_v2() {
    local sumstat_dir="$1"
    local prep_dir="$2"
    local posteriors_file="$3"
    local variant_map="${sumstat_dir}/variant_map.tsv"
    [[ ! -f "$variant_map" ]] && variant_map="${prep_dir}/variant_map.tsv"
    local eaf_file="${prep_dir}/references/ldref_eaf.tsv"
    local output_file="${sumstat_dir}/augmented_sumstat.gz"

    if [[ ! -f "$variant_map" ]]; then
        log_warn "variant_map.tsv not found, skipping augmented_sumstat.gz"
        return 0
    fi

    local tmpdir
    tmpdir=$(make_tmpdir "finalize_augmented_v2")

    local bench_dir="${sumstat_dir}/benchmark"
    local maf_file="${prep_dir}/references/maf_computed.tsv"

    migrate_sumstat_step_dir "$sumstat_dir" "filtered"
    local filtered_dir
    filtered_dir="$(get_sumstat_step_dir "$sumstat_dir" "filtered")"

    # Step 1: Extract base table from variant_map (10-col prep format).
    # Output: ldref_snpid(join_key), geno_snpid, EffectAllele, OtherAllele
    # Skip duplicate header rows that may exist in the concatenated file.
    awk -F'\t' -v OFS='\t' '
        NR==1 || $1=="chr" {next}
        $7 != "NA" && $7 != "" { print $7, $4, $8, $9 }
    ' "$variant_map" | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/vm_base.tsv"

    # Step 2: Extract B, SE, Z, P from per-chr matched files (header-aware).
    # Output: LDREF_SNPID, B, SE, Z, P
    # First, grab header from the first available matched file.
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
        log_error "No per-chr matched files (chr*_matched.tsv) found. Run filter-variants first."
        rm -rf "$tmpdir"; return 1
    fi

    # Concatenate all matched file bodies and extract columns.
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
        {
            print $lc, (bc ? $bc : "NA"), (sc ? $sc : "NA"), (zc ? $zc : "NA"), (pc ? $pc : "NA")
        }
    ' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/ss_sorted.tsv"

    # Step 3: Join variant_map with sumstat on ldref_snpid (-a 1 keeps all vm rows).
    # Result (8 cols): ldref(1), geno(2), ea(3), oa(4), B(5), SE(6), Z(7), P(8)
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/vm_base.tsv" "${tmpdir}/ss_sorted.tsv" > "${tmpdir}/j1.tsv"

    # Step 4: Join with posteriors on ldref_snpid (already in col 1).
    # Extract FREQ (EAF used for posterior calc) and EFFECT (posterior effect size).
    tail -n +2 "$posteriors_file" | awk -F'\t' -v OFS='\t' '{print $1, $5, $6}' | \
        LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/pp_sorted.tsv"
    # Result (10 cols): ldref(1), geno(2), ea(3), oa(4), B(5), SE(6), Z(7), P(8), EAF(9), postEffect(10)
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j1.tsv" "${tmpdir}/pp_sorted.tsv" > "${tmpdir}/j2.tsv"

    # Step 5: Join with benchmark effects (benchEffect).
    # Benchmark score inputs use genotype IDs; map to ldref_snpid via variant_map.
    if [[ -d "$bench_dir" ]]; then
        awk -F'\t' -v OFS='\t' '
            NR==1 || $1=="chr" {next}
            $4 != "NA" && $4 != "" && $7 != "NA" && $7 != "" { print $4, $7 }
        ' "$variant_map" | LC_ALL=C sort -t $'\t' -k1,1 -u > "${tmpdir}/geno_to_ldref.tsv"

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
            LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j2.tsv" "${tmpdir}/bench_sorted.tsv" > "${tmpdir}/j3.tsv"
        else
            awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "${tmpdir}/j2.tsv" > "${tmpdir}/j3.tsv"
        fi
    else
        awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "${tmpdir}/j2.tsv" > "${tmpdir}/j3.tsv"
    fi
    # j3: ldref(1), geno(2), ea(3), oa(4), B(5), SE(6), Z(7), P(8), EAF(9), postEffect(10), benchEffect(11)

    # Step 6: Join with MAF. Prefer genotype-based MAF, fall back to ldref EAF.
    if [[ -f "$maf_file" ]] && [[ $(wc -l < "$maf_file") -gt 1 ]]; then
        [[ ! -f "${tmpdir}/geno_to_ldref.tsv" ]] && \
            awk -F'\t' -v OFS='\t' '
                NR==1 || $1=="chr" {next}
                $4 != "NA" && $4 != "" && $7 != "NA" && $7 != "" { print $4, $7 }
            ' "$variant_map" | LC_ALL=C sort -t $'\t' -k1,1 -u > "${tmpdir}/geno_to_ldref.tsv"
        tail -n +2 "$maf_file" | awk -F'\t' -v OFS='\t' '{print $1, $2}' | \
            LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/geno_maf.tsv"
        LC_ALL=C join -t $'\t' -o 1.2,2.2 "${tmpdir}/geno_to_ldref.tsv" "${tmpdir}/geno_maf.tsv" | \
            LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/maf_sorted.tsv"
        LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j3.tsv" "${tmpdir}/maf_sorted.tsv" > "${tmpdir}/j4.tsv"
    elif [[ -f "$eaf_file" ]]; then
        tail -n +2 "$eaf_file" | awk -F'\t' -v OFS='\t' '{
            f = $4 + 0; maf = (f > 0.5 ? 1 - f : f); print $1, maf
        }' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/maf_sorted.tsv"
        LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j3.tsv" "${tmpdir}/maf_sorted.tsv" > "${tmpdir}/j4.tsv"
    else
        awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "${tmpdir}/j3.tsv" > "${tmpdir}/j4.tsv"
    fi
    # j4: ldref(1),geno(2),ea(3),oa(4),B(5),SE(6),Z(7),P(8),EAF(9),postEffect(10),benchEffect(11),MAF(12)

    # Step 7: Emit final output.
    {
        echo -e "RSID\tEffectAllele\tOtherAllele\tB\tSE\tZ\tP\tEAF\tMAF\tpostEffect\tbenchEffect"
        awk -F'\t' -v OFS='\t' '{
            print $1, $3, $4, $5, $6, $7, $8, $9, $12, $10, $11
        }' "${tmpdir}/j4.tsv"
    } | gzip -c > "$output_file"

    rm -rf "$tmpdir"

    local variant_count
    variant_count=$(zcat "$output_file" | wc -l)
    variant_count=$((variant_count - 1))
    log_info "Generated augmented sumstat with ${variant_count} variants"
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
        
        if [[ -f "${sumstat_dir}/scores.tsv.gz" ]]; then
            local score_count
            score_count=$(zcat "${sumstat_dir}/scores.tsv.gz" | wc -l)
            score_count=$((score_count - 1))
            echo "  - scores.tsv.gz: ${score_count} samples"
        fi
        if [[ -f "${sumstat_dir}/bench_score.gz" ]]; then
            local bench_count
            bench_count=$(zcat "${sumstat_dir}/bench_score.gz" | wc -l)
            bench_count=$((bench_count - 1))
            echo "  - bench_score.gz: ${bench_count} samples"
        fi
        if [[ -f "${sumstat_dir}/augmented_sumstat.gz" ]]; then
            local aug_count
            aug_count=$(zcat "${sumstat_dir}/augmented_sumstat.gz" | wc -l)
            aug_count=$((aug_count - 1))
            echo "  - augmented_sumstat.gz: ${aug_count} variants"
        fi
        if [[ -f "${sumstat_dir}/posteriors_combined.tsv" ]]; then
            local post_count
            post_count=$(wc -l < "${sumstat_dir}/posteriors_combined.tsv")
            post_count=$((post_count - 1))
            echo "  - posteriors_combined.tsv: ${post_count} posteriors"
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
    local score_file="${details_dir}/score_steps.tsv"

    migrate_sumstat_all_step_dirs "$sumstat_dir"
    local formatted_dir filtered_dir post_dir mapped_dir scores_dir
    formatted_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    filtered_dir=$(get_sumstat_step_dir "$sumstat_dir" "filtered")
    post_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors")
    mapped_dir=$(get_sumstat_step_dir "$sumstat_dir" "posteriors_mapped")
    scores_dir=$(get_sumstat_step_dir "$sumstat_dir" "scores")

    # Helpers (avoid hard failure if files missing/broken; report 0)
    local n_formatted n_filtered n_post n_mapped n_scores_files n_samples n_score_variants
    n_formatted=0
    n_filtered=0
    n_post=0
    n_mapped=0

    if [[ -f "${formatted_dir}/sumstat_formatted.tsv.gz" ]]; then
        n_formatted=$(gzip -cd "${formatted_dir}/sumstat_formatted.tsv.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_formatted" -gt 0 ]]; then n_formatted=$((n_formatted - 1)); else n_formatted=0; fi
    fi
    if [[ -f "${filtered_dir}/sumstat_filtered.tsv.gz" ]]; then
        n_filtered=$(gzip -cd "${filtered_dir}/sumstat_filtered.tsv.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_filtered" -gt 0 ]]; then n_filtered=$((n_filtered - 1)); else n_filtered=0; fi
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

    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "format-sumstat\t${n_formatted}\t${n_formatted}\tformatted sumstat rows"
        echo -e "filter-variants\t${n_formatted}\t${n_filtered}\tfiltered sumstat rows"
        echo -e "calc-posteriors\t${n_filtered}\t${n_post}\tposterior rows (failed_chr=${post_fail})"
        echo -e "format-posteriors\t${n_post}\t${n_mapped}\tmapped posterior rows (failed_chr=${mapped_fail})"
    } > "$steps_file"

    # Score summary
    n_scores_files=0
    n_samples=0
    n_score_variants=0
    if compgen -G "${scores_dir}/chr*.sscore" >/dev/null 2>&1; then
        n_scores_files=$(ls -1 "${scores_dir}"/chr*.sscore 2>/dev/null | wc -l | awk '{print $1}')
        # total variants scored: prefer work_chr*/variants.txt if present
        if compgen -G "${scores_dir}/work_chr*/variants.txt" >/dev/null 2>&1; then
            n_score_variants=$(for f in "${scores_dir}"/work_chr*/variants.txt; do wc -l < "$f"; done | awk '{s+=$1} END{print s+0}')
        fi
    fi
    if [[ -f "${sumstat_dir}/scores.tsv.gz" ]]; then
        n_samples=$(gzip -cd "${sumstat_dir}/scores.tsv.gz" 2>/dev/null | wc -l || true)
        if [[ "$n_samples" -gt 0 ]]; then n_samples=$((n_samples - 1)); else n_samples=0; fi
    fi

    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "calc-score\t${n_mapped}\t${n_score_variants}\tper-chr scoring ran (score_files=${n_scores_files}, failed_chr=${score_fail})"
        echo -e "combine-scores\t${n_score_variants}\t${n_score_variants}\tcombined score table written"
        echo -e "finalize-output\t${n_score_variants}\t${n_score_variants}\tfinal outputs written (samples=${n_samples})"
    } > "$score_file"

    log_debug "Wrote details: ${steps_file}, ${score_file}"
}
