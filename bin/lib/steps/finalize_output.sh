#!/bin/bash
# pgscalculator v2 - finalize-output step
# Generate final output files: sumstat_augmented.tsv.gz and copy config to details/

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
    
    # Step 2: Generate sumstat_augmented.tsv.gz (full) and augmented_sumstat.gz (v2 reduced schema)
    log_substep "Generating augmented sumstat"
    generate_augmented_sumstat "$sumstat_dir" "$prep_dir" "$posteriors_combined"
    write_augmented_sumstat_v2 "$sumstat_dir" "$prep_dir" "$posteriors_combined"
    
    # Step 3: Copy variant map to sumstat root (with rsid as col1)
    log_substep "Writing variant_map.tsv.gz"
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

generate_augmented_sumstat() {
    local sumstat_dir="$1"
    local prep_dir="$2"
    local posteriors_file="$3"
    local output_file="${sumstat_dir}/sumstat_augmented.tsv.gz"

    # Reuse per-chr matched files from filter-variants (Approach B).
    # These files contain all sumstat columns + LDREF_SNPID + GENO_ID, already restricted
    # to variant-map variants.  Concatenating them avoids re-sorting/re-joining the full
    # formatted sumstat (~17M rows); only the matched set (~1M rows) is sorted + joined
    # with posteriors.

    migrate_sumstat_step_dir "$sumstat_dir" "filtered"
    local filtered_dir
    filtered_dir="$(get_sumstat_step_dir "$sumstat_dir" "filtered")"

    local tmpdir
    tmpdir=$(make_tmpdir "finalize_output")

    # Concatenate per-chr matched files (header from first, body from all)
    local first_chr=true
    for chr in $(get_chromosomes); do
        local mf="${filtered_dir}/chr${chr}_matched.tsv"
        [[ -f "$mf" ]] || continue
        if [[ "$first_chr" == true ]]; then
            head -1 "$mf" > "${tmpdir}/header_raw.txt"
            tail -n +2 "$mf"
            first_chr=false
        else
            tail -n +2 "$mf"
        fi
    done > "${tmpdir}/matched_body.tsv"

    if [[ "$first_chr" == true ]]; then
        log_error "No per-chr matched files (chr*_matched.tsv) found in ${filtered_dir}. Run filter-variants first."
        rm -rf "$tmpdir"; return 1
    fi

    # Determine LDREF_SNPID and GENO_ID column positions from header
    local ldref_col geno_col
    ldref_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="LDREF_SNPID"){print i; exit}}' "${tmpdir}/header_raw.txt")
    geno_col=$(awk -F'\t' '{for(i=1;i<=NF;i++) if($i=="GENO_ID"){print i; exit}}' "${tmpdir}/header_raw.txt")

    # Rearrange: put LDREF_SNPID in col 1, strip it and GENO_ID from original positions,
    # then append GENO_ID at the end.  This gives a predictable layout for join.
    awk -F'\t' -v OFS='\t' -v lc="$ldref_col" -v gc="$geno_col" '{
        printf "%s", $lc
        for (i=1; i<=NF; i++) if (i != lc && i != gc) printf "%s%s", OFS, $i
        printf "%s%s\n", OFS, $gc
    }' "${tmpdir}/matched_body.tsv" | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/ldref_sorted.tsv"

    # Posteriors: RSID, EFFECT, PIP; sort by RSID.
    tail -n +2 "$posteriors_file" | awk -F'\t' -v OFS='\t' '{print $1, $6, $8}' | \
        LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/pp_sorted.tsv"

    # Join on LDREF_SNPID/RSID; -a 1 keeps all matched variants.
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/ldref_sorted.tsv" "${tmpdir}/pp_sorted.tsv" > "${tmpdir}/joined.tsv"

    # Emit: original sumstat cols (minus LDREF_SNPID) + GENO_ID + POST_EFFECT + POST_PIP + IN_ANALYSIS
    # joined layout: ldref_snpid(1), <sumstat cols>(2..NF-3), GENO_ID(NF-2), EFFECT(NF-1), PIP(NF)
    {
        # Header: strip LDREF_SNPID and GENO_ID, re-add GENO_ID with the augmented columns
        awk -F'\t' -v OFS='\t' -v lc="$ldref_col" -v gc="$geno_col" '{
            out = ""
            for (i=1; i<=NF; i++) {
                if (i == lc || i == gc) continue
                if (out != "") out = out OFS
                out = out $i
            }
            print out OFS "GENO_ID" OFS "POST_EFFECT" OFS "POST_PIP" OFS "IN_ANALYSIS"
        }' "${tmpdir}/header_raw.txt"

        # Data: skip col 1 (ldref_snpid); cols 2..NF-3 are sumstat; NF-2=GENO_ID, NF-1=EFFECT, NF=PIP
        awk -F'\t' -v OFS='\t' '{
            geno = $(NF-2); pe = $(NF-1); pp = $NF
            in_analysis = (pe != "NA" ? "Y" : "N")
            n = NF - 3
            for (i = 2; i <= n; i++) printf "%s%s", $i, (i < n ? OFS : "")
            printf "%s%s%s%s%s%s%s%s\n", OFS, geno, OFS, pe, OFS, pp, OFS, in_analysis
        }' "${tmpdir}/joined.tsv"
    } | gzip > "$output_file"

    rm -rf "$tmpdir"

    local variant_count
    variant_count=$(zcat "$output_file" | wc -l)
    variant_count=$((variant_count - 1))
    log_info "Generated augmented sumstat with ${variant_count} variants"
}

# v2 output: augmented_sumstat.gz with same row set as variant map; schema RSID, EffectAllele, OtherAllele, B, SE, Z, P, MAF, postEffect, benchEffect
write_augmented_sumstat_v2() {
    local sumstat_dir="$1"
    local prep_dir="$2"
    local posteriors_file="$3"
    local variant_map="${sumstat_dir}/variant_map.tsv"
    [[ ! -f "$variant_map" ]] && variant_map="${prep_dir}/variant_map.tsv"
    local full_augmented="${sumstat_dir}/sumstat_augmented.tsv.gz"
    local maf_file="${prep_dir}/references/maf_computed.tsv"
    local output_file="${sumstat_dir}/augmented_sumstat.gz"
    
    if [[ ! -f "$variant_map" ]]; then
        log_warn "variant_map.tsv not found, skipping augmented_sumstat.gz"
        return 0
    fi
    if [[ ! -f "$full_augmented" ]]; then
        log_warn "sumstat_augmented.tsv.gz not found, skipping augmented_sumstat.gz"
        return 0
    fi
    
    # Use sort+join pipeline — never loads the large sumstat into memory (avoids OOM).
    # Strategy: build a small base table from variant_map, then chain three joins:
    #   vm_base ⋈ sumstat (B,SE,Z,P)  ⋈ posteriors (postEffect)  ⋈ MAF
    # Each join rearranges so the next join key is in field 1, giving a predictable column layout.
    local tmpdir
    tmpdir=$(make_tmpdir "finalize_augmented_v2")
    
    # Detect variant map column count.
    local vm_ncol
    vm_ncol=$(awk -F'\t' 'NR==1{print NF; exit}' "$variant_map")
    
    # Step 1: Extract base table from variant_map.
    # Output 5 columns: join_key, ldref_snpid, geno_snpid, EffectAllele, OtherAllele
    # 13-col map: join_key=sumstat_snpid($4); ea/oa prefer sumstat alleles with ldref fallback
    # 10-col map: join_key=ldref_snpid($7) — works when sumstat RSIDs are rsIDs matching ldref
    tail -n +2 "$variant_map" | awk -F'\t' -v OFS='\t' -v nc="$vm_ncol" '
        nc >= 13 {
            key=$4; ldref=$10; geno=$7
            ea = ($5!="NA" && $5!="" ? $5 : $11)
            oa = ($6!="NA" && $6!="" ? $6 : $12)
            if (key != "NA" && key != "") print key, ldref, geno, ea, oa
            next
        }
        {
            key=$7; ldref=$7; geno=$4
            ea=$8; oa=$9
            if (key != "NA" && key != "") print key, ldref, geno, ea, oa
        }
    ' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/vm_base.tsv"
    
    # Step 2: Extract sumstat values from sumstat_augmented.tsv.gz (now correctly aligned).
    # Output 5 columns: RSID, B, SE, Z, P (SE may be NA if the sumstat lacks it).
    zcat "$full_augmented" | awk -F'\t' -v OFS='\t' '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i=="RSID"||$i=="rsid"||$i=="SNP"||$i=="ID") rsid_c=i
                if ($i=="B"||$i=="BETA") b_c=i
                if ($i=="SE") se_c=i
                if ($i=="Z") z_c=i
                if ($i=="P"||$i=="PVAL") p_c=i
            }
            next
        }
        {
            rsid = (rsid_c ? $rsid_c : "NA")
            b    = (b_c    ? $b_c    : "NA")
            se   = (se_c   ? $se_c   : "NA")
            z    = (z_c    ? $z_c    : "NA")
            p    = (p_c    ? $p_c    : "NA")
            print rsid, b, se, z, p
        }
    ' | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/ss_sorted.tsv"
    
    # Step 3: Join variant_map base with sumstat on join_key/RSID (-a 1 keeps all vm rows).
    # -o auto ensures unmatched rows get NA-filled columns from file 2.
    # Result (9 cols): key(1), ldref(2), geno(3), ea(4), oa(5), B(6), SE(7), Z(8), P(9)
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/vm_base.tsv" "${tmpdir}/ss_sorted.tsv" > "${tmpdir}/j1.tsv"
    
    # Step 4: Rearrange to put ldref in field 1 and join with posteriors.
    # posteriors: RSID(=ldref_snpid), EFFECT
    tail -n +2 "$posteriors_file" | awk -F'\t' -v OFS='\t' '{print $1, $6}' | \
        LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/pp_sorted.tsv"
    awk -F'\t' -v OFS='\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9}' "${tmpdir}/j1.tsv" | \
        LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/j1_ld.tsv"
    # Result (10 cols): ldref(1), key(2), geno(3), ea(4), oa(5), B(6), SE(7), Z(8), P(9), postEffect(10)
    LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j1_ld.tsv" "${tmpdir}/pp_sorted.tsv" > "${tmpdir}/j2.tsv"
    
    # Step 5: Rearrange to put geno in field 1 and join with MAF.
    # MAF file: GENO_ID(1), MAF(2)
    awk -F'\t' -v OFS='\t' '{print $3, $1, $2, $4, $5, $6, $7, $8, $9, $10}' "${tmpdir}/j2.tsv" | \
        LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/j2_geno.tsv"
    if [[ -f "${maf_file}" ]]; then
        tail -n +2 "${maf_file}" | LC_ALL=C sort -t $'\t' -k1,1 > "${tmpdir}/maf_sorted.tsv"
        LC_ALL=C join -t $'\t' -a 1 -e NA -o auto "${tmpdir}/j2_geno.tsv" "${tmpdir}/maf_sorted.tsv" > "${tmpdir}/j3.tsv"
    else
        awk -F'\t' -v OFS='\t' '{print $0, "NA"}' "${tmpdir}/j2_geno.tsv" > "${tmpdir}/j3.tsv"
    fi
    # j3 layout (11 cols): geno(1), ldref(2), key(3), ea(4), oa(5), B(6), SE(7), Z(8), P(9), postEffect(10), MAF(11)
    
    # Step 6: Emit final output.
    {
        echo -e "RSID\tEffectAllele\tOtherAllele\tB\tSE\tZ\tP\tMAF\tpostEffect\tbenchEffect"
        awk -F'\t' -v OFS='\t' '{
            print $2, $4, $5, $6, $7, $8, $9, $11, $10, "NA"
        }' "${tmpdir}/j3.tsv"
    } | gzip -c > "$output_file"
    rm -rf "$tmpdir"
    log_debug "Wrote v2 augmented sumstat (same row set as variant map): ${output_file}"
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
        if [[ -f "${sumstat_dir}/main_raw_score_all.gz" ]]; then
            echo "  - main_raw_score_all.gz: (v2 main score file)"
        fi
        if [[ -f "${sumstat_dir}/augmented_sumstat.gz" ]]; then
            local aug_count
            aug_count=$(zcat "${sumstat_dir}/augmented_sumstat.gz" | wc -l)
            aug_count=$((aug_count - 1))
            echo "  - augmented_sumstat.gz: ${aug_count} variants"
        fi
        if [[ -f "${sumstat_dir}/sumstat_augmented.tsv.gz" ]]; then
            local var_count
            var_count=$(zcat "${sumstat_dir}/sumstat_augmented.tsv.gz" | wc -l)
            var_count=$((var_count - 1))
            echo "  - sumstat_augmented.tsv.gz: ${var_count} variants"
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
