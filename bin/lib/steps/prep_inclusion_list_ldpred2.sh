#!/bin/bash
# pgscalculator v2 - prep-inclusion-list-ldpred2 step
# Build the LDpred2 variant map by joining the HM3/HM3+ LDpred2 LD-reference map
# (prep/ldref_ldpred2/map.tsv, produced by prep-ldref-ldpred2) with the genotype
# variant IDs (prep/genotypes/chr{N}_pvar_fmt, produced by prep-genotypes),
# build-aware on chr:pos + alleles. Carries the LD-score and LD-block-id columns
# through so the downstream LDpred2 R script can build its SFBM.
#
# This script is sourced by the main pgscalculator CLI (bin/pgscalculator) and
# by the step runner (bin/lib/steps/run_pipeline.sh).

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_prep_inclusion_list_ldpred2_deps() {
    require_command "awk" "awk is required for text processing"
    validate_required_config "CFG" "OUTDIR"

    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")

    require_dir "${prep_dir}/genotypes"           "Run 'pgscalculator prep-genotypes' first"
    require_file "${prep_dir}/genotypes/snplist_sorted" "Run 'pgscalculator prep-genotypes' first"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_prep_inclusion_list_ldpred2() {
    log_step "Running prep-inclusion-list-ldpred2"

    check_prep_inclusion_list_ldpred2_deps

    # Opt-in via config presence; the sBayesR-only path is the default. Mirrors
    # prep-ldref-ldpred2 (see bin/lib/steps/prep_ldref_ldpred2.sh).
    if [[ -z "${CFG_LDPRED2_LD_DIR:-}" ]]; then
        log_info "ldpred2.ld_dir is not configured; skipping prep-inclusion-list-ldpred2 (LDpred2 is opt-in)"
        return 0
    fi

    local outdir="${CFG_OUTDIR}"
    local prep_dir
    prep_dir=$(get_prep_dir "$outdir")
    local step_dir
    step_dir=$(get_step_dir "$outdir" "inclusion_list_ldpred2")
    ensure_dir "$step_dir"

    local ldpred2_map="${prep_dir}/ldref_ldpred2/map.tsv"
    require_file "$ldpred2_map" "Run 'pgscalculator prep-ldref-ldpred2' first (missing ${ldpred2_map})"

    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi

    local geno_dir="${prep_dir}/genotypes"
    local out_file="${prep_dir}/variant_map_ldpred2.tsv"
    local genotype_build="${CFG_GENOTYPE_BUILD:-GRCh37}"

    log_substep "Joining LDpred2 LD reference with genotype pvars (build: ${genotype_build})"

    echo -e "chr\tpos_b37\tpos_b38\tgeno_snpid\tgeno_a1\tgeno_a2\tldref_snpid\tldref_a1\tldref_a2\tldref_a2freq\tld\tblock_id" \
        > "$out_file"

    local total_rows=0
    local total_matched=0

    for chr in $(get_chromosomes); do
        local pvar_fmt="${geno_dir}/chr${chr}_pvar_fmt"
        if [[ ! -f "$pvar_fmt" ]]; then
            log_warn "Missing genotype data for chr${chr}: ${pvar_fmt}; LDpred2 map will skip this chromosome"
            continue
        fi
        log_debug "Building LDpred2 variant map for chr${chr}"

        local chr_tmp="${step_dir}/chr${chr}.tmp"

        # Awk join semantics mirror create_chr_variant_map() in prep_inclusion_list.sh:
        # left-join from the LD reference (one row per LD-ref variant, NA geno cols on miss),
        # tries direct + swap + complement on both alleles (handles strand ambiguity).
        #
        # Allele convention reconciliation:
        #   bigsnpr a1 (ALT, effect)  -> our ldref_a1 (effect allele)
        #   bigsnpr a0 (REF, other)   -> our ldref_a2 (other allele)
        #   bigsnpr af_UKBB (freq of a1) -> 1 - af_UKBB == freq of a0 == our ldref_a2freq
        # See docs/plans/ldpred2-integration.md §5.3 for rationale.
        awk -F'\t' -v OFS='\t' -v build="$genotype_build" -v want_chr="$chr" '
        BEGIN {
            c["A"]="T"; c["T"]="A"; c["G"]="C"; c["C"]="G"
        }
        # First file: ldpred2 map.tsv (whole-genome). Filter to want_chr and load.
        # Map schema: chr(1) pos_b37(2) pos_b38(3) a0(4) a1(5) rsid(6) af_UKBB(7) ld(8) block_id(9)
        FNR==NR {
            if (FNR == 1) next
            if ($1 != want_chr) next
            la1 = toupper($5)
            la2 = toupper($4)
            af1 = $7 + 0
            n_ld++
            ld_pos_b37[n_ld] = $2
            ld_pos_b38[n_ld] = $3
            ld_a1[n_ld]      = la1
            ld_a2[n_ld]      = la2
            ld_a2freq[n_ld]  = 1 - af1
            ld_rsid[n_ld]    = $6
            ld_ld[n_ld]      = $8
            ld_block[n_ld]   = $9
            if (build == "GRCh38") {
                match_key[n_ld] = $3
            } else {
                match_key[n_ld] = $2
            }
            next
        }
        # Second file: chr{N}_pvar_fmt. Schema: chrpos(1) a1(2) a2(3) snpid(4)
        # Both files are chr-scoped, so position alone is a unique join key here.
        FNR != NR {
            split($1, p, ":")
            pos = p[2]
            ga1 = toupper($2); ga2 = toupper($3); gid = $4
            k1 = pos SUBSEP ga1 SUBSEP ga2
            k2 = pos SUBSEP ga2 SUBSEP ga1
            fa1 = c[ga1]; fa2 = c[ga2]
            k3 = pos SUBSEP fa1 SUBSEP fa2
            k4 = pos SUBSEP fa2 SUBSEP fa1
            geno_id[k1] = gid; geno_a1[k1] = ga1; geno_a2[k1] = ga2
            geno_id[k2] = gid; geno_a1[k2] = ga1; geno_a2[k2] = ga2
            geno_id[k3] = gid; geno_a1[k3] = fa1; geno_a2[k3] = fa2
            geno_id[k4] = gid; geno_a1[k4] = fa1; geno_a2[k4] = fa2
            next
        }
        END {
            for (i = 1; i <= n_ld; i++) {
                mk  = match_key[i]
                la1 = ld_a1[i]; la2 = ld_a2[i]
                k1 = mk SUBSEP la1 SUBSEP la2
                k2 = mk SUBSEP la2 SUBSEP la1
                fa1 = c[la1]; fa2 = c[la2]
                k3 = mk SUBSEP fa1 SUBSEP fa2
                k4 = mk SUBSEP fa2 SUBSEP fa1
                gid = "NA"; ga1 = "NA"; ga2 = "NA"
                if      (k1 in geno_id) { gid = geno_id[k1]; ga1 = geno_a1[k1]; ga2 = geno_a2[k1] }
                else if (k2 in geno_id) { gid = geno_id[k2]; ga1 = geno_a1[k2]; ga2 = geno_a2[k2] }
                else if (k3 in geno_id) { gid = geno_id[k3]; ga1 = geno_a1[k3]; ga2 = geno_a2[k3] }
                else if (k4 in geno_id) { gid = geno_id[k4]; ga1 = geno_a1[k4]; ga2 = geno_a2[k4] }
                print want_chr, ld_pos_b37[i], ld_pos_b38[i], gid, ga1, ga2, \
                      ld_rsid[i], ld_a1[i], ld_a2[i], ld_a2freq[i], ld_ld[i], ld_block[i]
            }
        }
        ' "$ldpred2_map" "$pvar_fmt" > "$chr_tmp"

        local chr_n chr_matched
        chr_n=$(wc -l < "$chr_tmp")
        chr_matched=$(awk -F'\t' '$4 != "NA"' "$chr_tmp" | wc -l)
        total_rows=$((total_rows + chr_n))
        total_matched=$((total_matched + chr_matched))
        log_debug "chr${chr}: ${chr_n} LDpred2 variants, ${chr_matched} matched to genotypes"

        cat "$chr_tmp" >> "$out_file"
        rm -f "$chr_tmp"
    done

    local match_pct="0.00"
    if [[ "$total_rows" -gt 0 ]]; then
        match_pct=$(awk "BEGIN {printf \"%.2f\", 100 * ${total_matched} / ${total_rows}}")
    fi

    local details_dir="${outdir}/prep/details"
    ensure_dir "$details_dir"
    local n_map_input
    n_map_input=$(($(wc -l < "$ldpred2_map") - 1))
    {
        echo -e "STEP\tN_BEFORE\tN_AFTER\tDESC"
        echo -e "ldpred2-map-load\t${n_map_input}\t${total_rows}\tLDpred2 LD-ref map rows (filtered to active autosomes)"
        echo -e "ldpred2-geno-join\t${total_rows}\t${total_matched}\tLD ref ↔ genotype join (match_rate=${match_pct}%)"
    } > "${details_dir}/prep_inclusion_list_ldpred2_steps.tsv"

    mark_step_completed "$step_dir"

    log_info "LDpred2 variant map: ${total_rows} rows, ${total_matched} matched to genotypes (${match_pct}%)"
    log_info "Variant map saved to: ${out_file}"
}
