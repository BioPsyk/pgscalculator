#!/bin/bash
# pgscalculator v2 - prep-ldref-ldpred2 step
# Validate an LDpred2 LD reference directory and convert its map.rds to TSV.
#
# This script is sourced by the main pgscalculator CLI (bin/pgscalculator)
# and by the step runner (bin/lib/steps/run_pipeline.sh).

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_prep_ldref_ldpred2_deps() {
    require_command "Rscript" "Rscript is required for the LDpred2 map conversion"
    validate_required_config "CFG" "OUTDIR"
}

# Default map basename inside CFG_LDPRED2_LD_DIR, derived from the configured
# variant set. Users can override by setting ldpred2.ld_meta_file explicitly.
ldpred2_default_map_basename() {
    local set="${CFG_LDPRED2_LD_VARIANT_SET:-hm3_plus}"
    case "$set" in
        hm3)      echo "map_hm3.rds" ;;
        hm3_plus) echo "map_hm3_plus.rds" ;;
        *)
            log_error "ldpred2.ld_variant_set: '${set}' is not recognised (expected hm3 or hm3_plus)"
            return 1
            ;;
    esac
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_prep_ldref_ldpred2() {
    log_step "Running prep-ldref-ldpred2"

    check_prep_ldref_ldpred2_deps

    # Opt-in via config presence: skip cleanly when LDpred2 is not configured.
    # This keeps sBayesR-only configs working with no change.
    if [[ -z "${CFG_LDPRED2_LD_DIR:-}" ]]; then
        log_info "ldpred2.ld_dir is not configured; skipping prep-ldref-ldpred2 (LDpred2 is opt-in)"
        return 0
    fi

    local outdir="${CFG_OUTDIR}"
    local step_dir
    step_dir=$(get_step_dir "$outdir" "ldref_ldpred2")
    ensure_dir "$step_dir"

    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi

    local ld_dir="${CFG_LDPRED2_LD_DIR}"
    require_dir "$ld_dir" "LDpred2 LD reference directory not found"

    local map_file="${CFG_LDPRED2_LD_META_FILE:-}"
    if [[ -z "$map_file" ]]; then
        local map_basename
        map_basename=$(ldpred2_default_map_basename) || return 1
        map_file="${ld_dir}/${map_basename}"
        log_debug "ldpred2.ld_meta_file not set; defaulting to ${map_file}"
    fi
    require_file "$map_file" "LDpred2 LD map file not found"

    local out_tsv="${step_dir}/map.tsv"
    local r_helper="${SCRIPT_DIR}/R/ldpred2_map_to_tsv.R"
    require_file "$r_helper" "Bundled helper bin/R/ldpred2_map_to_tsv.R not found"

    log_substep "Validating LD blocks and converting map -> TSV"
    log_debug "ld_dir   = ${ld_dir}"
    log_debug "map_file = ${map_file}"
    log_debug "out_tsv  = ${out_tsv}"

    if ! Rscript "$r_helper" \
            --ld-dir   "$ld_dir" \
            --map-file "$map_file" \
            --out-tsv  "$out_tsv"; then
        log_error "ldpred2_map_to_tsv.R failed"
        return 1
    fi

    # Provenance row (mirrors the prep_ldref_steps.tsv pattern used elsewhere
    # under prep/details/ — same audit story across methods).
    local details_dir="${outdir}/prep/details"
    ensure_dir "$details_dir"
    local variant_set="${CFG_LDPRED2_LD_VARIANT_SET:-hm3_plus}"
    local n_variants=$(($(wc -l < "$out_tsv") - 1))
    {
        echo -e "FIELD\tVALUE\tDESC"
        echo -e "ld_dir\t${ld_dir}\tLDpred2 LD reference dir"
        echo -e "map_file\t${map_file}\tLDpred2 LD map .rds"
        echo -e "ld_variant_set\t${variant_set}\tConfigured variant set"
        echo -e "variants\t${n_variants}\tRows in map.tsv (excluding header)"
    } > "${details_dir}/prep_ldref_ldpred2.tsv"

    mark_step_completed "$step_dir"
    log_info "LDpred2 LD map written: ${out_tsv} (${n_variants} variants)"
}
