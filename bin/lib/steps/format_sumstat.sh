#!/bin/bash
# pgscalculator v2 - format-sumstat step
# Format summary statistics: add GRCh37 build coordinates only
# (N/EAF/B/SE derivation happens AFTER filtering in filter-variants step)

# This script is sourced by the main pgscalculator CLI

# =============================================================================
# DEPENDENCIES CHECK
# =============================================================================

check_format_sumstat_deps() {
    require_command "awk" "awk is required for text processing"
    require_command "zcat" "zcat is required for reading gzipped files"
    require_command "gzip" "gzip is required for compressing output"
    
    # Check config variables
    validate_required_config "CFG" "INPUT" "OUTDIR"
    
    local input_dir="${CFG_INPUT}"
    
    # Check for cleaned sumstat files from cleansumstats
    require_file "${input_dir}/cleaned_GRCh38.gz" "Cleaned sumstat (GRCh38) not found in input directory"
    require_file "${input_dir}/cleaned_GRCh37.gz" "Cleaned sumstat (GRCh37 map) not found in input directory"
}

# =============================================================================
# MAIN STEP FUNCTION
# =============================================================================

run_format_sumstat() {
    local sumstat_name="$1"
    
    log_step "Running format-sumstat for: $sumstat_name"
    
    # Check dependencies
    check_format_sumstat_deps
    
    # Set up directories
    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    migrate_sumstat_step_dir "$sumstat_dir" "formatted"
    local step_dir
    step_dir=$(get_sumstat_step_dir "$sumstat_dir" "formatted")
    ensure_dir "$step_dir"
    
    # Check if already completed
    if check_step_completed "$step_dir"; then
        log_info "Step already completed. Use --force to re-run."
        return 0
    fi
    
    local input_dir="${CFG_INPUT}"
    local input_grch38="${input_dir}/cleaned_GRCh38.gz"
    local input_grch37="${input_dir}/cleaned_GRCh37.gz"
    
    # Single pass: add GRCh37 coordinates, filter NA, and split by chromosome
    log_substep "Adding GRCh37 build coordinates and splitting by chromosome"
    add_build_coordinates_and_split "$input_grch38" "$input_grch37" "$step_dir"
    
    # Mark step as completed
    mark_step_completed "$step_dir"
    
    # Report results
    local total_count=0
    local chr_count=0
    for chr in $(get_chromosomes); do
        local chr_file="${step_dir}/chr${chr}.tsv"
        if [[ -f "$chr_file" ]]; then
            local count
            count=$(wc -l < "$chr_file")
            count=$((count - 1))
            total_count=$((total_count + count))
            chr_count=$((chr_count + 1))
            log_debug "chr${chr}: ${count} variants"
        fi
    done
    
    log_info "Formatted sumstat with GRCh37 coordinates: ${total_count} variants across ${chr_count} chromosomes"
    log_info "Output directory: ${step_dir}"
    log_info "Note: N/EAF/B/SE derivation will happen after filtering"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

add_build_coordinates_and_split() {
    local input_grch38="$1"
    local input_grch37="$2"
    local outdir="$3"
    
    # GRCh37 file has: CHR, POS, RSID (3 columns)
    # GRCh38 file has: CHR, POS, 0, RSID, EffectAllele, ...
    # After paste (b37 first, then b38):
    #   1: CHR_b37, 2: POS_b37, 3: RSID_b37, 4: CHR_b38, 5: POS_b38, 6+: rest
    # We want: CHR_b37, POS_b37, POS_b38, 0, RSID, ... (drop RSID_b37 and CHR_b38)
    # Use cut -f1-2,5- to skip columns 3 and 4
    #
    # This function combines:
    # - Adding GRCh37 coordinates
    # - Filtering NA coordinates (b37 liftover failures)
    # - Splitting output by chromosome
    # All in a single pass for efficiency.

    require_file "$input_grch37" "Cleaned sumstat (GRCh37 map) not found"
    require_file "$input_grch38" "Cleaned sumstat (GRCh38) not found"

    local tmpdir fifo37 fifo38
    tmpdir=$(make_tmpdir "format_sumstat")
    fifo37="${tmpdir}/grch37.fifo"
    fifo38="${tmpdir}/grch38.fifo"
    mkfifo "$fifo37" "$fifo38"

    # Start streaming decompress in background
    zcat "$input_grch37" > "$fifo37" & local pid37=$!
    zcat "$input_grch38" > "$fifo38" & local pid38=$!

    # Single pass: paste, reorder columns, filter NA, split by chromosome
    paste "$fifo37" "$fifo38" | cut -f1-2,5- | \
        awk -F'\t' -v OFS='\t' -v outdir="$outdir" '
        NR == 1 {
            # Store header for per-chromosome files
            header = $0
            # Find CHR column index (should be column 1 after cut)
            for(i=1; i<=NF; i++) {
                if($i == "CHR" || $i == "chr" || $i == "#CHR") chr_col = i
                if($i == "POS" || $i == "pos") pos_col = i
            }
            if (!chr_col) chr_col = 1
            if (!pos_col) pos_col = 2
            next
        }
        {
            chr = $chr_col
            pos = $pos_col
            
            # Skip rows with empty/NA coordinates (b37 liftover failed)
            if (chr == "" || chr == "NA" || pos == "" || pos == "NA") next
            
            # Only process valid chromosomes (1-22)
            if (chr < 1 || chr > 22) next
            
            # Write to per-chromosome file
            outfile = outdir "/chr" chr ".tsv"
            if (!(chr in seen)) {
                print header > outfile
                seen[chr] = 1
            }
            print > outfile
        }
        END {
            # Report count of chromosomes written
            for (c in seen) count++
            if (count == 0) {
                print "ERROR: No valid variants written to any chromosome file" > "/dev/stderr"
                exit 1
            }
        }
        '

    # Ensure both zcat processes succeeded
    wait "$pid37" || { log_error "Failed to decompress GRCh37 file"; rm -rf "$tmpdir"; exit 1; }
    wait "$pid38" || { log_error "Failed to decompress GRCh38 file"; rm -rf "$tmpdir"; exit 1; }

    rm -rf "$tmpdir"

    # Sanity-check at least one chromosome file was created
    local chr_files_count
    chr_files_count=$(ls "${outdir}"/chr*.tsv 2>/dev/null | wc -l)
    if [[ "$chr_files_count" -eq 0 ]]; then
        log_error "add_build_coordinates_and_split produced no chromosome files in: ${outdir}"
        exit 1
    fi
}
