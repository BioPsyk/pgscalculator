#!/bin/bash
# pgscalculator v2 - format-sumstat step
# Format summary statistics: split by chromosome (GRCh38 native coordinates)
# (N/EAF/B/SE derivation happens AFTER filtering in filter-variants step)
#
# Note: With the dual-position variant_map (pos_b37 + pos_b38), we no longer need
# to paste GRCh37 coordinates. The variant_map provides both positions from the
# pre-augmented LD reference.

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
    
    # Only require GRCh38 file (variant_map has both positions from prep)
    require_file "${input_dir}/cleaned_GRCh38.gz" "Cleaned sumstat (GRCh38) not found in input directory"
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
    
    # Split sumstat by chromosome (GRCh38 native coordinates)
    # No paste with GRCh37 needed - variant_map has both positions from prep
    log_substep "Splitting sumstat by chromosome (GRCh38 native coordinates)"
    split_by_chromosome "$input_grch38" "$step_dir"
    
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
    
    log_info "Formatted sumstat: ${total_count} variants across ${chr_count} chromosomes (GRCh38 coordinates)"
    log_info "Output directory: ${step_dir}"
    log_info "Note: Matching to variant_map uses pos_b38; N/EAF/B/SE derivation happens after filtering"
}

# =============================================================================
# PROCESSING FUNCTIONS
# =============================================================================

split_by_chromosome() {
    local input_grch38="$1"
    local outdir="$2"
    
    # Simplified: just split GRCh38 sumstat by chromosome
    # No paste with GRCh37 needed - variant_map has both positions from prep
    
    require_file "$input_grch38" "Cleaned sumstat (GRCh38) not found"
    
    # Single pass: decompress, filter NA chr/pos, split by chromosome
    zcat "$input_grch38" | \
        awk -F'\t' -v OFS='\t' -v outdir="$outdir" '
        NR == 1 {
            header = $0
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
            
            # Skip rows with empty/NA coordinates
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
            for (c in seen) count++
            if (count == 0) {
                print "ERROR: No valid variants written to any chromosome file" > "/dev/stderr"
                exit 1
            }
        }
        '
    
    # Sanity-check at least one chromosome file was created
    local chr_files_count
    chr_files_count=$(ls "${outdir}"/chr*.tsv 2>/dev/null | wc -l)
    if [[ "$chr_files_count" -eq 0 ]]; then
        log_error "split_by_chromosome produced no chromosome files in: ${outdir}"
        exit 1
    fi
}
