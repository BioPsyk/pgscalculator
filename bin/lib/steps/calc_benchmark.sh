#!/bin/bash
# pgscalculator v2 - calc-benchmark step
# Calculate benchmark scores using MAF filtering and LD pruning

check_calc_benchmark_deps() {
    require_command "plink2" "plink2 is required for benchmark calculation"
    require_command "awk" "awk is required for text processing"
    validate_required_config "CFG" "OUTDIR" "GENODIR" "GENOFILE"
    require_file "${CFG_GENOFILE}" "Genotype manifest file not found"
    require_dir "${CFG_GENODIR}" "Genotype directory not found"
}

run_calc_benchmark() {
    local sumstat_name="$1"
    local maf_threshold="${2:-0.05}"
    
    log_step "Running calc-benchmark for: $sumstat_name"
    check_calc_benchmark_deps
    
    local outdir="${CFG_OUTDIR}"
    local prep_dir=$(get_prep_dir "$outdir")
    local sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")
    local filtered_dir="${sumstat_dir}/filtered"
    local step_dir="${sumstat_dir}/benchmark"
    ensure_dir "$step_dir"
    
    require_dir "$filtered_dir" "Run 'pgscalculator filter-variants' first"
    
    if check_step_completed "$step_dir"; then
        log_info "Step already completed."
        return 0
    fi
    
    local genodir="${CFG_GENODIR}"
    local genofile="${CFG_GENOFILE}"
    local whitelist_file="${prep_dir}/whitelist/variant_whitelist.tsv"
    
    log_info "MAF threshold: ${maf_threshold}"
    
    local success_count=0
    for chr in $(get_chromosomes); do
        local filtered_file="${filtered_dir}/chr${chr}_filtered.tsv"
        [[ ! -f "$filtered_file" ]] && continue
        log_substep "Processing chromosome ${chr}"
        if process_benchmark_chr "$chr" "$filtered_file" "$genodir" "$genofile" "$step_dir" "$maf_threshold" "$whitelist_file"; then
            ((success_count++))
        fi
    done
    
    [[ $success_count -gt 0 ]] && combine_benchmark_scores "$step_dir"
    mark_step_completed "$step_dir"
    log_info "Output directory: ${step_dir}"
}

process_benchmark_chr() {
    local chr="$1" filtered_file="$2" genodir="$3" genofile="$4"
    local step_dir="$5" maf_threshold="$6" whitelist_file="$7"
    
    local chr_workdir="${step_dir}/work_chr${chr}"
    mkdir -p "$chr_workdir"
    
    # Prepare benchmark sumstat with genotype IDs
    local bench_sumstat="${chr_workdir}/bench_sumstat.tsv"
    awk -F'\t' 'NR > 1 {print $1, $2}' "$whitelist_file" > "${chr_workdir}/rsid_map.txt"
    
    awk -F'\t' -v OFS='\t' '
        ARGIND == 1 { rsid_to_geno[$1] = $2; next }
        NR == 1 { for(i=1;i<=NF;i++){if($i=="SNP"||$i=="RSID")sc=i;if($i=="A1")ac=i;if($i=="B"||$i=="BETA")bc=i} print "ID","A1","BETA"; next }
        { if($sc in rsid_to_geno && $bc!="NA" && $bc!="") print rsid_to_geno[$sc],$ac,$bc }
    ' "${chr_workdir}/rsid_map.txt" "$filtered_file" > "$bench_sumstat"
    
    local variant_count=$(tail -n +2 "$bench_sumstat" | wc -l)
    [[ $variant_count -lt 10 ]] && return 0
    
    # Get genotype prefix
    local pgen=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "pgen")
    local geno_opt geno_prefix
    if [[ -n "$pgen" && -f "$pgen" ]]; then
        geno_prefix="${pgen%.pgen}"
        geno_opt="--pfile"
    else
        local bed=$(get_geno_files_for_chr "$genofile" "$genodir" "$chr" "bed")
        [[ -z "$bed" || ! -f "$bed" ]] && return 1
        geno_prefix="${bed%.bed}"
        geno_opt="--bfile"
    fi
    
    # Extract variants
    awk -F'\t' 'NR > 1 {print $1}' "$bench_sumstat" > "${chr_workdir}/variants.txt"
    
    # LD pruning
    plink2 $geno_opt "$geno_prefix" --extract "${chr_workdir}/variants.txt" \
        --maf "$maf_threshold" --indep-pairwise 250 50 0.25 \
        --out "${chr_workdir}/pruned" --threads 1 > "${chr_workdir}/prune.log" 2>&1 || return 1
    
    [[ ! -f "${chr_workdir}/pruned.prune.in" ]] && return 0
    
    # Score with pruned variants
    awk -F'\t' -v OFS='\t' 'ARGIND==1{k[$1]=1;next} NR==1||$1 in k' \
        "${chr_workdir}/pruned.prune.in" "$bench_sumstat" > "${chr_workdir}/score_input.tsv"
    
    plink2 $geno_opt "$geno_prefix" --extract "${chr_workdir}/pruned.prune.in" \
        --score "${chr_workdir}/score_input.tsv" 1 2 3 header cols=scoresums ignore-dup-ids \
        --out "${chr_workdir}/bench" --threads 1 > "${chr_workdir}/score.log" 2>&1 || return 1
    
    [[ -f "${chr_workdir}/bench.sscore" ]] && {
        sed -i '1s/^#//' "${chr_workdir}/bench.sscore"
        mv "${chr_workdir}/bench.sscore" "${step_dir}/chr${chr}_bench.sscore"
        return 0
    }
    return 1
}

combine_benchmark_scores() {
    local step_dir="$1"
    local files=()
    for chr in $(get_chromosomes); do
        [[ -f "${step_dir}/chr${chr}_bench.sscore" ]] && files+=("${step_dir}/chr${chr}_bench.sscore")
    done
    [[ ${#files[@]} -eq 0 ]] && return 0
    
    # Build reference
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="IID"||$i=="#IID")c=i;next}!s[$c]++{print $c}' "${files[0]}" > "${step_dir}/ref.txt"
    
    # Combine
    awk -F'\t' '
        ARGIND==1{r[FNR]=$0;n=FNR;next}
        FNR==1{for(i=1;i<=NF;i++){if($i=="IID"||$i=="#IID")ic=i;if($i=="SCORE1_SUM")sc=i;if($i=="ALLELE_CT")ac=i}next}
        {s[$ic]+=$sc;a[$ic]+=$ac}
        END{print "IID\tALLELE_CT\tSCORE1_SUM";for(i=1;i<=n;i++){id=r[i];print id"\t"a[id]+0"\t"s[id]+0}}
    ' "${step_dir}/ref.txt" "${files[@]}" > "${step_dir}/benchmark.sscore"
}
