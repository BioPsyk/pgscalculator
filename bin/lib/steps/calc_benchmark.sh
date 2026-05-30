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

# Compute a benchmark (MAF filter + LD prune + P+T scoring) per method. Each
# method's benchmark is restricted to that method's filtered (LD-ref) variant
# set, so the P+T baseline is judged within the same variant universe as the
# method it is compared against (§11). A single calc-benchmark invocation
# covers every active method; the SLURM driver dispatches one method per task
# via CFG_METHOD.
run_calc_benchmark() {
    local sumstat_name="$1"

    log_step "Running calc-benchmark for: $sumstat_name"
    check_calc_benchmark_deps

    local outdir="${CFG_OUTDIR}"
    local sumstat_dir
    sumstat_dir=$(get_sumstat_dir "$outdir" "$sumstat_name")

    local methods
    if [[ -n "${CFG_METHOD:-}" ]]; then
        methods="${CFG_METHOD}"
    else
        methods="${CFG_METHODS:-sbayesr}"
    fi

    local method ran=0
    for method in $methods; do
        case "$method" in
            sbayesr|ldpred2) ;;
            *) log_warn "calc-benchmark: skipping unknown method '${method}'"; continue ;;
        esac
        migrate_sumstat_step_dir "$sumstat_dir" "$(method_filtered_dir_name "$method")"
        local filtered_dir
        filtered_dir=$(get_method_filtered_dir "$sumstat_dir" "$method")
        if [[ ! -d "$filtered_dir" ]]; then
            log_info "calc-benchmark: no filtered variants for method '${method}' (skipping)"
            continue
        fi
        run_calc_benchmark_method "$sumstat_dir" "$method" "$filtered_dir" && ran=1
    done

    [[ $ran -eq 0 ]] && log_warn "calc-benchmark produced no benchmarks (no filtered variants found)"
    return 0
}

run_calc_benchmark_method() {
    local sumstat_dir="$1" method="$2" filtered_dir="$3"

    # Config from benchmark section (CFG_BENCHMARK_MAF_THRESHOLD, CFG_BENCHMARK_INDEP_PAIRWISE)
    local maf_threshold="${CFG_BENCHMARK_MAF_THRESHOLD:-0.05}"
    local indep_raw="${CFG_BENCHMARK_INDEP_PAIRWISE:-250 50 0.25}"
    local indep_pairwise
    indep_pairwise=$(echo "$indep_raw" | sed 's/[][]//g; s/,/ /g' | awk '{$1=$1;print}')
    [[ -z "$indep_pairwise" ]] && indep_pairwise="250 50 0.25"

    local step_dir
    step_dir=$(get_method_benchmark_dir "$sumstat_dir" "$method")
    ensure_dir "$step_dir"

    if check_step_completed "$step_dir"; then
        log_info "calc-benchmark (${method}): already completed."
        return 0
    fi

    local genodir="${CFG_GENODIR}"
    local genofile="${CFG_GENOFILE}"

    log_substep "Benchmark for method '${method}' (MAF ${maf_threshold}; indep-pairwise ${indep_pairwise})"

    local success_count=0
    for chr in $(get_chromosomes); do
        local filtered_file="${filtered_dir}/chr${chr}_filtered.tsv"
        [[ ! -f "$filtered_file" ]] && continue
        log_substep "[${method}] Processing chromosome ${chr}"
        if process_benchmark_chr "$chr" "$filtered_file" "$genodir" "$genofile" "$step_dir" "$maf_threshold" "$indep_pairwise"; then
            ((success_count++))
        fi
    done

    [[ $success_count -gt 0 ]] && combine_benchmark_scores "$step_dir"
    mark_step_completed "$step_dir"
    log_info "Benchmark (${method}) output directory: ${step_dir}"
    return 0
}

process_benchmark_chr() {
    local chr="$1" filtered_file="$2" genodir="$3" genofile="$4"
    local step_dir="$5" maf_threshold="$6" indep_pairwise="${7:-250 50 0.25}"
    
    local chr_workdir="${step_dir}/work_chr${chr}"
    mkdir -p "$chr_workdir"
    
    # Prepare benchmark sumstat with genotype IDs.
    # Use GENO_ID directly from filtered file (resolved via chr:pos+alleles in variant_map).
    local bench_sumstat="${chr_workdir}/bench_sumstat.tsv"
    awk -F'\t' -v OFS='\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "GENO_ID") gc = i
                if ($i == "EffectAllele" || $i == "A1") ac = i
                if ($i == "B" || $i == "BETA") bc = i
            }
            print "ID", "A1", "BETA"
            next
        }
        gc && ac && bc && $gc != "NA" && $gc != "" && $bc != "NA" && $bc != "" {
            print $gc, $ac, $bc
        }
    ' "$filtered_file" > "$bench_sumstat"
    
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
    
    local plink_threads="${CFG_PLINK_THREADS:-1}"
    if ! [[ "$plink_threads" =~ ^[0-9]+$ ]] || [[ "$plink_threads" -lt 1 ]]; then
        plink_threads=1
    fi

    # LD pruning (window step r2 from config: benchmark.indep_pairwise)
    plink2 $geno_opt "$geno_prefix" --extract "${chr_workdir}/variants.txt" \
        --maf "$maf_threshold" --indep-pairwise $indep_pairwise \
        --out "${chr_workdir}/pruned" --threads "${plink_threads}" > "${chr_workdir}/prune.log" 2>&1 || return 1
    
    [[ ! -f "${chr_workdir}/pruned.prune.in" ]] && return 0
    
    # Score with pruned variants
    awk -F'\t' -v OFS='\t' 'ARGIND==1{k[$1]=1;next} NR==1||$1 in k' \
        "${chr_workdir}/pruned.prune.in" "$bench_sumstat" > "${chr_workdir}/score_input.tsv"
    
    plink2 $geno_opt "$geno_prefix" --extract "${chr_workdir}/pruned.prune.in" \
        --score "${chr_workdir}/score_input.tsv" 1 2 3 header cols=nallele,scoresums ignore-dup-ids \
        --out "${chr_workdir}/bench" --threads "${plink_threads}" > "${chr_workdir}/score.log" 2>&1 || return 1
    
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
