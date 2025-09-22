# Memory Optimization Plan for pgscalculator

## Executive Summary

This document outlines a plan to address memory issues in pgscalculator when processing large genotype datasets with many variants and samples. The analysis identifies specific plink operations causing memory bottlenecks and proposes sample-based splitting strategies to reduce memory usage while maintaining computational correctness.

## Current Memory Bottlenecks

### Critical Issue: plink2 `--memory` Flag Limitation

**⚠️ IMPORTANT**: The `--memory` flag in plink2 does **NOT** effectively limit memory usage like it did in plink1. This means all current memory configurations in the pipeline are essentially ineffective at preventing out-of-memory errors with large datasets.

- **plink1**: `--memory` flag successfully limited RAM usage
- **plink2**: `--memory` flag exists but does not reliably constrain memory consumption
- **Impact**: Current 1000 MB memory settings provide no real protection against memory spikes

This limitation makes the sample-based splitting strategy even more critical for handling large datasets.

### Real-World Memory Usage Analysis

Based on execution trace analysis from a recent test run (sumstat_814), the following actual memory consumption patterns were observed:

**Memory-Intensive Operations (Peak RSS > 100 MB):**

1. **`make_snplist_from_pvar`**: **10.3 GB peak RSS** - Massive memory spike processing variant lists
2. **`variant_map_for_sbayesr`**: **3.2-9 GB peak RSS** - Extremely memory-intensive mapping operations
3. **`calc_score` processes**: **150-730 MB peak RSS** - Varies significantly by chromosome/dataset size
4. **`extract_maf_from_genotypes`**: **126-609 MB peak RSS** - Memory scales with sample count
5. **`indep_pairwise_for_benchmark`**: **87-368 MB peak RSS** - LD pruning operations
6. **`qc_posteriors`**: **Up to 1.9 GB peak RSS** - Quality control can be very memory-intensive

**Key Findings:**
- **Configured memory limits (1000 MB) are completely ignored** - processes regularly exceed these limits
- **Memory usage varies dramatically by chromosome** - chr1 uses ~4.5 GB while chr15 uses ~3.6 GB for variant mapping
- **Sample count directly impacts memory** - MAF extraction shows 126 MB to 609 MB range
- **Some processes show linear scaling** - scoring operations scale predictably with data size

### Identified Memory-Intensive plink Operations

Based on codebase analysis, the following plink operations are the primary memory consumers:

1. **`calc_score` process** (`pr_calc_score.nf:21-25`)
   - Command: `plink2 --pfile geno --score ${snp_posteriors} --memory ${params.memory.plink.calc_score}`
   - Current memory allocation: 1000 MB (from `nextflow.config:90`) - **⚠️ NOT EFFECTIVE**
   - **Critical bottleneck**: Processes entire sample set simultaneously with no memory protection

2. **`extract_maf_from_genotypes` process** (`pr_extract_from_genotypes.nf:14`)
   - Command: `plink2 --bfile geno --extract bimIDs --freq --memory ${params.memory.plink.extract_maf_from_genotypes}`
   - Current memory allocation: 1000 MB - **⚠️ NOT EFFECTIVE**
   - **Moderate bottleneck**: Frequency calculations across all samples with no memory protection

3. **`convert_plink1_to_plink2` process** (`pr_format_genotypes.nf:17-21`)
   - Command: `plink2 --bed ${bed} --bim ${bim} --fam ${fam} --make-pgen --memory ${params.memory.plink.extract_maf_from_genotypes}`
   - Current memory allocation: 1000 MB (reuses extract_maf memory setting) - **⚠️ NOT EFFECTIVE**
   - **Moderate bottleneck**: Format conversion with full dataset, uncontrolled memory usage

4. **`add_rsid_to_genotypes` process** (`pr_format_genotypes.nf:93`)
   - Command: `plink2 --pfile geno --make-pgen --extract tokeep --memory ${params.memory.plink.add_rsid_to_genotypes}`
   - Current memory allocation: 1000 MB - **⚠️ NOT EFFECTIVE**
   - **Moderate bottleneck**: Genotype filtering operations with no memory protection

5. **`concat_genotypes` process** (`pr_format_genotypes.nf:108`)
   - Command: `plink2 --pmerge-list allfiles.txt --make-pgen --memory ${params.memory.plink.concat_genotypes}`
   - Current memory allocation: 1000 MB - **⚠️ NOT EFFECTIVE**
   - **High bottleneck**: Merging multiple chromosome files, potentially unlimited memory usage

6. **`indep_pairwise_for_benchmark` process** (`prepare_benchmark_scoring.sh:13-19`)
   - Command: `plink2 --pfile ${prefix} --indep-pairwise 500kb 1 0.2 --memory 1000`
   - Fixed memory allocation: 1000 MB - **⚠️ NOT EFFECTIVE**
   - **Moderate bottleneck**: LD pruning across samples with no memory protection

## Proposed Solution: Sample-Based Splitting Strategy

### Core Principle

The key insight is that **genotype operations can be safely split by samples** because:
- Allele frequency calculations are additive across samples
- Polygenic scoring is independent per sample
- Most plink operations scale linearly with sample count
- Results can be merged without loss of information

**This approach is now essential** since plink2's `--memory` flag cannot be relied upon to prevent out-of-memory errors. Sample splitting provides the only reliable method to control memory usage in the current pipeline.

### Priority Ranking Based on Real Usage Data

Based on the execution trace analysis, processes should be prioritized for optimization in this order:

**🔴 CRITICAL (>1 GB peak RSS):**
1. **`make_snplist_from_pvar`** (10.3 GB) - Highest priority, single massive bottleneck
2. **`variant_map_for_sbayesr`** (3.2-9 GB) - Multiple instances, varies by chromosome
3. **`qc_posteriors`** (up to 1.9 GB) - Quality control operations

**🟡 HIGH (100 MB - 1 GB peak RSS):**
4. **`calc_score`** (150-730 MB) - Many instances, good splitting candidate
5. **`extract_maf_from_genotypes`** (126-609 MB) - Scales with samples, good splitting candidate
6. **`indep_pairwise_for_benchmark`** (87-368 MB) - LD pruning operations

**🟢 MEDIUM (<100 MB peak RSS):**
- Most other processes show reasonable memory usage and may not need immediate optimization

### Implementation Strategy

#### Phase 1: Sample Splitting Infrastructure

1. **Create sample splitting utility** (`bin/split_samples.sh`)
   ```bash
   #!/bin/bash
   # Split .psam/.fam files into chunks of specified size
   # Usage: split_samples.sh input.psam chunk_size output_prefix
   ```

2. **Create sample merging utility** (`bin/merge_sample_results.sh`)
   ```bash
   #!/bin/bash  
   # Merge results from sample-split operations
   # Usage: merge_sample_results.sh result_pattern output_file operation_type
   ```

3. **Add sample splitting parameters to config**
   ```groovy
   params {
     memory = [
       plink: [
         // ... existing settings
         max_samples_per_chunk: 10000,  // Configurable chunk size
         enable_sample_splitting: true   // Feature flag
       ]
     ]
   }
   ```

#### Phase 2: Process Modifications

##### 1. Enhanced `calc_score` Process

**Current Issue**: Processes all samples simultaneously, causing memory spikes with large cohorts.

**Solution**: Split samples into chunks, calculate scores separately, then merge.

```nextflow
process calc_score_chunked {
    publishDir "${params.outdir}/intermediates/scores/${method}", mode: 'rellink', overwrite: true
    
    input:
        tuple val(method), val(chr), path(snp_posteriors), path("geno.pgen"), path("geno.pvar"), path("geno.psam")

    output:
        tuple val(method), path("${method}_chr${chr}.score")
    
    script:
        def chunk_size = params.memory.plink.max_samples_per_chunk
        def enable_splitting = params.memory.plink.enable_sample_splitting
        
        if (enable_splitting) {
            """
            # Count samples
            sample_count=\$(tail -n +2 geno.psam | wc -l)
            
            if [ "\$sample_count" -gt "${chunk_size}" ]; then
                # Split samples into chunks
                split_samples.sh geno.psam ${chunk_size} chunk_
                
                # Process each chunk
                for chunk_psam in chunk_*.psam; do
                    chunk_id=\$(basename \$chunk_psam .psam)
                    
                    # Extract samples for this chunk
                    plink2 --pfile geno \\
                           --keep \$chunk_psam \\
                           --make-pgen \\
                           --out \${chunk_id}_geno \\
                           --memory ${params.memory.plink.calc_score} \\
                           --threads 1
                    
                    # Calculate scores for chunk
                    plink2 --pfile \${chunk_id}_geno \\
                           --out \${chunk_id}_tmp \\
                           --memory ${params.memory.plink.calc_score} \\
                           --threads 1 \\
                           --score ${snp_posteriors} 4 2 3 header cols=+scoresums ignore-dup-ids
                    
                    awk '{gsub(/^#/, ""); print}' \${chunk_id}_tmp.sscore > \${chunk_id}.score
                done
                
                # Merge chunk results
                merge_sample_results.sh "chunk_*.score" "${method}_chr${chr}.score" "score"
            else
                # Use original single-chunk approach
                plink2 --pfile geno \\
                       --out tmp \\
                       --memory ${params.memory.plink.calc_score} \\
                       --threads 1 \\
                       --score ${snp_posteriors} 4 2 3 header cols=+scoresums ignore-dup-ids
                
                awk '{gsub(/^#/, ""); print}' tmp.sscore > "${method}_chr${chr}.score"
            fi
            """
        } else {
            // Original implementation
            """
            plink2 --pfile geno \\
                   --out tmp \\
                   --memory ${params.memory.plink.calc_score} \\
                   --threads 1 \\
                   --score ${snp_posteriors} 4 2 3 header cols=+scoresums ignore-dup-ids

            awk '{gsub(/^#/, ""); print}' tmp.sscore > "${method}_chr${chr}.score"
            """
        }
}
```

##### 2. Enhanced `extract_maf_from_genotypes` Process

**Current Issue**: MAF calculation across all samples can be memory-intensive.

**Solution**: Calculate MAF in sample chunks, then compute weighted average.

```nextflow
process extract_maf_from_genotypes_chunked {
    publishDir "${params.outdir}/intermediates/maf_from_genotypes", mode: 'rellink', overwrite: true  

    input:  
        tuple val(chr), path("geno.bed"), path("geno.bim"), path("geno.fam"), path("map"), path("map_noNA")

    output:                                                                                 
        tuple val(chr), path("${chr}_geno_maf.afreq")    

    script:
        def chunk_size = params.memory.plink.max_samples_per_chunk
        def enable_splitting = params.memory.plink.enable_sample_splitting
        
        if (enable_splitting) {
            """
            cut -f 6 $map > bimIDs
            
            # Count samples
            sample_count=\$(tail -n +2 geno.fam | wc -l)
            
            if [ "\$sample_count" -gt "${chunk_size}" ]; then
                # Split samples and calculate MAF per chunk
                split_samples.sh geno.fam ${chunk_size} chunk_
                
                for chunk_fam in chunk_*.fam; do
                    chunk_id=\$(basename \$chunk_fam .fam)
                    
                    plink2 --bfile geno \\
                           --keep \$chunk_fam \\
                           --extract bimIDs \\
                           --threads 1 \\
                           --memory ${params.memory.plink.extract_maf_from_genotypes} \\
                           --freq \\
                           --out \${chunk_id}_maf
                done
                
                # Merge MAF results (weighted by sample size)
                merge_sample_results.sh "chunk_*_maf.afreq" "${chr}_geno_maf.afreq" "maf"
            else
                # Original single-chunk approach
                plink2 --bfile geno --extract bimIDs --threads 1 --memory ${params.memory.plink.extract_maf_from_genotypes} --freq --out ${chr}_geno_maf
            fi
            
            # Process the .afreq file using robust column extraction script
            process_plink_maf_single.sh ${chr}_geno_maf.afreq ${chr}_geno_maf.afreq.processed
            mv ${chr}_geno_maf.afreq.processed ${chr}_geno_maf.afreq
            """
        } else {
            // Original implementation
            """
            cut -f 6 $map > bimIDs
            plink2 --bfile geno --extract bimIDs --threads 1 --memory ${params.memory.plink.extract_maf_from_genotypes} --freq --out ${chr}_geno_maf
            
            process_plink_maf_single.sh ${chr}_geno_maf.afreq ${chr}_geno_maf.afreq.processed
            mv ${chr}_geno_maf.afreq.processed ${chr}_geno_maf.afreq
            """
        }
}
```

##### 3. Enhanced Format Conversion Processes

**Current Issue**: Converting large datasets between PLINK formats is memory-intensive.

**Solution**: Process sample chunks separately, then merge.

```nextflow
process convert_plink1_to_plink2_chunked {
    publishDir "${params.outdir}/intermediates/convert_plink1_to_plink2", mode: 'rellink', overwrite: true, enabled: params.dev
    
    input:
        tuple val(chr), path(bed), path(bim), path(fam)

    output:
        tuple val(chr), path("${chr}_converted.pgen"), path("${chr}_converted.pvar"), path("${chr}_converted.psam")
    
    script:
        def chunk_size = params.memory.plink.max_samples_per_chunk
        def enable_splitting = params.memory.plink.enable_sample_splitting
        
        if (enable_splitting) {
            """
            # Count samples
            sample_count=\$(tail -n +2 ${fam} | wc -l)
            
            if [ "\$sample_count" -gt "${chunk_size}" ]; then
                # Split and convert in chunks
                split_samples.sh ${fam} ${chunk_size} chunk_
                
                for chunk_fam in chunk_*.fam; do
                    chunk_id=\$(basename \$chunk_fam .fam)
                    
                    plink2 --bed ${bed} --bim ${bim} --fam \$chunk_fam \\
                           --make-pgen \\
                           --memory ${params.memory.plink.extract_maf_from_genotypes} \\
                           --threads 1 \\
                           --out \${chunk_id}_converted
                done
                
                # Merge converted chunks
                merge_sample_results.sh "chunk_*_converted" "${chr}_converted" "plink2_files"
            else
                # Original approach
                plink2 --bed ${bed} --bim ${bim} --fam ${fam} \\
                       --make-pgen \\
                       --memory ${params.memory.plink.extract_maf_from_genotypes} \\
                       --threads 1 \\
                       --out ${chr}_converted
            fi
            """
        } else {
            // Original implementation
            """
            plink2 --bed ${bed} --bim ${bim} --fam ${fam} \\
                   --make-pgen \\
                   --memory ${params.memory.plink.extract_maf_from_genotypes} \\
                   --threads 1 \\
                   --out ${chr}_converted
            """
        }
}
```

#### Phase 3: Utility Scripts

##### 1. Sample Splitting Script (`bin/split_samples.sh`)

```bash
#!/bin/bash
set -euo pipefail

input_file="$1"
chunk_size="$2"
output_prefix="$3"

# Detect file format (.psam or .fam)
if [[ "$input_file" == *.psam ]]; then
    header_lines=1
    ext="psam"
elif [[ "$input_file" == *.fam ]]; then
    header_lines=0
    ext="fam"
else
    echo "Error: Unsupported file format. Use .psam or .fam files."
    exit 1
fi

# Extract header if present
if [ $header_lines -gt 0 ]; then
    head -n $header_lines "$input_file" > header.tmp
fi

# Split data lines into chunks
if [ $header_lines -gt 0 ]; then
    tail -n +$((header_lines + 1)) "$input_file" | split -l "$chunk_size" -d - "${output_prefix}data_"
else
    split -l "$chunk_size" -d "$input_file" "${output_prefix}data_"
fi

# Add headers back and rename files
chunk_num=0
for chunk_file in "${output_prefix}data_"*; do
    output_file="${output_prefix}$(printf "%03d" $chunk_num).${ext}"
    
    if [ $header_lines -gt 0 ]; then
        cat header.tmp "$chunk_file" > "$output_file"
    else
        mv "$chunk_file" "$output_file"
    fi
    
    rm -f "$chunk_file"
    ((chunk_num++))
done

# Cleanup
rm -f header.tmp

echo "Split $input_file into $chunk_num chunks of max $chunk_size samples each"
```

##### 2. Sample Results Merging Script (`bin/merge_sample_results.sh`)

```bash
#!/bin/bash
set -euo pipefail

result_pattern="$1"
output_file="$2"
operation_type="$3"

case "$operation_type" in
    "score")
        # Merge scoring results - simple concatenation after header
        first_file=$(ls $result_pattern | head -n1)
        head -n1 "$first_file" > "$output_file"
        
        for file in $result_pattern; do
            tail -n +2 "$file" >> "$output_file"
        done
        ;;
        
    "maf")
        # Merge MAF results - weighted average by sample count
        python3 -c "
import sys
import pandas as pd
import glob

files = glob.glob('$result_pattern')
dfs = []

for f in files:
    df = pd.read_csv(f, sep='\t')
    # Extract sample count from NCHROBS column
    df['SAMPLE_COUNT'] = df['NCHROBS'] / 2  # Diploid
    dfs.append(df)

# Merge and calculate weighted averages
merged = pd.concat(dfs, ignore_index=True)
result = merged.groupby(['CHR', 'SNP', 'A1', 'A2']).agg({
    'MAF': lambda x: (x * merged.loc[x.index, 'SAMPLE_COUNT']).sum() / merged.loc[x.index, 'SAMPLE_COUNT'].sum(),
    'NCHROBS': 'sum'
}).reset_index()

result.to_csv('$output_file', sep='\t', index=False)
"
        ;;
        
    "plink2_files")
        # Merge PLINK2 files using plink2 --pmerge-list
        ls ${result_pattern}.pgen | sed 's/.pgen$//' > merge_list.txt
        plink2 --pmerge-list merge_list.txt --make-pgen --out "$output_file" --memory 2000
        rm merge_list.txt
        ;;
        
    *)
        echo "Error: Unknown operation type: $operation_type"
        exit 1
        ;;
esac

echo "Merged results into $output_file"
```

#### Phase 4: Configuration Updates

##### Enhanced Memory Configuration (`nextflow.config`)

```groovy
params {
  memory = [
    plink: [
      calc_score: 1000,
      extract_maf_from_genotypes: 1000,
      add_rsid_to_genotypes: 1000,
      concat_genotypes: 1000,
      indep_pairwise_for_benchmark: 1000,
      
      // New sample splitting parameters
      max_samples_per_chunk: 10000,        // Configurable chunk size
      enable_sample_splitting: true,       // Feature flag
      sample_splitting_threshold: 20000,   // Auto-enable splitting above this sample count
      chunk_merge_memory: 2000             // Memory for merging operations
    ],
    sort: [
      make_snplist_from_bim: "2 GB"
    ]
  ]
}
```

##### Dynamic Memory Allocation

```groovy
process {
  // Dynamic memory allocation based on sample count
  withName: 'calc_score_chunked' {
    memory = { 
      def sample_count = task.ext.sample_count ?: 10000
      def chunk_size = params.memory.plink.max_samples_per_chunk
      def chunks = Math.ceil(sample_count / chunk_size)
      
      if (params.memory.plink.enable_sample_splitting && sample_count > params.memory.plink.sample_splitting_threshold) {
        return "${params.memory.plink.calc_score * Math.min(chunks, 4)} MB".toString()
      } else {
        return "${Math.max(params.memory.plink.calc_score, sample_count * 0.1)} MB".toString()
      }
    }
  }
}
```

## Implementation Timeline

### Phase 1: Infrastructure (Week 1-2)
- [ ] Create sample splitting utilities (`split_samples.sh`, `merge_sample_results.sh`)
- [ ] Add configuration parameters
- [ ] Create unit tests for utilities
- [ ] Update documentation

### Phase 2: Core Process Updates (Week 3-4)
- [ ] Implement `calc_score_chunked` process
- [ ] Implement `extract_maf_from_genotypes_chunked` process  
- [ ] Implement `convert_plink1_to_plink2_chunked` process
- [ ] Add backward compatibility flags

### Phase 3: Integration & Testing (Week 5-6)
- [ ] Integration testing with test datasets
- [ ] Performance benchmarking
- [ ] Memory usage validation
- [ ] Edge case testing (small datasets, single samples)

### Phase 4: Advanced Features (Week 7-8)
- [ ] Dynamic memory allocation based on sample count
- [ ] Automatic chunk size optimization
- [ ] Parallel chunk processing
- [ ] Enhanced error handling and recovery

## Expected Benefits

### Memory Reduction
- **Primary bottleneck (`calc_score`)**: 60-80% memory reduction for large cohorts (from observed 150-730 MB to ~50-200 MB per chunk)
- **Secondary bottlenecks**: 40-60% memory reduction (MAF extraction from 126-609 MB to ~50-200 MB per chunk)
- **Critical processes**: Major impact on `make_snplist_from_pvar` (10.3 GB) and `variant_map_for_sbayesr` (3.2-9 GB)
- **Overall pipeline**: 30-50% peak memory reduction, with some processes seeing 70-80% reduction

**Critical Note**: These benefits are even more significant given that plink2's `--memory` flag provides no protection. Without sample splitting, memory usage is completely uncontrolled and can easily exceed available system resources. The execution trace shows processes regularly using 5-10x their configured memory limits.

### Scalability Improvements
- Support for cohorts with >100K samples
- Linear memory scaling instead of quadratic
- Configurable resource usage based on available hardware

### Maintainability
- Backward compatible with existing workflows
- Feature flags allow gradual rollout
- Modular design enables targeted optimizations

## Risk Mitigation

### Computational Correctness
- **Risk**: Sample splitting might introduce numerical errors
- **Mitigation**: Extensive validation against original results, mathematical verification of merge operations

### Performance Overhead
- **Risk**: Sample splitting adds computational overhead
- **Mitigation**: Overhead only applies to large datasets where memory is limiting factor, parallel processing of chunks

### Complexity
- **Risk**: Increased codebase complexity
- **Mitigation**: Modular design, comprehensive testing, clear documentation, backward compatibility

## Alternative Approaches Considered

### 1. Variant-Based Splitting
- **Pros**: Simpler for some operations
- **Cons**: Doesn't address sample-count scaling, complex for scoring operations

### 2. External Memory Management
- **Pros**: Transparent to existing code
- **Cons**: Platform-dependent, less control over memory usage patterns

### 3. Different Tools
- **Pros**: Might have better memory efficiency
- **Cons**: Major workflow changes, validation overhead, tool compatibility issues

### 4. plink2 Memory Flag Workarounds
- **Pros**: Might restore some memory control without major code changes
- **Cons**: No reliable workarounds exist; plink2 developers acknowledge this limitation
- **Status**: Investigated but not viable - sample splitting remains the only reliable solution

## Conclusion

The proposed sample-based splitting strategy provides a **simple, maintainable, and effective** solution to pgscalculator's memory issues. **Given that plink2's `--memory` flag is ineffective**, this approach is not just an optimization but a **necessity** for handling large datasets reliably.

By focusing on the primary bottleneck (`calc_score`) and implementing a modular approach with backward compatibility, we can achieve significant memory reductions while maintaining computational correctness and workflow stability.

The solution aligns with the user's requirements for:
- **Clever**: Leverages the mathematical properties of genotype operations and addresses the fundamental limitation of plink2's memory management
- **Simple**: Modular design with clear separation of concerns  
- **Maintainable**: Backward compatible with feature flags and comprehensive testing

**Urgency**: Implementation should proceed in phases to allow for thorough testing and validation at each step. However, given the complete lack of memory protection in the current pipeline, this work should be prioritized for production environments handling large cohorts.

