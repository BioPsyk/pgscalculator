# SBayesRC Performance and Computational Requirements

This document provides performance metrics and computational requirements for running SBayesRC in the pgscalculator pipeline.

## Overview

SBayesRC is a method that incorporates functional genomic annotations with high-density SNPs (> 7 million) for polygenic prediction. The computational requirements vary significantly based on the annotation data and LD reference used.

*Source: [SBayesRC GitHub Repository](https://github.com/zhilizheng/SBayesRC)*

## Annotation Data Options

### 1. Baseline Model 2.2 (Recommended)
- **SNPs**: ~8 million SNPs
- **File size**: ~293MB compressed, ~2GB uncompressed
- **Description**: Functional annotation information from baseline model 2.2
- **Reference**: Márquez-Luna 2021, DOI: 10.1038/s41467-021-25171-9
- **Download**: Available from SBayesRC resources

### 2. Custom Annotations
- Users can provide custom functional annotations
- Must match the SNP set used in LD reference
- Format requirements detailed in SBayesRC documentation

## LD Reference Options

### HapMap3 SNPs (Faster, Smaller)
- **SNPs**: ~1 million SNPs
- **File size**: ~2.9GB compressed
- **Computational time**: Significantly faster
- **Memory requirements**: Lower
- **Use case**: Quick testing, proof of concept, resource-constrained environments

### Imputed SNPs (Full Resolution)
- **SNPs**: ~7 million SNPs  
- **File size**: Much larger (exact size varies by ancestry)
- **Computational time**: Much longer
- **Memory requirements**: Higher
- **Use case**: Production runs, maximum accuracy

## Computational Requirements

### Minimum Requirements
- **CPU cores**: 6 cores minimum
- **Memory**: 10GB minimum
- **Storage**: 50GB+ for data and intermediate files
- **Time**: 1-4 hours (HapMap3), 8-24 hours (Imputed)

### Recommended Requirements  
- **CPU cores**: 22 cores
- **Memory**: 40GB
- **Storage**: 100GB+ 
- **Time**: 30 minutes - 2 hours (HapMap3), 2-8 hours (Imputed)

### SLURM Resource Allocation
Based on successful runs in this pipeline:
```bash
# Recommended SLURM parameters
--mem=40g 
--ntasks=1 
--cpus-per-task=22 
--time=1:00:00  # HapMap3
--time=4:00:00  # Imputed (adjust based on data size)
```

## Performance Comparison

| LD Reference | SNPs | Typical Runtime* | Memory Usage | Accuracy |
|--------------|------|------------------|--------------|----------|
| HapMap3      | ~1M  | 30min - 2h      | 10-20GB     | Good     |
| Imputed      | ~7M  | 2-8h            | 20-40GB     | Best     |

*Runtime varies significantly based on:
- Number of CPU cores
- Memory availability  
- I/O performance
- Summary statistics size
- Annotation complexity

## System Compatibility

### Tested Operating Systems
- CentOS > 7
- Debian > 9  
- Ubuntu > 20.04
- macOS > 11

### Known Issues
- **Ubuntu 22.04**: Has broken openBLAS version 0.3.20 affecting matrix operations
- **CentOS 7**: Default R may have issues with RcppEigen package (gcc 4.8)
- **Solution**: Use R from anaconda for problematic systems

## Resource Planning

### For Testing/Development
- Use **HapMap3 LD reference** 
- Expect **30 minutes to 2 hours** runtime
- Allocate **20GB memory, 6+ cores**

### For Production
- Use **Imputed LD reference** for maximum accuracy
- Expect **2-8 hours** runtime  
- Allocate **40GB memory, 22 cores**
- Plan for **100GB+ storage**

### Download Times
- **Annotations**: ~30 minutes (293MB)
- **HapMap3 LD**: ~2-4 hours (2.9GB) 
- **Imputed LD**: ~8-12 hours (varies by ancestry)

## Optimization Tips

1. **Use screen/tmux** for long-running downloads and jobs
2. **Pre-download data** during off-peak hours
3. **Use HapMap3** for initial testing and validation
4. **Scale up to Imputed** only for final production runs
5. **Monitor memory usage** - increase if jobs fail with OOM errors
6. **Use parallel processing** where supported (eigen decomposition)

## Container Performance

The pgscalculator uses containerized SBayesRC which:
- Ensures consistent environment across systems
- Includes optimized BLAS libraries
- Handles dependency management automatically
- May have slight performance overhead vs native installation

## References

- SBayesRC Method: Zheng Z, Liu S, Sidorenko J, et al. (2024) Nature Genetics, doi:10.1038/s41588-024-01704-y
- Implementation: [SBayesRC GitHub Repository](https://github.com/zhilizheng/SBayesRC)
- Resources: [SBayesRC Data Hub](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/)
