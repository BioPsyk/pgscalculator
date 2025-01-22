# Frequently Asked Questions

## Configuration Files

### Q: There are two config files, which one is used?

The configuration system in pgscalculator follows a hierarchical structure:

1. **Base Configuration**: The `nextflow.config` file serves as the base configuration, containing default settings for all environments.

2. **Method-Specific Configuration**: Files in the `conf/` directory (e.g., `sbayesr.config`, `prscs.config`) contain method-specific settings.

When running the pipeline, any parameters defined in your chosen method-specific config file will override the corresponding settings in the base `nextflow.config`. This allows for:

- Maintaining default settings in one central location
- Customizing specific parameters for different analysis methods
- Flexibility in configuration without modifying the base settings

For example, if `nextflow.config` sets `params.method="sbayesr"` but your chosen config file sets `params.method="prscs"`, the latter will take precedence. 

## Variant Selection

### Q: Is it possible to add a variant exclude file?

No, but you can achieve the same result using an include file, which will exclude everything not present in that file. Here's how to handle variant exclusion:

1. First, extract the list of variants from your genotype files (not summary statistics)
2. Remove the variants you want to exclude from this list 
3. Use the resulting list as an include file with the `--snplist` parameter

This approach effectively excludes unwanted variants by only including the ones you want to keep. It's more efficient as it allows the pipeline to process only the relevant variants from the start.

Important note: The `--snplist` parameter is designed to filter variants in the genotype files only. The variant IDs in your include list must match those in your genotype files. Currently, this feature does not filter variants in the summary statistics file - this may be added as functionality in a future release.

## Analysis Methods

### Q: What is the purpose of the benchmark calculation?

The benchmark calculation serves as a reference point to evaluate the performance of more sophisticated PRS methods. Here's why it's important:

1. **Baseline Comparison**: Simple pruning and thresholding (P+T) serves as a baseline method that is:
   - Well-established
   - Computationally efficient
   - Easy to interpret

2. **Performance Evaluation**: By comparing advanced methods (like SBayesR) against this benchmark, you can:
   - Validate that more complex models provide meaningful improvements
   - Quantify the gain in predictive power
   - Justify the additional computational resources required

3. **Quality Control**: If an advanced method performs worse than the benchmark, it might indicate:
   - Potential issues with model parameters
   - Data quality problems
   - Need for further optimization

This comparison helps ensure that the additional complexity of advanced methods translates into tangible improvements in prediction accuracy.


