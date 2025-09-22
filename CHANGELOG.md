# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.2] - 2025-09-22
### Fixed
- **Critical sorting issue in variant_map_for_sbayesr.sh that caused incomplete variant mappings**
- Missing sort operation for snp2 file before first join operation, which led to join failures
- Incomplete or incorrect results when processing variant mapping for SBayesR workflow

### Added
- Comprehensive sorting validation tests to detect and prevent sorting issues
- Dedicated test suite for variant mapping sorting validation (`test_variant_map_sorting_validation.sh`)
- Helper functions for sorting validation in unit tests (`_check_file_sorted`, `_validate_join_prerequisites`)
- Detailed technical documentation of sorting issues and fixes
- Performance testing for large datasets (100+ variants) in sorting validation

### Changed
- Enhanced existing variant mapping tests with sorting validation capabilities
- Improved test coverage for edge cases (empty files, single records, different genome builds)

## [1.3.1] - 2025-09-13
### Added
- Robust column-name-based MAF extraction script for better maintainability
- PSAM file standardization to IID-only format to ensure consistent score calculations
- Comprehensive unit tests for MAF extraction functionality

### Fixed
- Zero scores in combined output due to FID/IID format inconsistencies between genotype files
- Incorrect MAF values and missing NCHROBS in raw_maf_chrall output
- Pipeline sensitivity to varying plink2 .afreq output column formats
- Variant ID mangling when using genotype files with missing rsIDs

### Changed
- MAF extraction now uses column names instead of positional indices for robustness
- Score calculation simplified after implementing consistent PSAM standardization

## [1.3.0] - 2025-06-29
### Added
- **PLINK2 dosage format support as default genotype input format**
- **Automatic PLINK1 to PLINK2 conversion with format detection**
- Multi-architecture support for Docker images (amd64 and arm64)
- Docker manifest support for seamless cross-platform deployment
- Enhanced genotype file processing with improved variant ID handling
- Comprehensive unit tests for genotype format conversion and processing

### Changed
- **Default genotype format changed from PLINK1 to PLINK2 dosage format**
- Updated Docker build and push scripts to handle multi-arch builds
- Improved Docker image distribution with platform-specific tags
- Enhanced variant mapping and benchmark scoring for PLINK2 format
- Improved debugging and channel tracing capabilities

## [1.2.7] - 2025-05-19
### Fixed
- Intermediate scores to be in method specific folders

## [1.2.6] - 2024-11-29
### Added
- Enhanced FAQ documentation with detailed explanation of benchmark calculations

## [1.2.5] - 2024-11-29
### Fixed

- renaming a duplicate header in augmentet output

## [1.2.4] - 2024-09-02
### Changed

- remove memory upper limits in nextflow.config. Replaced by setting plink and sort memory variables for each specific process in nextflow.config

### Fixed
- Remove the tmp in mount init script

## [1.2.3] - 2024-08-23
### Changed

- Removed memory restriction set in nextflow.config, and removed all labels from all processes

## [1.2.2] - 2024-08-16
### Changed

- README.md to keep only the absolutely most important to run using singularity
- Extra documentation moved to their own doc files in docs/

### Fixed
- Missing mount dir for supplied snplists in pgscalculator.sh

## [1.2.1] - 2024-04-30
### Added

- Memory constraints on plink2, so that the overall memory footprint will be lower.
- Checker for missing or duplicate variant ids in the genotype data, which will fill in with chr:pos:a1:a2 if missing, and add an extra number if still a duplicate.

## [1.2.0] - 2024-04-24
### Added

- b38 support and improved allele flip management. Also added an early filter on b37 NA coordinates.

## [1.1.0] - 2024-04-10
### Added

- The sbayesr workflow has been improved in all aspects: a companion variant map file, an augmentated sumstat file, and in general much better channeling and output structure.

## [1.0.0] - 2023-11-10
### Added

- Everything to run an sbayesr workflow.

## [0.1.0] - 2023-10-02
### Added

- Basic containerized setup, and a general reformatter for the different prs softwares

