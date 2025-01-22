# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

