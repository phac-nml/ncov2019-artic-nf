# phac-nml/ncov2019-artic-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v3.0.0 - [2025-02-25]
Version 3.0.0 is a major, breaking release that updates to running `clair3` over `medaka` or `nanopolish` along with fixing some of the bugs that had been identified over the years.

`Added`:
- Modules:
    - `checkFastqForModel`: To check the fastq header for the model to use
    - `articDownloadModels`: To download `clair3` models
- Parameters:
    - `clair3_model`: Optional parameter to specify `clair3` model to use

`Changed`: 
- Updated artic `1.4.x` --> `1.6.1`
    - Changes variant calling to `clair3`
    - Fixes some of the bugs reported
        - BCFTools one specifically
- Formatting modules files to separate out processes a bit more
- Formatting environment files to remove the `defaults` branch

`Removed`:
- Medaka and Nanopolish parameters and pipeline
    - Parameters:
        - `--nanopolish`
        - `--medaka`
        - `--fast5_pass`
        - `--sequencing_summary`
        - `--medaka_model`
        - `--bwa`
        - `--no_longshot`

- Fast5 uploads script

## v2.0.2 - [2024-07-23]
Version 2.0.2 adjusts internal retry resources

## v2.0.1 - [2024-07-15]
Version 2.0.1 fixes a minor bug in which all numeric sample names (ex. 231301) were not being compared to the nextclade output for frameshift reporting and adjustments along with updating the nextclade dataset to default to 'latest'

### `Changed`:
- Bugfix in `qc.py` such to allow numeric sample names to be properly compared to the nextclade output for frameshift reporting
- Default nextclade dataset set to 'latest'

## v2.0.0 - [2024-04-25]
Overall version 2.0.0 has the all the same outputs as version 1.1.0 but with some adjustments to the output locations and the input parameter names. This makes this release incompatible with previous automation unfortunately but it is ultimately for a more robust  and easier to run/develop pipeline

### `Added`
- Nextclade for an additional validation to frameshifts
    - Run [nextclade](https://github.com/nextstrain/nextclade) with specified args to create TSV output
        - The TSV output is used to help validate frameshifts by looking at the consensus sequence and determining if there are any frameshifts and if they are known. Known frameshifts can be used to help correct the `qc_pass` column
    - Nextclade specific args added:
        - `--nextclade_dataset`: Name of nextclade dataset to pull. Default: 'sars-cov-2'
        - `--nextclade_tag`: Tag of nextclade dataset to pull. Default: '2024-04-15--15-08-22Z'
- New Params
    - Max resource args for some customizable resource management
    - `local_scheme` arg to make it easier to provide a local primer scheme
    - `--no_longshot` arg to allow skipping longshot in artic medaka pipeline
- Scheme validation workflow
    - Checks that the folders are found with the needed files
    - Creates amplicon bed file for ncov-tools
    - Detects primer names

### `Changed`
- Input param arguments
    - No longer both snake_case and camelCase args
    - All args have been set to snake_case
        - `--medakaModel` to `--medaka_model`
        - `--schemeRepoURL` to `--scheme_repo_url`
        - `--schemeVersion` to `--scheme_version`
        - `--minReadsPerBarcode` to `--min_reads_per_barcode`
        - `--minReadsGuppyPlex` to `--min_reads_guppyplex`
        - `--correctN` to `--skip_correct_n` and logic changes
        - `--sequencingTechnology` to `--sequencing_technology`
        - `--csqAfThreshold` to `--csq_af_threshold`
        - `--csqDpThreshold` to `--csq_dp_threshold`
- Output file locations
    - Output locations are still based on the process name. These have been deduplicated to make development and updates easier leading to a slight difference.
        - `articNcovNanopore_sequenceAnalysisMedaka_articDownloadScheme` to `articNcovNanopore_schemeValidate_validateScheme/`
        - `articNcovNanopore_sequenceAnalysisMedaka_articGuppyPlex` to `articNcovNanopore_articGuppyPlex`
        - `articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish` to `articNcovNanopore_articMinION`
        - `articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka` to `articNcovNanopore_articMinION`
        - `articNcovNanopore_sequenceAnalysisMedaka_renameSamples` to `articNcovNanopore_renameBarcodeSamples`
- Nextflow config utilization
    - Removed params from the `base.config` config and placed them in the main `nextflow.config`
    - `base.config` now for resource setting
- Ncov-tools config
    - Cleaned up ncov-tools config comments
    - Allow run name and primer names to be set by pipeline instead of being locked in by config

### `Removed`
- All Illumina steps and args
- All climb upload steps and args
- `--schemeDir` parameter
    - Now scheme directory is automatically adjusted to `primer-schemes`

### `Developer Changes`
- Test dataset fastq files were gzipped to save more space
- Tests that were not being run were removed
- Formatting/organization of nextflow workflow code was revisited and standardized
- Formatting/organization of most python code was revisited
- Removal of duplication of workflows and modules
    - Ex. `articMinionNanopolish` module and `articMinionMedaka` module consolidated to one `articMinion` module
    - Individual workflows for `--nanopolish` and `--medaka` were combined into one workflow as the steps were almost all the same
- Better passing of variables/channels to processes
- Updated channel creation and handling of empty channels
