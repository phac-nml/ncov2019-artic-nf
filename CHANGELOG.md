# phac-nml/ncov2019-artic-nf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [2024-04-25]
Overall version 2.0.0 has the all the same outputs as version 1.1.0 but with some adjustments to the output locations and the input parameter names. This makes this release incompatible with previous automation unfortunately but it is ultimately for a more robust  and easier to run/develop pipeline

### `Added`
- Nextclade for an additional validation to frameshifts
    - Nextclade specific args
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
- Output file locations
    - Output locations are still based on the process name. These have been deduplicated to make development and updates easier leading to a slight difference.
        - Ex. `articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish` to `articNcovNanopore_articMinION`
- Nextflow config utilization
    - Removed params from the `base.config` config and placed them in the main `nextflow.config`
    - `base.config` now for resource setting
- Ncov-tools config
    - Cleaned up ncov-tools config comments
    - Allow run name and primer names to be set by pipeline instead of being locked in by config

### `Removed`
- All Illumina steps and args
- All climb upload steps and args

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
