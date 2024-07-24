# ncov2019-artic-nf
A Nextflow pipeline for running the ARTIC network's fieldbioinformatics tools (https://github.com/artic-network/fieldbioinformatics) on nanopore SARS-CoV-2 data to create pileups, variant calls, consensus sequences, plots and stats

### Introduction

------------

This Nextflow pipeline automates the ARTIC network [nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html "nCoV-2019 novel coronavirus bioinformatics protocol").

**This fork** rewrites the pipeline to reorganize the analysis steps to better follow more current nextflow practices along with optimizing its usage for the NML Canada's *nanopore* infrastructure. Changes include:
- Renaming input barcoded samples using a samplesheet
- Allowing flat fastq directory input when running with Medaka
- Better tracking of samples with no/too few reads that are filtered out
- Scheme input corrections/validation
- Bumping `artic` to v1.2.4
    - For Medaka version updates and models when running with Medaka
    - Default medaka version is specified as "r941_min_hac_g507"
- Running [ncov-tools](https://github.com/jts/ncov-tools) for plots, sequence statistics, and sequence quality status
- Addition of nextclade
- Additional steps and checks for to the final output QC file
- Automating the upload of output data to [IRIDA](https://github.com/phac-nml/irida) for storage
- Step tool version tracking

### Release Notes
For full changes visit the [CHANGELOG](CHANGELOG.md)

#### *v2.0.2
Version 2.0.2 adjusts internal retry resources

#### *v2.0.1*
Version 2.0.1 fixes a minor bug in which all numeric sample names (ex. 231301) were not being compared to the nextclade output for frameshift reporting and adjustments along with updating the nextclade dataset to 'latest'

`Changed`:
- Bugfix in `qc.py` such to allow numeric sample names to be properly compared to the nextclade output for frameshift reporting
- Default nextclade dataset set to 'latest'

#### *v2.0.0*
Overall version 2.0.0 has the all the same outputs as version 1.1.0 but with some adjustments to the output locations and the input parameter names. This makes this release incompatible with previous automation unfortunately but it is ultimately for a more robust  and easier to run/develop pipeline

`Added`:
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

`Changed`:
- Input param arguments
    - No longer both snake_case and camelCase args
    - All args have been set to snake_case
- Output file locations
    - Output locations are still based on the process name. These have been deduplicated to make development and updates easier leading to a slight difference.
        - Ex. `articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish` to `articNcovNanopore_articMinION`
- Nextflow config utilization

`Removed`:
- All Illumina steps and args
- All climb upload steps and args

### Installation

------------

An up-to-date version of Nextflow (`version >=22.10.0`) and Conda is required to run the pipeline.

You can either:
1. Download and install [conda](https://docs.conda.io/en/latest/miniconda.html) and then install Nextflow and Mamba using conda (Recommended)
    - Conda command: `conda create -n nextflow -c conda-forge -c bioconda nextflow mamba`
2. Follow the instructions at https://www.nextflow.io/ to download and install Nextflow
    - You will still have to install other dependencies after

Nextflow controls the execution of the pipeline steps while conda controls the installation of tools/depedencies required to run the pipeline

### Running Nanopore Data

------------

*Note 1*: Only `conda` is supported for dependency installation/containerization. It is required to have conda installed and use that when running this pipeline

*Note 2:* When running this pipeline it is best to specify your conda cache directory with `--cache 'path/to/cacheDir'` as generation of the ncov-tools env can be extremly slow even with mamba. All environments should resolve with no issue provided there are enough resources for the main nextflow process to do so but if not, you will get an error that displays the command and env that failed. Using that command, you can manually create the env which can be reused if giving the cache directory.

#### **Method 1 - Nanopolish with Barcoded Reads**

*Running*

Base Minimum Command:
```bash
nextflow run phac-nml/ncov2019-artic-nf \
    -profile conda \
    --nanopolish \
    --prefix "prefix_str" \
    --basecalled_fastq /path/to/fastq_pass_dir \
    --fast5_pass /path/to/fast5_pass_dir \ 
    --sequencing_summary /path/to/sequencing_summary.txt
```
This command will run on the input fastq directories using the default `freed` primer scheme to create output data named as `prefix_barcode##`.

Recommended Command:
```bash
nextflow run phac-nml/ncov2019-artic-nf \
    -profile conda \
    --cache /path/to/conda_cache_dir/ \
    --nanopolish \
    --prefix "prefix_str" \
    --basecalled_fastq /path/to/fastq_pass_dir \
    --fast5_pass /path/to/fast5_pass_dir \
    --sequencing_summary /path/to/sequencing_summary.txt \
    --irida /path/to/samplesheet_with_names.tsv \
    --scheme_version "Wanted_Scheme" \
    --min_length Min_Scheme_Length_INT \
    --max_length Max_Scheme_Length_INT
```
This command will run on the input fastq directories using the provided scheme and filtering size options. It will rename the inputs based on the samplesheet to allow better data tracking. Input metadata will be added to the final results

------------

#### **Method 2 - Medaka with Barcoded Reads**

*Running*

Basic command:
```bash
nextflow run phac-nml/ncov2019-artic-nf \
    -profile conda \
    --medaka \
    --medaka_model "Model" \
    --prefix "prefix_str" \
    --basecalled_fastq /path/to/fastq_pass_dir
```
This command will run medaka on the input fastq_pass_dir and create the same outputs as the nanopolish side. Files will still be output named as `prefix_barcode##`

Recommended Command:
```bash
nextflow run phac-nml/ncov2019-artic-nf \
    -profile conda \
    --cache /path/to/conda_cache_dir/ \
    --medaka \
    --medaka_model "Model" \
    --prefix "prefix_str" \
    --basecalled_fastq /path/to/fastq_pass_dir \
    --irida /path/to/samplesheet_with_names.tsv \
    --scheme_version "Wanted_Scheme" \
    --min_length Min_Scheme_Length_INT \
    --max_length Max_Scheme_Length_INT
```
This command will run on the input fastq directories using the provided scheme and filtering size options. It will rename the inputs based on the samplesheet to allow better data tracking. Input metadata will be added to the final results

*EXTREMELY IMPORTANT NOTE:* The medaka model should be set to the one that most closely matches your data otherwise the variant calls may be inaccurate. More info available here: https://github.com/nanoporetech/medaka#models. The default model is currently "r941_min_hac_g507" which should work for anyone still running the R9.4.1 flowcells

#### **Method 3 - Medaka with Flat Directory Input**

*Running*

Basic command:
```bash
nextflow run phac-nml/ncov2019-artic-nf \
    -profile conda \
    --medaka \
    --medaka_model "Model"
    --prefix "prefix_str" \
    --basecalled_fastq /path/to/flat_fastq_dir
```
This command will run medaka on the input flat_fastq_dir and create the same outputs as the other running methods. These files will retain their basename from the input directory.

Recommended Command:
```bash
nextflow run phac-nml/ncov2019-artic-nf \
    -profile conda \
    --cache /path/to/conda_cache_dir/ \
    --medaka \
    --medaka_model "Model" \
    --prefix "prefix_str" \
    --basecalled_fastq /path/to/fastq_pass_dir \
    --irida /path/to/samplesheet_with_names.tsv \
    --scheme_version "Wanted_Scheme" \
    --min_length Min_Scheme_Length_INT \
    --max_length Max_Scheme_Length_INT
```
This command will run on the input flat directory of fastq files using the provided scheme and filtering size options and add the input metadata to the final results.

*EXTREMELY IMPORTANT NOTE:* The medaka model should be set to the one that most closely matches your data otherwise the variant calls may be inaccurate. More info available here: https://github.com/nanoporetech/medaka#models. The default model is currently "r941_min_hac_g507" which should work for anyone still running the R9.4.1 flowcells

#### **Optional Args**
Arguments available for different aspects of the pipeline. The mandatory args are found above in the example commands

#### Help
Run `nextflow run phac-nml/ncov2019-artic-nf --help` to get a list of all the supported pipeline arguments

#### Schemes
`--scheme_repo_url <STR>`: Url to the repo that you want to pull schemes from. Default: 'https://github.com/phac-nml/primer-schemes.git'

`--scheme_version <STR>`: Scheme version to use for the pipeline. Default: 'freed'

`--scheme <STR>`: Name of the virus scheme to use. Default: 'nCoV-2019'

`--local_scheme <PATH>`: Path to directory formatted with a local scheme

Scheme directory formatting looks as such:
```
primer-schemes
└── <SCHEME>
    └── <SCHEME VERSION>
        ├── *.primer.bed
        ├── *.reference.fasta
        └── *.scheme.bed
```

#### Read Filtering
`--min_reads_per_barcode <INT>`: Minimum number of reads required for input barcode directory or flat fastq file to be run in the pipeline. Default: 100

`--min_reads_guppyplex <INT>`: Minimum number of guppyplex size selected reads required to be run in subsequent steps. Default: 10

`--min_length <INT>`: Minimum read length that should be set based on amplicon size. Default: 800

`--max_length <INT>`: Maximum read length that should be set based on amplicon size. Default: 1600

#### Artic Minion
`--normalise <INT>`: Normalize input reads for each amplicon to this value. Default: 500

`--bwa`: Run BWA to align reads to reference sequence instead of minimap2. Default: False

`--no_frameshift`: Add in VCF filter --no-frameshift argument that checks for %3==0 allele calls. Default: False

`--no_longshot`: Run medaka variant over longshot for variant analysis in the medaka pipeline. Default: False

#### Metadata
`--irida <PATH>`: Path to IRIDA samplesheet.tsv file to be used to rename samples and assist in IRIDA data uploads. The required columns are `sample`, `run`, and `barcode`. For uploads a `project_id` column is also required.

Below is an example samplesheet:

| sample | run | barcode | project_id | ct | date |
|-|-|-|-|-|-|
| Sample_name1 | run_name | 14 | 1422 | 23.33 | 2020-08-22 |
| Sample_name2 | run_name | 15 | 1422 | 22.53 | 2020-08-22 |
| Sample_name3 | run_name | 45 | 1422 | NA | NA |
| Sample_name6 | run_name | 65 | 1422 | NA | 2020-08-22 |

#### Uploads
`--upload_irida <PATH>`: Path to IRIDA instance config file to which the data will be uploaded to. Must be set with a --irida samplesheet.tsv file as well to know which project each sample goes to.

Format:
```conf
[Settings]
client_id = uploader
client_secret = <secret from IRIDA>
username = test
password = unsecure_password
base_url = https://<your_irida_instance>/irida/api
parser = directory
```

#### Nextclade
Nextclade is currently run and output as a separate results TSV along with being used to confirm SnpEFF (from ncov-tools) frameshift calls

`--nextclade_dataset <STR>`: Name of nextclade dataset to use. Default: 'sars-cov-2'

`--nextclade_tag <STR>`: Nextclade dataset tag to use. Default: '2024-04-15--15-08-22Z'

#### QC
`--sequencing_technology <STR>`: Name of the sequencing technology used. Default: Nanopore

`--pcr_primers <PATH>`: Path to qPCR primers to check for any distuptive mutations in. Default: 'baseDir/extra_data/pcr_primers.bed'

`--skip_correct_n`: Skip running custom python/bcftools N correction step

The correctN program double checks the failed vcf variant file positions (masked as Ns) to see if there sufficient evidence in the pileup to call a reference base at that location instead. This is written separately

#### **Negative Controls**

Negative controls can be used if `--irida` is given and one of the sample names has any of the following within it `negative`, `water`, `blank`, `ntc` in any capitalizations or spot.

This will run the negative control module of ncov-tools to check if there is any contamination in the run

*Negative Controls Limitations*

Connecting ncov-tools to the pipeline means that the ncov-tools config file has to be mostly set before running. In an attempt to make automation easier, the negative control names are set based on the key words mentioned (`negative`, `water`, `blank`, `ntc`) and passed by the samplesheet before. In the future, we may set-up a flag where you can set your negative control names but for now this is the easiest way to get all our users to run the pipeline.

---------------

### Executors
By default, the pipeline just runs on the local machine. You can specify `-profile slurm` to use a SLURM cluster, or `-profile lsf` to use an LSF cluster. It is recommended that you use your own config files though to allow better resource allocations

### Profiles
You can use multiple profiles at once, separating them with a comma. This is described in the Nextflow [documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) 

### Configs
Common configuration options are set in `conf/base.config`. Workflow specific configuration options are set in `conf/nanopore.config` and `conf/illumina.config` They are described and set to sensible defaults (as suggested in the [nCoV-2019 novel coronavirus bioinformatics protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html "nCoV-2019 novel coronavirus bioinformatics protocol"))

---------------
