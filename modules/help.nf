
def printHelp() {
  log.info """
  Usage:
    nextflow run connor-lab/ncov2019-artic-nf -profile (conda,mamba,slurm,lsf) --prefix [prefix] [workflow-options]

  Description:
    Turn Nanopore SARS-CoV2 sequencing reads generated by tiling amplification protocols (with Primal Bed file schemes)
    into consensus sequences with ARTIC (https://github.com/artic-network/fieldbioinformatics)

    All options set via CLI can be found below or looked at in the "nextflow.config" file

  Nextflow arguments (single DASH):
    -profile                  Allowed values: conda, mamba, slurm, lsf, nml

  Mandatory workflow arguments:
    --prefix                  A (unique) string prefix for output files.
                               Sequencing run name is a good choice e.g DDMMYY_MACHINEID_RUN_FLOWCELLID.
    --basecalled_fastq        The output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass 
                               This directory should contain barcodeXXX subdirectories which are auto-detected

    Optional:
      ## Basic Args ##
      --outdir                Output directory (Default: ./results)
      --cache                 Path to directory to create and reuse conda envs

      ## Scheme Related ##
      --scheme_repo_url       Repo to download your primer scheme from (Default: 'https://github.com/phac-nml/primer-schemes.git')
      --scheme_version        ARTIC scheme version (Default: 'freed')
      --scheme                Scheme name (Default: 'nCoV-2019')
      --local_scheme          Path to local scheme directory with format matching that found in the description

      ## Read Filtering ##
      --min_reads_per_barcode  Minimum number of reads accepted for a single sample when supplying demultiplexed Fastq
                                files as input. Barcodes having fewer reads are ignored. (Default: 100)
      --min_reads_guppyplex    Minimun number of reads accepted for a single sample after artic guppyplex read filtering
                                (Default: 10)
      --min_length             Minimum read length for artic guppyplex (Default: 800)
      --max_length             Maximum read length for artic guppyplex (Default: 1600)

      ## Artic Minion Args ##
      --clair3_model          Specify the string name or directory of the Clair3 model to use
                                Clair3 checks the fastq header for model info which if not found AND this parameter is not set will error
                                (Default: Null)
      --normalise             Normalise input reads for each amplicon. Set to 0 to not run normalization
                               (Default: 500)
      --no_frameshift         Add in VCF filter --no-frameshift argument that checks for %3==0 allele calls
                               (Default: false)

      ## Ncov-tools Specific ##
      --ncov                  Path to ncov-tools config file (Default: baseDir/extra_data/config.yaml)
                               Likely not one to modify but just in case

      ## Nextclade ##
      --nextclade_dataset     Name of nextclade dataset to use (Default: sars-cov-2)
      --nextclade_tag         Nextclade dataset tag to use (Default: latest)

      ## Metadata, QC, and IRIDA Uploads ##
      --skip_correct_n         Skip running custom python/bcftools N correction step (Default: False)
      --pcr_primers            Path to qPCR primers to check for mutations in
                                (Default: baseDir/extra_data/pcr_primers.bed)
      --sequencing_technology  Add in the sequencing technology used when running the pipeline
                                (Default: Nanopore)
      --irida                  Path to Irida samplesheet.tsv file to: provide metadata for final output and ncov-tools,
                                rename samples (if they are barcodeXX), and provide needed info for IRIDA Uploads
                                (Default: '')
                                samplesheet.tsv should be formated as [sample, run, barcode, project_id, ct]
      --upload_irida           Path to IRIDA instance config file to which the data will be uploaded to. Must be set
                                with a --irida samplesheet.tsv to know which project the sample goes to
                                (Default: '')

      ## Other Args ##
      --gff                   Path to annotation gff for variant consequence calling and typing. (Default: unset, don't run typing unless set)
      --yaml                  Path to YAML file with typing schemes.
                               Format: { <typing_scheme_name> : { coverage: <float>, variants: <gene_name>: <[ D614G, IHV68I ]> }}
      --max_memory            Maximum memory to allow to be allocated when running processes
      --max_cpus              Maximum cpus to allow to be allocated when running processes
      --max_time              Maximum time to allow to be allocated when running processes
  """.stripIndent()
}
