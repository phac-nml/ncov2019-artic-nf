
def printHelp() {
  log.info"""
  Usage:
    nextflow run connor-lab/ncov2019-artic-nf -profile (singularity,docker,conda,nml) ( --illumina | --nanpolish | --medaka ) --prefix [prefix] [workflow-options]

  Description:
    Turn Nanopore or Illumina SARS-CoV2 sequencing reads generated by ARTIC tiling amplification protocols into consensus sequences.
      - Nanopore: ARTIC (https://github.com/artic-network/fieldbioinformatics)
      - Illumina: iVar (https://github.com/andersen-lab/ivar)

    All options set via CLI can be set in conf directory

  Nextflow arguments (single DASH):
    -profile                  Allowed values: conda, singularity, docker, nml, [COG-UK institutional profile]

  Mandatory workflow arguments (mutually exclusive):
    --illumina                Run the Illumina workflow
    --nanopolish              Run the Nanopore/Nanopolish workflow (https://github.com/jts/nanopolish)
    --medaka                  Run the Nanopore/Medaka workflow (https://github.com/nanoporetech/medaka)

  Nanopore workflow options:
    Mandatory:
      --prefix                A (unique) string prefix for output files.
                              Sequencing run name is a good choice e.g DDMMYY_MACHINEID_RUN_FLOWCELLID.
      --basecalled_fastq      The output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. 
                              This should contain barcodeXXX directories (NML instance), which are 
                              auto-detected and analysed in parallel.
      --fast5_pass            Directory containing fast5 files - usually fast5_pass. NOT REQUIRED FOR MEDAKA WORKFLOW.
      --sequencing_summary    Path to sequencing_summary.txt. NOT REQUIRED FOR MEDAKA WORKFLOW.

    Optional:
      --outdir                Output directory (Default: ./results)
 
      --schemeVersion         ARTIC scheme version (Default: '2kb_resende')
      --schemeRepoURL         Repo to download your primer scheme from (Default: 'https://github.com/phac-nml/resende-ncov2019.git')
      --schemeDir             Directory within schemeRepoURL that contains primer schemes (Default: 'primer_schemes')
      --scheme                Scheme name (Default: 'nCoV-2019')
 
      --min_length            Minimum read length for artic guppyplex (Default: 1600)
      --max_length            Maximum read length for artic guppyplex (Default: 2400)
      --bwa                   Use BWA for mapping Nanopore reads (Default: false, uses Minimap2)
      --outCram               Output cram instead of bam files (Default: false)
      --minReadsPerBarcode    Minimum number of reads accepted for a single barcode when supplying deplexed Fastq
                              files as input. Barcodes having fewer reads are ignored. (Default: 100)
                              
      --ncov                  Path to ncov-tools config file (Default: baseDir/extra_data/config.yaml)
      --irida                 Path to Irida sample_list.tsv file to rename data and upload to Irida (Default: false)
                              Sample_list.tsv should be formated as [sample, run, barcode, project_id, ct]
      --upload_irida          Path to IRIDA instance config file to which the data will be uploaded to. Must be set
                              with a --irida sample sheet to know which project the sample goes to
 
  Illumina workflow options:
    Mandatory:
      --prefix                A (unique) string prefix for output files.
                              Sequencing run name is a good choice e.g DDMMYY_MACHINEID_RUN_FLOWCELLID.
      --directory             Path to a directory containing paired-end Illumina reads. 
                              Reads will be found and paired RECURSIVELY beneath this directory.
    Optional:
      --outdir                Output directory (Default: ./results)

      --schemeVersion         ARTIC scheme version (Default: 'V2')
      --schemeRepoURL         Repo to download your primer scheme from (Default: 'https://github.com/artic-network/artic-ncov2019')
      --schemeDir             Directory within schemeRepoURL that contains primer schemes (Default: 'primer_schemes')
      --scheme                Scheme name (Default: 'nCoV-2019')
 
      --bed                   Path to iVar-compatible bed file, also requires --ref
                              Overrides --scheme* options. (Default: unset, download scheme from git)
      --ref                   Path to iVar-compatible reference fasta file, also requires --bed 
                              Overrides --scheme* options. (Default: unset, download scheme from git)

      --allowNoprimer         Allow reads that don't have primer sequence? 
                              Depends on your library prep method: ligation == false, tagmentation == true (Default: true)
      --illuminaKeepLen       Length (bp) of reads to keep after primer trimming (Default: 20)
      --illuminaQualThreshold Sliding window quality threshold for keeping 
                              reads after primer trimming (Default: 20)
      --mpileupDepth          Mpileup depth (Default: 100000)
      --ivarFreqThreshold     iVar frequency threshold for consensus variant (ivar consensus -t, Default: 0.75)
      --ivarMinDepth          Minimum coverage depth to call variant (ivar consensus -m; ivar variants -m, Default: 10)
      --ivarMinFreqThreshold  iVar frequency threshold to call variant (ivar variants -t, Default: 0.25)
      --ivarMinVariantQuality iVar minimum mapQ to call variant (ivar variants -q, Default: 20)

      --cram                  Input is CRAM files not fastq (Default: false)
      --outCram               Output cram instead of bam files (Default: false)
  """.stripIndent()
}
