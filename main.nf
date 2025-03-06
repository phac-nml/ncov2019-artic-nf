#!/usr/bin/env nextflow

// enable dsl2
// ===============================
nextflow.enable.dsl = 2

// include modules
// ===============================
include { printHelp } from './modules/help.nf'

// import main workflow
// ===============================
include { articNcovNanopore } from './workflows/articNcovNanopore.nf' 

// Help
// ===============================
if ( params.help ) {
    printHelp()
    exit 0
}

// Previous param warnings/exit
def previousParams = [
    'medakaModel',
    'schemeRepoURL',
    'schemeVersion',
    'minReadsPerBarcode',
    'minReadsGuppyPlex',
    'correctN',
    'sequencingTechnology',
    'csqAfThreshold',
    'csqDpThreshold',
    'nanopolish',
    'medaka',
    'medaka_model',
    'fast5_pass',
    'sequencing_summary',
    'bwa',
    'no_longshot'
]

def overlap = previousParams.intersect(params.keySet())
if ( overlap != [] ) {
    log.error """The following previously used params were given: "$overlap"
    Please resubmit with the new parameters (if there is an equivalent) that can be found when running the --help command
    """
    System.exit(1)
}

// Pipeline required input checks
// ===============================
if (! params.basecalled_fastq ) {
    log.error """Please supply a directory containing basecalled fastqs with `--basecalled_fastq`
    This can either be:
        - The output directory from guppy or dorado, usually fastq_pass, with directories called barcodeXX containing fastq-formatted
        - A flat directory of named ".fastq" or ".fastq.gz" files
    """.stripIndent()
    System.exit(1)
}

// Prefix existance and formatting check
if ( ! params.prefix ) {
    log.error """Please supply a run prefix for your output files with "--prefix <prefix>" """
    System.exit(1)
} else if ( params.prefix =~ /\// ) {
    log.error "The --prefix that you supplied contains a \"/\", please replace it with another character"
    System.exit(1)
} 

// Main Workflow
// ===============================
workflow {
    // ===============================
    // Input barcodes and fastqs check
    // ===============================
    nanoporeBarcodeDirs = file("${params.basecalled_fastq}/barcode*", type: 'dir', maxdepth: 1 )
    nanoporeFastqs = file("${params.basecalled_fastq}/*.fastq*", type: 'file', maxdepth: 1)

    // If barcodes, check number of reads and then go to pass/filtered branches
    if ( nanoporeBarcodeDirs ) {
        Channel.fromPath( nanoporeBarcodeDirs )
            .filter( ~/.*barcode[0-9]{1,4}$/ )
            .map{ dir -> 
                def count = 0
                for ( x in dir.listFiles() ) {
                    if ( x.isFile() ) {
                        count += x.countFastq()
                    }
                }
                return [ dir.baseName, file(dir), count ]
            }.set{ ch_initial_fastq_dirs }

    // No barcodes, check number of reads in single file and then go to pass/filtered branches
    } else if ( nanoporeFastqs ) {
        Channel.fromPath( nanoporeFastqs )
            .map{ in_f -> 
                def count = 0
                if (in_f.isFile()) {
                    count = in_f.countFastq()
                }
                return [ in_f.baseName.replaceAll(~/\.fastq.*$/, ''), file(in_f), count ]
            }.set{ ch_initial_fastq_dirs }

    } else {
        log.error"""Couldn't detect whether your Nanopore data was in barcoded directories or flat directories
        Please use "--basecalled_fastq <DIR>" to point to either a barcoded nanopore directory or a directory of flat fastq files
        See the README for examples of each if needed
        """
        System.exit(1)
    }

    // Final input pathes and filtering
    ch_initial_fastq_dirs
        .branch{ it ->
            pass: it[2] > params.min_reads_per_barcode
                return [ it[0], it[1] ]
            filtered: it[2] <= params.min_reads_per_barcode
                return [ it[0], it[1] ]
        }
        .set{ ch_fastqs }

    // Execute Main Named process
    main:
    articNcovNanopore(
        ch_fastqs.pass,
        ch_fastqs.filtered
    )
}
// :)
