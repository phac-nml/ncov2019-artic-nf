#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include modules
include {printHelp} from './modules/help.nf'

// import subworkflows
include {articNcovNanopore} from './workflows/articNcovNanopore.nf' 

// Help
if ( params.help ) {
    printHelp()
    exit 0
}

// Simple profile warning
if ( params.profile ) {
    log.error("ERROR: Profile should have a single dash: -profile")
    System.exit(1)
}

// Pipeline required input checks
if ( params.nanopolish ) {
    if (! params.basecalled_fastq ) {
        log.error("ERROR: Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
        System.exit(1)
    } else if (! params.fast5_pass ) {
        log.error("ERROR: Please supply a directory containing fast5 files with --fast5_pass (this is the fast5_pass directory)")
        System.exit(1)
    } else if (! params.sequencing_summary ) {
        log.error("ERROR: Please supply the path to the sequencing_summary.txt file from your run with --sequencing_summary")
        System.exit(1)
    }
} else if ( params.medaka ) {
    if (! params.basecalled_fastq ) {
        log.error("ERROR: Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
        System.exit(1)
    }
} else {
    log.error("ERROR: Please select a workflow with --nanopolish or --medaka, or use --help to print help")
    System.exit(1)
}

if ( ! params.prefix ) {
    log.error("ERROR: Please supply a prefix for your output files with --prefix")
    log.error("ERROR: For more information use --help to print help statement")
    System.exit(1)
} else {
    if ( params.prefix =~ /\// ){
        log.error("The --prefix that you supplied contains a \"/\", please replace it with another character")
        System.exit(1)
    }
} 

// Main Workflow
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
        log.error("ERROR: Couldn't detect whether your Nanopore run was barcoded or not. Use --basecalled_fastq to point to the unmodified guppy output directory.")
        System.exit(1)
    }

    // Final input pathes and filtering
    ch_initial_fastq_dirs
        .branch{ it ->
            pass: it[2] > params.minReadsPerBarcode
                return [ it[0], it[1] ]
            filtered: it[2] <= params.minReadsPerBarcode
                return [ it[0], it[1] ]
        }
        .set{ ch_fastqs }

    // Execute Main Named process
    main:
    if ( params.nanopolish || params.medaka ) {
        articNcovNanopore(ch_fastqs.pass, ch_fastqs.filtered)
    } else {
        println("Please select a workflow with --nanopolish or --medaka")
    }
}
