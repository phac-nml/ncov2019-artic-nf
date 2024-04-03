/*
    Main Pipeline Workflows
        - ARTIC ncov nanopore nanopolish workflow
        - ARTIC ncov nanopore medata workflow
*/
// Modules to include
include {
    articDownloadScheme ;
    articGuppyPlex ;
    articMinION ;
    articRemoveUnmappedReads
} from '../modules/artic.nf' 

include {
    makeQCCSV ;
    writeQCSummaryCSV ;
    correctQCSummaryCSV
} from '../modules/qc.nf'

include {
    renameSamples ;
    accountNoReadsInput ;
    accountReadFilterFailures ;
    generateFastqIridaReport ;
    generateFastaIridaReport ;
    generateFast5IridaReport ;
    correctFailNs ;
    runNcovTools ;
    snpDists ;
    uploadIridaFiles ;
    uploadCorrectN ;
    outputVersions
} from '../modules/nml.nf'

// Subworkflows to include
include {Genotyping} from './typing.nf'

/*
    Initialize channels from params
*/
// Nanopolish required channels, will be ignored when running medaka but still passed to the process
ch_fast5s = params.fast5_pass ? file(params.fast5_pass, type: 'dir', checkIfExists: true) : []
ch_seqsum = params.sequencing_summary ? file(params.sequencing_summary, type: 'file', checkIfExists: true) : []

// ncov-tools config file channel
ch_ncov_config = params.ncov ? file(params.ncov, type: 'file', checkIfExists: true) : []

// Optional channels
ch_local_scheme = params.local_scheme ? file(params.local_scheme, type: 'dir', checkIfExists: true) : []
ch_irida_metadata = params.irida ? file(params.irida, type: 'file', checkIfExists: true) : []
ch_irida_upload_conf = params.upload_irida ? file(params.upload_irida, type: 'file', checkIfExists: true) : []

workflow articNanopore {
    take:
    ch_fastqs
    ch_filtered_fastqs

    main:
    // Tool version tracking
    ch_versions = Channel.empty()

    // =============================== //
    // Logic Checks - Probably unneeded here
    // =============================== //
    if ( params.nanopolish ) { 
        if ( (! ch_fast5s) || (! ch_seqsum) ) {
            log.error("ERROR: Nanopolish given as input pipeline but the fast5 pass directory or sequencing_summary file weren't found")
            System.exit(1)
        }
    }

    // =============================== //
    // Scheme and Reference
    // =============================== //
    if (! ch_local_scheme) {
        articDownloadScheme()
        // Correct scheme and create files subworkflow to-do
    } else {
        articDownloadScheme() // remove later for just the check below
        // Correct scheme and create files subworkflow to-do
    }

    // =============================== //
    // Artic Size Filtering and Renaming
    // =============================== //
    articGuppyPlex(
        ch_fastqs
    )
    ch_versions = ch_versions.mix(articGuppyPlex.out.versions.first())
    ch_fastqs = articGuppyPlex.out.fastq

    // Rename if we have IRIDA metadata
    if ( ch_irida_metadata ) {
        renameSamples(
            ch_fastqs,
            ch_irida_metadata
        )

        // Re-map to create samplename path tuple again
        renameSamples.out.fastq
          .map{ fastq -> [ fastq.baseName.replaceAll(~/\.fastq.*$/, ''), fastq] }
          .set{ ch_fastqs }
    }

    // =============================== //
    // Failed Sample Tracking
    // =============================== //
    // Failing input read count check
    accountNoReadsInput(
        ch_filtered_fastqs.collect{ it[1] },
        ch_irida_metadata
    )
    accountNoReadsInput.out.count_filter
        .ifEmpty(file('placeholder_accountNoReadsInput.txt'))
        .set{ ch_noReadsTracking }
    
    // Failing size selection read count check
    ch_fastqs
        .branch{
            pass: it[1].countFastq() > params.minReadsArticGuppyPlex
            filtered: it[1].countFastq() <= params.minReadsArticGuppyPlex
        }.set{ ch_fastqs }
    accountReadFilterFailures(
        ch_fastqs.filtered
            .collect{ it[1] },
        ch_irida_metadata
    )
    accountReadFilterFailures.out.size_filter
        .ifEmpty(file('placeholder_accountReadFilterFailures.txt'))
        .set{ ch_filterReadsTracking }

    // =============================== //
    // Artic nCoV Minion Pipleine
    // =============================== //
    articMinION(
        ch_fastqs.pass,
        ch_fast5s,
        ch_seqsum,
        articDownloadScheme.out.scheme
    )
    ch_versions = ch_versions.mix(articMinION.out.versions.first())

    articRemoveUnmappedReads(
        articMinION.out.mapped
    )
    ch_versions = ch_versions.mix(articRemoveUnmappedReads.out.versions.first())

    // =============================== //
    // Failing N Position Adjustment
    // =============================== //
    ch_corrected_fasta = Channel.empty()
    if ( params.correctN ) {
        correctFailNs(
            articMinION.out.ptrim
                .join(articMinION.out.ptrimbai, by:0)
                .join(articMinION.out.consensus_fasta, by:0)
                .join(articMinION.out.fail_vcf, by:0),
            articDownloadScheme.out.reffasta
        )
        ch_versions = ch_versions.mix(correctFailNs.out.versions.first())
        ch_corrected_fasta = correctFailNs.out.corrected_consensus
    }

    // =============================== //
    // Run ncov-tools
    // =============================== //
    runNcovTools(
        ch_ncov_config, 
        articDownloadScheme.out.reffasta, 
        articDownloadScheme.out.ncov_amplicon, 
        articMinION.out.all
            .collect(),
        articDownloadScheme.out.bed,
        ch_irida_metadata,
        ch_corrected_fasta
            .collect{ it[1] }
            .ifEmpty([])
    )
    ch_versions = ch_versions.mix(runNcovTools.out.versions.first())

    // =============================== //
    // Run QC
    // =============================== //
    snpDists(runNcovTools.out.aligned)
    ch_versions = ch_versions.mix(snpDists.out.versions.first())

    // Making final CSV file with a few steps
    makeQCCSV(
        articMinION.out.ptrim
            .join(articMinION.out.consensus_fasta, by: 0)
            .join(articMinION.out.vcf, by: 0)
            .combine(articDownloadScheme.out.reffasta)
            .combine(runNcovTools.out.lineage)
            .combine(runNcovTools.out.ncovtools_qc)
            .combine(runNcovTools.out.ncovtools_negative)
            .combine(runNcovTools.out.snpeff_path)
            .combine(articDownloadScheme.out.bed),
        ch_irida_metadata,
        params.pcr_primers
    )

    // Adding pass/fail column
    makeQCCSV.out.csv.splitCsv()
        .unique()
        .branch {
            header: it[-1] == 'nextflow_qc_pass'
            fail: it[-1] == 'FALSE'
            pass: it[-1] == 'TRUE'
        }.set { qc }

    // Concat to final CSV
    writeQCSummaryCSV(
        qc.header
            .concat(qc.pass)
            .concat(qc.fail)
            .toList()
    )

    // Fix final CSV
    //  Null values, neg controls, read tracking, etc.
    correctQCSummaryCSV(
        writeQCSummaryCSV.out,
        ch_noReadsTracking,
        ch_filterReadsTracking
    )

    // =============================== //
    // Uploads and Final Tracking
    // =============================== //
    // IRIDA Upload Samplesheets made even if not uploading (for tracking)
    if ( ch_irida_metadata ) {
        // Adding size filter for IRIDA uploads, 1kB needed or it breaks
        generateFastqIridaReport(
            ch_fastqs.pass
                .concat(ch_fastqs.filtered)
                .filter{ it[1].size() > 1024 }
                .collect{ it[1] },
            ch_irida_metadata
        )
        generateFastaIridaReport(
            articMinION.out.consensus_fasta
              .collect{ it[1] },
            ch_irida_metadata
        )

        // We upload now if given a config
        if ( ch_irida_upload_conf ) {
            // Fast5s made here as its slow and not needed for normal run tracking
            ch_fast5_upload = Channel.empty()
            if ( params.fast5_pass ) {
                generateFast5IridaReport(
                    ch_fast5Pass,
                    ch_irida_metadata
                )
                ch_fast5_upload = generateFast5IridaReport.out.fast5_dir
            }

            uploadIridaNanopolish(
                generateFastqIridaReport.out.fastq_dir,
                generateFastaIridaReport.out.fasta_dir,
                ch_fast5_upload
                    .ifEmpty([]),
                ch_irida_upload_conf,
                correctQCSummaryCSV.out
            )
            ch_versions = ch_versions.mix(uploadIridaNanopolish.out.versions.first())

            // Upload corrected fasta files if any
            if ( params.correctN ) {
                uploadCorrectN(
                    ch_corrected_fasta
                      .collect{ it[1] },
                    ch_irida_upload_conf,
                    ch_irida_metadata
                )
            }
        }
    }

    outputVersions(ch_versions.collect())

    // =============================== //
    // Emit for Genotyping Subworkflow
    // =============================== //
    emit:
        reffasta = articDownloadScheme.out.reffasta
        vcf = articMinION.out.vcf
}

/*
  Control Workflow
    Controls what pipeline to utilize to get results
    Could likely be removed as the workflow has been generalized
*/
workflow articNcovNanopore {
    take:
    ch_fastq_pass
    ch_fastq_filtered

    main:
    articNanopore(
        ch_fastq_pass,
        ch_fastq_filtered
    )
    // For typing workflow
    ch_nanopore_vcf = articNanopore.out.vcf
    ch_nanopore_reffasta = articNanopore.out.reffasta

    // Typing workflow that isn't really used but can be kept in
    if ( params.gff ) {
        Channel.fromPath("${params.gff}")
            .set{ ch_refGff }

        Channel.fromPath("${params.yaml}")
            .set{ ch_typingYaml }

        Genotyping(
            ch_nanopore_vcf,
            ch_refGff,
            ch_nanopore_reffasta,
            ch_typingYaml
        )
    }
}
