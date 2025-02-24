/*
    Main Pipeline Workflows
        - ARTIC ncov nanopore clair3 workflow
*/
// Modules to include
include {
    checkFastqForModel ;
    articDownloadModels ;
    articGuppyPlex ;
    articMinION
} from '../modules/artic.nf' 
include {
    makeQCCSV ;
    writeQCSummaryCSV ;
    correctQCSummaryCSV
} from '../modules/qc.nf'
include {
    renameBarcodeSamples ;
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
include {
    nextcladeDatasetGet ;
    nextcladeRun
} from '../modules/nextclade.nf'

// Subworkflows to include
include { schemeValidate } from './schemeValidate.nf'
include { Genotyping } from './typing.nf'

/*
    Initialize channels from params
*/
// ncov-tools config file channel
ch_ncov_config = params.ncov ? file(params.ncov, type: 'file', checkIfExists: true) : []

// Optional channels
ch_irida_metadata = params.irida ? file(params.irida, type: 'file', checkIfExists: true) : []
ch_irida_upload_conf = params.upload_irida ? file(params.upload_irida, type: 'file', checkIfExists: true) : []
ch_pcr_primers = params.pcr_primers ? file(params.pcr_primers, type: 'file', checkIfExists: true) : []

/*
  Main Workflow
*/
workflow articNcovNanopore {
    take:
    ch_fastqs
    ch_filtered_fastqs

    main:
    // Tool version tracking
    ch_versions = Channel.empty()

    // =============================== //
    // Scheme and Reference
    // =============================== //
    schemeValidate()
    ch_scheme = schemeValidate.out.scheme               // channel: [ val(scheme_version), path(scheme) ]
    ch_reference = schemeValidate.out.reference         // channel: [ path(reference.fasta) ]
    ch_primer_bed = schemeValidate.out.primer_bed       // channel: [ path(primer_bed) ]
    ch_amplicon = schemeValidate.out.amplicon_bed       // channel: [ path(amplicon_bed) ]
    ch_primer_prefix = schemeValidate.out.primer_prefix // channel: [ val(primer_prefix) ]

    // Version tracking update from subworkflow
    ch_versions = ch_versions.mix(schemeValidate.out.versions)

    // =============================== //
    // Artic Size Filtering and Renaming
    // =============================== //
    // Check if the clair3 model has been selected, or can be selected from fastqs
    checkFastqForModel( 
        ch_fastqs
    )
    ch_fastqs_checked = check_model.out.check_done.join(ch_fastqs)

    articGuppyPlex(
        ch_fastqs_checked
    )
    ch_versions = ch_versions.mix(articGuppyPlex.out.versions)
    ch_fastqs = articGuppyPlex.out.fastq

    // Rename if we have IRIDA metadata
    if ( ch_irida_metadata ) {
        // Renames samples that are from barcode directories
        renameBarcodeSamples(
            ch_fastqs,
            ch_irida_metadata
        )
        ch_versions = ch_versions.mix(renameBarcodeSamples.out.versions)

        // Re-map to create samplename path tuple again
        renameBarcodeSamples.out.fastq
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
        .ifEmpty([])
        .set{ ch_no_reads_tracking }
    
    // Failing size selection read count check
    //  Remap fastqs channel to branch the pass/filtered fastq files and use those accordingly
    ch_fastqs
        .branch{
            pass: it[1].countFastq() > params.min_reads_guppyplex
            filtered: it[1].countFastq() <= params.min_reads_guppyplex
        }.set{ ch_fastqs }
    accountReadFilterFailures(
        ch_fastqs.filtered
            .collect{ it[1] },
        ch_irida_metadata
    )
    accountReadFilterFailures.out.size_filter
        .ifEmpty([])
        .set{ ch_filtered_reads_tracking }

    // =============================== //
    // Artic nCoV Minion Pipleine
    // =============================== //
    // Fetch the R10 models
    articDownloadModels() 
    
    articMinION(
        ch_fastqs.pass,
        ch_reference,
        ch_primer_bed,
    )
    ch_versions = ch_versions.mix(articMinION.out.versions)

    // =============================== //
    // Failing N Position Adjustment
    // =============================== //
    ch_corrected_fasta = Channel.empty()
    if ( ! params.skip_correct_n ) {
        correctFailNs(
            articMinION.out.ptrim
                .join(articMinION.out.ptrimbai, by:0)
                .join(articMinION.out.consensus_fasta, by:0)
                .join(articMinION.out.fail_vcf, by:0),
            ch_reference
        )
        ch_versions = ch_versions.mix(correctFailNs.out.versions)
        ch_corrected_fasta = correctFailNs.out.corrected_consensus
    }

    // =============================== //
    // Run ncov-tools
    // =============================== //
    runNcovTools(
        ch_ncov_config, 
        ch_reference, 
        ch_amplicon, 
        articMinION.out.all
            .collect(),
        ch_primer_bed,
        ch_irida_metadata,
        ch_corrected_fasta
            .collect{ it[1] }
            .ifEmpty([]),
        ch_primer_prefix
    )
    ch_versions = ch_versions.mix(runNcovTools.out.versions)

    // =============================== //
    // Run nextclade
    // =============================== //
    nextcladeDatasetGet(
        params.nextclade_dataset,
        params.nextclade_tag
    )
    ch_versions = ch_versions.mix(nextcladeDatasetGet.out.versions)

    nextcladeRun(
        articMinION.out.consensus_fasta,
        nextcladeDatasetGet.out.dataset
    )
    ch_versions = ch_versions.mix(nextcladeRun.out.versions)

    // =============================== //
    // Run QC
    // =============================== //
    snpDists(
        runNcovTools.out.aligned
    )
    ch_versions = ch_versions.mix(snpDists.out.versions)

    // Making final CSV file with a few steps
    makeQCCSV(
        articMinION.out.ptrim
            .join(articMinION.out.consensus_fasta, by: 0)
            .join(articMinION.out.vcf, by: 0)
            .join(nextcladeRun.out.tsv, by: 0),
        ch_reference,
        runNcovTools.out.lineage,
        runNcovTools.out.ncovtools_qc,
        runNcovTools.out.ncovtools_negative,
        runNcovTools.out.snpeff_path,
        ch_primer_bed,
        ch_irida_metadata,
        ch_pcr_primers,
        params.sequencing_technology
    )
    ch_versions = ch_versions.mix(makeQCCSV.out.versions)

    // Adding pass/fail column
    //  Not used anymore really but it is still here for now
    makeQCCSV.out.csv
        .splitCsv()
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
        ch_no_reads_tracking,
        ch_filtered_reads_tracking
    )
    ch_versions = ch_versions.mix(correctQCSummaryCSV.out.versions)

    // =============================== //
    // Uploads and Final Tracking
    // =============================== //
    // IRIDA Upload Samplesheets made even if not uploading (for tracking) or after run-check upload
    if ( ch_irida_metadata ) {
        // Size filter for IRIDA fastq uploads, 1kB needed or it breaks the upload
        generateFastqIridaReport(
            ch_fastqs.pass
                .concat(ch_fastqs.filtered)
                .filter{ it[1].size() > 1024 }
                .collect{ it[1] },
            ch_irida_metadata
        )
        ch_versions = ch_versions.mix(generateFastqIridaReport.out.versions)

        generateFastaIridaReport(
            articMinION.out.consensus_fasta
              .collect{ it[1] },
            ch_irida_metadata
        )
        ch_versions = ch_versions.mix(generateFastaIridaReport.out.versions)

        // Upload now if given a config
        if ( ch_irida_upload_conf ) {
            uploadIridaFiles(
                generateFastqIridaReport.out.fastq_dir,
                generateFastaIridaReport.out.fasta_dir,
                ch_irida_upload_conf,
                correctQCSummaryCSV.out.final_csv
            )
            ch_versions = ch_versions.mix(uploadIridaFiles.out.versions)

            // Upload corrected fasta files if any
            if ( ! params.skip_correct_n ) {
                uploadCorrectN(
                    ch_corrected_fasta
                      .collect{ it[1] },
                    ch_irida_upload_conf,
                    ch_irida_metadata
                )
                ch_versions = ch_versions.mix(uploadCorrectN.out.versions)
            }
        }
    }

    // =============================== //
    // Tool versions for tracking
    // =============================== //
    // Using nf-core version output
    outputVersions(
        ch_versions
            .unique()
            .collectFile(name: 'collated_versions.yml')
    )

    // =============================== //
    // Typing workflow that isn't really used but can be kept in
    // =============================== //
    if ( params.gff ) {
        Channel.fromPath("${params.gff}")
            .set{ ch_ref_gff }

        Channel.fromPath("${params.yaml}")
            .set{ ch_typing_yaml }

        Genotyping(
            articMinION.out.vcf,
            ch_ref_gGff,
            ch_reference,
            ch_typing_yaml
        )
    }

    // =============================== //
    // Completion
    // =============================== //
    workflow.onComplete {
        log.info "Workflow complete"
    }
}
