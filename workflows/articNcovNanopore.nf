// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' 
include {articGuppyPlex} from '../modules/artic.nf' 
include {articGuppyPlexFlat} from '../modules/artic.nf' 
include {articMinIONNanopolish} from  '../modules/artic.nf' 
include {articMinIONMedaka} from  '../modules/artic.nf'
include {articRemoveUnmappedReads} from '../modules/artic.nf' 

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'
include {correctQCSummaryCSV} from '../modules/qc.nf'

include {bamToCram} from '../modules/out.nf'

include {collateSamples} from '../modules/upload.nf'

include {accountNoReadsInput} from '../modules/nml.nf'
include {renameSamples} from '../modules/nml.nf'
include {accountReadFilterFailures} from '../modules/nml.nf'
include {generateFastqIridaReport} from '../modules/nml.nf'
include {generateFastaIridaReport} from '../modules/nml.nf'
include {generateFast5IridaReport} from '../modules/nml.nf'
include {correctFailNs} from '../modules/nml.nf'
include {runNcovTools} from '../modules/nml.nf'
include {snpDists} from '../modules/nml.nf'
include {uploadIridaNanopolish} from '../modules/nml.nf'
include {uploadIridaMedaka} from '../modules/nml.nf'
include {uploadCorrectN} from '../modules/nml.nf'
include {outputVersions} from '../modules/nml.nf'


// import subworkflows
include {CLIMBrsync} from './upload.nf'
include {Genotyping} from './typing.nf'

// workflow component for artic pipeline
workflow sequenceAnalysisNanopolish {
    take:
      ch_runFastqDirs
      ch_badFastqDirs
      ch_fast5Pass
      ch_seqSummary
      ch_irida
    
    main:

      ch_versions = Channel.empty()


      articDownloadScheme()

      accountNoReadsInput(ch_badFastqDirs.collect(),
                          ch_irida)
      accountNoReadsInput.out.count_filter
        .ifEmpty(file('placeholder_accountNoReadsInput.txt'))
        .set{ ch_noReadsTracking }
      
      articGuppyPlex(ch_runFastqDirs.flatten())
      ch_versions = ch_versions.mix(articGuppyPlex.out.versions.first())

      if (params.irida) {

       renameSamples(articGuppyPlex.out.fastq
                                       .combine(ch_irida))

       accountReadFilterFailures(renameSamples.out.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect(),
                                  ch_irida)
       accountReadFilterFailures.out.size_filter
        .ifEmpty(file('placeholder_accountReadFilterFailures.txt'))
        .set{ ch_filterReadsTracking }

       articMinIONNanopolish(renameSamples.out
                                          .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                          .combine(articDownloadScheme.out.scheme)
                                          .combine(ch_fast5Pass)
                                          .combine(ch_seqSummary))
       ch_versions = ch_versions.mix(articMinIONNanopolish.out.versions.first())

       generateFastqIridaReport(renameSamples.out.collect(),
                                ch_irida)

       generateFastaIridaReport(articMinIONNanopolish.out.consensus_fasta.collect(),
                                ch_irida)
      }
      else {

       accountReadFilterFailures(articGuppyPlex.out.fastq.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect(),
                                  ch_irida)
       accountReadFilterFailures.out.size_filter
        .ifEmpty(file('placeholder_accountReadFilterFailures.txt'))
        .set{ ch_filterReadsTracking }

       articMinIONNanopolish(articGuppyPlex.out.fastq
                                          .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                          .combine(articDownloadScheme.out.scheme)
                                          .combine(ch_fast5Pass)
                                          .combine(ch_seqSummary))
       ch_versions = ch_versions.mix(articMinIONNanopolish.out.versions.first())
      }

      articRemoveUnmappedReads(articMinIONNanopolish.out.mapped)
      ch_versions = ch_versions.mix(articRemoveUnmappedReads.out.versions.first())

      if (params.correctN) {
        correctFailNs(articMinIONNanopolish.out.ptrim
                          .join(articMinIONNanopolish.out.ptrimbai, by:0)
                          .join(articMinIONNanopolish.out.consensus_fasta, by:0)
                          .join(articMinIONNanopolish.out.fail_vcf, by:0),
                          articDownloadScheme.out.reffasta)
        ch_versions = ch_versions.mix(correctFailNs.out.versions.first())

        correctFailNs.out.corrected_consensus.collect()
              .ifEmpty(file('placeholder.txt'))
              .set{ ch_corrected }
      }

      Channel.fromPath("${params.ncov}")
             .set{ ch_ncov }

      runNcovTools(ch_ncov, 
                      articDownloadScheme.out.reffasta, 
                      articDownloadScheme.out.ncov_amplicon, 
                      articMinIONNanopolish.out[0].collect(),
                      articDownloadScheme.out.bed,
                      ch_irida,
                      ch_corrected)
      ch_versions = ch_versions.mix(runNcovTools.out.versions.first())
      
      snpDists(runNcovTools.out.aligned)
      ch_versions = ch_versions.mix(snpDists.out.versions.first())

      makeQCCSV(articMinIONNanopolish.out.ptrim
                                     .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                                     .join(articMinIONNanopolish.out.vcf, by: 0)
                                     .combine(articDownloadScheme.out.reffasta)
                                     .combine(runNcovTools.out.lineage)
                                     .combine(runNcovTools.out.ncovtools_qc)
                                     .combine(runNcovTools.out.ncovtools_negative)
                                     .combine(ch_irida)
                                     .combine(runNcovTools.out.snpeff_path)
                                     .combine(articDownloadScheme.out.bed),
                params.pcr_primers)

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'nextflow_qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     correctQCSummaryCSV(writeQCSummaryCSV.out,
                        ch_noReadsTracking,
                        ch_filterReadsTracking)

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out.mapped_bam))

     if (params.outCram) {
        bamToCram(articMinIONNanopolish.out.ptrim.map{ it[0] } 
                        .join (articDownloadScheme.out.reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }

     if (params.irida) {
       if (params.upload_irida) {
         Channel.fromPath("${params.upload_irida}")
             .set{ ch_upload }

         generateFast5IridaReport(ch_fast5Pass, ch_irida)

         uploadIridaNanopolish(generateFastqIridaReport.out, generateFastaIridaReport.out, generateFast5IridaReport.out, ch_upload, correctQCSummaryCSV.out)
         ch_versions = ch_versions.mix(uploadIridaNanopolish.out.versions.first())

         if (params.correctN) {
           uploadCorrectN(correctFailNs.out.corrected_consensus.collect(),
                          ch_upload,
                          ch_irida)
         }
       }
     }

     outputVersions(ch_versions.collect())

    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONNanopolish.out.vcf
}

workflow sequenceAnalysisMedaka {
    take:
      ch_runFastqDirs
      ch_badFastqDirs
      ch_irida

    main:
    
      ch_versions = Channel.empty()

      articDownloadScheme()

      accountNoReadsInput(ch_badFastqDirs.collect(),
                          ch_irida)
      accountNoReadsInput.out.count_filter
        .ifEmpty(file('placeholder_accountNoReadsInput.txt'))
        .set{ ch_noReadsTracking }

      articGuppyPlex(ch_runFastqDirs.flatten())
      ch_versions = ch_versions.mix(articGuppyPlex.out.versions.first())

      if (params.irida) {

       renameSamples(articGuppyPlex.out.fastq
                                       .combine(ch_irida))

       accountReadFilterFailures(renameSamples.out.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect(),
                                  ch_irida)
       accountReadFilterFailures.out.size_filter
        .ifEmpty(file('placeholder_accountReadFilterFailures.txt'))
        .set{ ch_filterReadsTracking }

       articMinIONMedaka(renameSamples.out
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(articDownloadScheme.out.scheme))
       ch_versions = ch_versions.mix(articMinIONMedaka.out.versions.first())

       generateFastqIridaReport(renameSamples.out.collect(),
                                ch_irida)

       generateFastaIridaReport(articMinIONMedaka.out.consensus_fasta.collect(),
                                ch_irida)
      } else {

       accountReadFilterFailures(articGuppyPlex.out.fastq.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect(),
                                  ch_irida)
       accountReadFilterFailures.out.size_filter
        .ifEmpty(file('placeholder_accountReadFilterFailures.txt'))
        .set{ ch_filterReadsTracking }

       articMinIONMedaka(articGuppyPlex.out.fastq
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(articDownloadScheme.out.scheme))
       ch_versions = ch_versions.mix(articMinIONMedaka.out.versions.first())
      }

      articRemoveUnmappedReads(articMinIONMedaka.out.mapped)
      ch_versions = ch_versions.mix(articRemoveUnmappedReads.out.versions.first())

      if (params.correctN) {
        correctFailNs(articMinIONMedaka.out.ptrim
                          .join(articMinIONMedaka.out.ptrimbai, by:0)
                          .join(articMinIONMedaka.out.consensus_fasta, by:0)
                          .join(articMinIONMedaka.out.fail_vcf, by:0),
                          articDownloadScheme.out.reffasta)
        ch_versions = ch_versions.mix(correctFailNs.out.versions.first())

        correctFailNs.out.corrected_consensus.collect()
              .ifEmpty(file('placeholder.txt'))
              .set{ ch_corrected }
      }

      Channel.fromPath("${params.ncov}")
             .set{ ch_ncov }

      runNcovTools(ch_ncov, 
                      articDownloadScheme.out.reffasta, 
                      articDownloadScheme.out.ncov_amplicon, 
                      articMinIONMedaka.out[0].collect(),
                      articDownloadScheme.out.bed,
                      ch_irida,
                      ch_corrected)
      ch_versions = ch_versions.mix(runNcovTools.out.versions.first())
      
      snpDists(runNcovTools.out.aligned)
      ch_versions = ch_versions.mix(snpDists.out.versions.first())

      makeQCCSV(articMinIONMedaka.out.ptrim
                                     .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                                     .join(articMinIONMedaka.out.vcf, by: 0)
                                     .combine(articDownloadScheme.out.reffasta)
                                     .combine(runNcovTools.out.lineage)
                                     .combine(runNcovTools.out.ncovtools_qc)
                                     .combine(runNcovTools.out.ncovtools_negative)
                                     .combine(ch_irida)
                                     .combine(runNcovTools.out.snpeff_path)
                                     .combine(articDownloadScheme.out.bed),
                params.pcr_primers)

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'nextflow_qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     correctQCSummaryCSV(writeQCSummaryCSV.out,
                        ch_noReadsTracking,
                        ch_filterReadsTracking)

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out.mapped_bam))

     if (params.outCram) {
        bamToCram(articMinIONMedaka.out.ptrim.map{ it[0] } 
                        .join (articDownloadScheme.out.reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }

     if (params.irida) {
       if (params.upload_irida) {
         Channel.fromPath("${params.upload_irida}")
             .set{ ch_upload }

         uploadIridaMedaka(generateFastqIridaReport.out, generateFastaIridaReport.out, ch_upload, correctQCSummaryCSV.out)
         ch_versions = ch_versions.mix(uploadIridaMedaka.out.versions.first())

         if (params.correctN) {
           uploadCorrectN(correctFailNs.out.corrected_consensus.collect(),
                          ch_upload,
                          ch_irida)
         }
       }
     }

     outputVersions(ch_versions.collect())

    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONMedaka.out.vcf
}

// Write new process for analyzing flat fastq/fastq.gz files
workflow sequenceAnalysisMedakaFlat {
    take:
      ch_fastqs
      ch_irida

    main:
    
      ch_versions = Channel.empty()
      // Won't check this but still need it to get the final value
      Channel.fromPath('placeholder_accountNoReadsInput.txt')
             .set{ ch_noReadsTracking }

      articDownloadScheme()

      articGuppyPlexFlat(ch_fastqs)
      ch_versions = ch_versions.mix(articGuppyPlexFlat.out.versions.first())

      accountReadFilterFailures(articGuppyPlexFlat.out.fastq.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect(),
                                  ch_irida)
      accountReadFilterFailures.out.size_filter
        .ifEmpty(file('placeholder_accountReadFilterFailures.txt'))
        .set{ ch_filterReadsTracking }

      articMinIONMedaka(articGuppyPlexFlat.out.fastq
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(articDownloadScheme.out.scheme))
      ch_versions = ch_versions.mix(articMinIONMedaka.out.versions.first())

      articRemoveUnmappedReads(articMinIONMedaka.out.mapped)
      ch_versions = ch_versions.mix(articRemoveUnmappedReads.out.versions.first())

      if (params.correctN) {
        correctFailNs(articMinIONMedaka.out.ptrim
                          .join(articMinIONMedaka.out.ptrimbai, by:0)
                          .join(articMinIONMedaka.out.consensus_fasta, by:0)
                          .join(articMinIONMedaka.out.fail_vcf, by:0),
                          articDownloadScheme.out.reffasta)
        ch_versions = ch_versions.mix(correctFailNs.out.versions.first())

        correctFailNs.out.corrected_consensus.collect()
              .ifEmpty(file('placeholder.txt'))
              .set{ ch_corrected }
      }

      // Add ncov-tools config file
      Channel.fromPath("${params.ncov}")
             .set{ ch_ncov }

      runNcovTools(ch_ncov, 
                      articDownloadScheme.out.reffasta, 
                      articDownloadScheme.out.ncov_amplicon, 
                      articMinIONMedaka.out[0].collect(),
                      articDownloadScheme.out.bed,
                      ch_irida,
                      ch_corrected)
      ch_versions = ch_versions.mix(runNcovTools.out.versions.first())
      
      snpDists(runNcovTools.out.aligned)
      ch_versions = ch_versions.mix(snpDists.out.versions.first())

      makeQCCSV(articMinIONMedaka.out.ptrim
                                     .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                                     .join(articMinIONMedaka.out.vcf, by: 0)
                                     .combine(articDownloadScheme.out.reffasta)
                                     .combine(runNcovTools.out.lineage)
                                     .combine(runNcovTools.out.ncovtools_qc)
                                     .combine(runNcovTools.out.ncovtools_negative)
                                     .combine(ch_irida)
                                     .combine(runNcovTools.out.snpeff_path)
                                     .combine(articDownloadScheme.out.bed),
                params.pcr_primers)

      makeQCCSV.out.csv.splitCsv()
                       .unique()
                       .branch {
                           header: it[-1] == 'nextflow_qc_pass'
                           fail: it[-1] == 'FALSE'
                           pass: it[-1] == 'TRUE'
                       }
                       .set { qc }

     writeQCSummaryCSV(qc.header.concat(qc.pass).concat(qc.fail).toList())

     correctQCSummaryCSV(writeQCSummaryCSV.out,
                        ch_noReadsTracking,
                        ch_filterReadsTracking)

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out.mapped_bam))

     // Upload Data
     if (params.irida) {
       if (params.upload_irida) {
         Channel.fromPath("${params.upload_irida}")
             .set{ ch_upload }
        
         // Generate upload datasets
         generateFastqIridaReport(articGuppyPlexFlat.out.fastq.collect(), ch_irida)
         generateFastaIridaReport(articMinIONMedaka.out.consensus_fasta.collect(), ch_irida)

         // Upload
         uploadIridaMedaka(generateFastqIridaReport.out, generateFastaIridaReport.out, ch_upload, correctQCSummaryCSV.out)
         ch_versions = ch_versions.mix(uploadIridaMedaka.out.versions.first())

         if (params.correctN) {
           uploadCorrectN(correctFailNs.out.corrected_consensus.collect(),
                          ch_upload,
                          ch_irida)
         }
       }
     }

     outputVersions(ch_versions.collect())

    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONMedaka.out.vcf
}

// Process that controls what pipeline to utilize to get results
workflow articNcovNanopore {
    take:
      ch_fastqDirs
      ch_badFastqDirs
    
    main:
      // Add metadata (if there is any, otherwise this will just become false which will work when passed to different processes)
      Channel.fromPath("${params.irida}")
             .set{ ch_irida }

      // Actually run the different pipeline processes
      if ( params.nanopolish ) {
          Channel.fromPath( "${params.fast5_pass}" )
                 .set{ ch_fast5Pass }

          Channel.fromPath( "${params.sequencing_summary}" )
                 .set{ ch_seqSummary }

          sequenceAnalysisNanopolish(ch_fastqDirs, ch_badFastqDirs, ch_fast5Pass, ch_seqSummary, ch_irida)

          sequenceAnalysisNanopolish.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisNanopolish.out.reffasta.set{ ch_nanopore_reffasta }

      } else if ( params.flat ) {
          Channel.fromPath( "${params.basecalled_fastq}/*.fastq*", type: 'file', maxDepth: 1 )
                 .set{ ch_fastqs }

          sequenceAnalysisMedakaFlat(ch_fastqs, ch_irida)

          sequenceAnalysisMedakaFlat.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisMedakaFlat.out.reffasta.set{ ch_nanopore_reffasta }
      
      } else if ( params.medaka ) {
          sequenceAnalysisMedaka(ch_fastqDirs, ch_badFastqDirs, ch_irida)

          sequenceAnalysisMedaka.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisMedaka.out.reffasta.set{ ch_nanopore_reffasta }
      }

      if ( params.gff ) {
          Channel.fromPath("${params.gff}")
                 .set{ ch_refGff }

          Channel.fromPath("${params.yaml}")
                 .set{ ch_typingYaml }

          Genotyping(ch_nanopore_vcf, ch_refGff, ch_nanopore_reffasta, ch_typingYaml)

      }
      if ( params.upload ) {

        Channel.fromPath("${params.CLIMBkey}")
               .set{ ch_CLIMBkey }

        CLIMBrsync(sequenceAnalysis.out.qc_pass, ch_CLIMBkey )
      }
}
