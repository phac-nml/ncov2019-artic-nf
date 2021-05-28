// ARTIC ncov workflow

// enable dsl2
nextflow.preview.dsl = 2

// import modules
include {articDownloadScheme} from '../modules/artic.nf' 
include {articGuppyPlex} from '../modules/artic.nf' 
include {articMinIONNanopolish} from  '../modules/artic.nf' 
include {articMinIONMedaka} from  '../modules/artic.nf'
include {articRemoveUnmappedReads} from '../modules/artic.nf' 

include {makeQCCSV} from '../modules/qc.nf'
include {writeQCSummaryCSV} from '../modules/qc.nf'
include {correctQCSummaryCSV} from '../modules/qc.nf'

include {bamToCram} from '../modules/out.nf'

include {collateSamples} from '../modules/upload.nf'

include {renameSamples} from '../modules/nml.nf'
include {accountReadFilterFailures} from '../modules/nml.nf'
include {generateFastqIridaReport} from '../modules/nml.nf'
include {generateFastaIridaReport} from '../modules/nml.nf'
include {generateFast5IridaReport} from '../modules/nml.nf'
include {correctFailNs} from '../modules/nml.nf'
include {runNcovTools} from '../modules/nml.nf'
include {snpDists} from '../modules/nml.nf'
include {uploadIrida} from '../modules/nml.nf'
include {uploadCorrectN} from '../modules/nml.nf'


// import subworkflows
include {CLIMBrsync} from './upload.nf'
include {Genotyping} from './typing.nf'

// workflow component for artic pipeline
workflow sequenceAnalysisNanopolish {
    take:
      ch_runFastqDirs
      ch_fast5Pass
      ch_seqSummary
    
    main:
      articDownloadScheme()
      
      articGuppyPlex(ch_runFastqDirs.flatten())

      if (params.irida) {
       Channel.fromPath("${params.irida}")
              .set{ ch_irida }

       renameSamples(articGuppyPlex.out.fastq
                                       .combine(ch_irida))

       accountReadFilterFailures(renameSamples.out.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect())

       articMinIONNanopolish(renameSamples.out
                                          .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                          .combine(articDownloadScheme.out.scheme)
                                          .combine(ch_fast5Pass)
                                          .combine(ch_seqSummary))

       generateFastqIridaReport(articGuppyPlex.out.fastq.toList(), ch_irida)

       generateFastaIridaReport(articMinIONNanopolish.out.consensus_fasta.collect(),
                                ch_irida)
      }
      else {
       Channel.fromPath("${params.irida}")
              .set{ ch_irida }

       accountReadFilterFailures(articGuppyPlex.out.fastq.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect())

       articMinIONNanopolish(articGuppyPlex.out.fastq
                                          .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                          .combine(articDownloadScheme.out.scheme)
                                          .combine(ch_fast5Pass)
                                          .combine(ch_seqSummary))
      }

      articRemoveUnmappedReads(articMinIONNanopolish.out.mapped)

      if (params.correctN) {
        correctFailNs(articMinIONNanopolish.out.ptrim
                          .join(articMinIONNanopolish.out.ptrimbai, by:0)
                          .join(articMinIONNanopolish.out.consensus_fasta, by:0)
                          .join(articMinIONNanopolish.out.fail_vcf, by:0),
                          articDownloadScheme.out.reffasta)

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
      
      snpDists(runNcovTools.out.aligned)

      makeQCCSV(articMinIONNanopolish.out.ptrim
                                     .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                                     .join(articMinIONNanopolish.out.vcf, by: 0)
                                     .combine(articDownloadScheme.out.reffasta)
                                     .combine(runNcovTools.out.lineage)
                                     .combine(runNcovTools.out.ncovtools_qc)
                                     .combine(runNcovTools.out.ncovtools_negative)
                                     .combine(ch_irida)
                                     .combine(runNcovTools.out.snpeff_path),
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

     correctQCSummaryCSV(writeQCSummaryCSV.out)

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONNanopolish.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))

     if (params.outCram) {
        bamToCram(articMinIONNanopolish.out.ptrim.map{ it[0] } 
                        .join (articDownloadScheme.out.reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }

     if (params.irida) {
       if (params.upload_irida) {
         Channel.fromPath("${params.upload_irida}")
             .set{ ch_upload }

         generateFast5IridaReport(ch_fast5Pass, ch_irida)

         uploadIrida(generateFastqIridaReport.out, generateFastaIridaReport.out, generateFast5IridaReport.out, ch_upload, correctQCSummaryCSV.out)

         if (params.correctN) {
           uploadCorrectN(correctFailNs.out.corrected_consensus.collect(),
                          ch_upload,
                          ch_irida)
         }
       }
     }

    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONNanopolish.out.vcf

}

workflow sequenceAnalysisMedaka {
    take:
      ch_runFastqDirs

    main:
      articDownloadScheme()

      articGuppyPlex(ch_runFastqDirs.flatten())

      if (params.irida) {
       Channel.fromPath("${params.irida}")
              .set{ ch_irida }

       renameSamples(articGuppyPlex.out.fastq
                                       .combine(ch_irida))

       accountReadFilterFailures(renameSamples.out.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect())

       articMinIONMedaka(renameSamples.out
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(articDownloadScheme.out.scheme))

       generateFastqIridaReport(articGuppyPlex.out.fastq.toList(), ch_irida)

       generateFastaIridaReport(articMinIONNanopolish.out.consensus_fasta.collect(),
                                ch_irida)
      }
      else {
       Channel.fromPath("${params.irida}")
              .set{ ch_irida }

       accountReadFilterFailures(articGuppyPlex.out.fastq.filter{ it.countFastq() <= params.minReadsArticGuppyPlex }.collect())

       articMinIONMedaka(articGuppyPlex.out.fastq
                                      .filter{ it.countFastq() > params.minReadsArticGuppyPlex }
                                      .combine(articDownloadScheme.out.scheme))
      }

      articRemoveUnmappedReads(articMinIONMedaka.out.mapped)

      if (params.correctN) {
        correctFailNs(articMinIONMedaka.out.ptrim
                          .join(articMinIONMedaka.out.ptrimbai, by:0)
                          .join(articMinIONMedaka.out.consensus_fasta, by:0)
                          .join(articMinIONMedaka.out.fail_vcf, by:0),
                          articDownloadScheme.out.reffasta)

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
      
      snpDists(runNcovTools.out.aligned)

      makeQCCSV(articMinIONMedaka.out.ptrim
                                     .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                                     .join(articMinIONMedaka.out.vcf, by: 0)
                                     .combine(articDownloadScheme.out.reffasta)
                                     .combine(runNcovTools.out.lineage)
                                     .combine(runNcovTools.out.ncovtools_qc)
                                     .combine(runNcovTools.out.ncovtools_negative)
                                     .combine(ch_irida)
                                     .combine(runNcovTools.out.snpeff_path),
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

     correctQCSummaryCSV(writeQCSummaryCSV.out)

     collateSamples(qc.pass.map{ it[0] }
                           .join(articMinIONMedaka.out.consensus_fasta, by: 0)
                           .join(articRemoveUnmappedReads.out))

     if (params.outCram) {
        bamToCram(articMinIONMedaka.out.ptrim.map{ it[0] } 
                        .join (articDownloadScheme.out.reffasta.combine(ch_preparedRef.map{ it[0] })) )

      }
    emit:
      qc_pass = collateSamples.out
      reffasta = articDownloadScheme.out.reffasta
      vcf = articMinIONMedaka.out.vcf

}


workflow articNcovNanopore {
    take:
      ch_fastqDirs
    
    main:
      if ( params.nanopolish ) {
          Channel.fromPath( "${params.fast5_pass}" )
                 .set{ ch_fast5Pass }

          Channel.fromPath( "${params.sequencing_summary}" )
                 .set{ ch_seqSummary }

          sequenceAnalysisNanopolish(ch_fastqDirs, ch_fast5Pass, ch_seqSummary)

          sequenceAnalysisNanopolish.out.vcf.set{ ch_nanopore_vcf }

          sequenceAnalysisNanopolish.out.reffasta.set{ ch_nanopore_reffasta }

      } else if ( params.medaka ) {
          sequenceAnalysisMedaka(ch_fastqDirs)

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
