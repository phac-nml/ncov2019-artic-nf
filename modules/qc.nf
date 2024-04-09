process makeQCCSV {
    // Make singular sample qc.csv files
    tag { sampleName }
    label 'conda_extra'
    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(fasta), path(vcf), path(ref), path(lineage), path(ncov_summary), path(ncov_negative), path(snp_eff_path), path(scheme_bed)
    path sample_sheet
    path pcr_bed

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

    script:
    qcSetting = "--nanopore"

    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    if ( params.irida )
        """
        qc.py ${qcSetting} \
        --outfile ${params.prefix}.${sampleName}.qc.csv \
        --sample ${sampleName} \
        --ref ${ref} \
        --bam ${bam} \
        --fasta ${fasta} \
        --pangolin ${lineage} \
        --ncov_summary ${ncov_summary} \
        --ncov_negative ${ncov_negative} \
        --vcf ${vcf} \
        --sample_sheet ${sample_sheet} \
        --revision ${rev} \
        --scheme ${params.schemeVersion} \
        --scheme_bed ${scheme_bed} \
        --script_name 'nml-ncov2019-artic-nf' \
        --sequencing_technology ${params.sequencingTechnology} \
        --snpeff_tsv ${snp_eff_path}/${sampleName}_aa_table.tsv \
        --pcr_bed ${pcr_bed}
        """
    else
        """
        qc.py ${qcSetting} \
        --outfile ${params.prefix}.${sampleName}.qc.csv \
        --sample ${sampleName} \
        --ref ${ref} \
        --bam ${bam} \
        --fasta ${fasta} \
        --pangolin ${lineage} \
        --ncov_summary ${ncov_summary} \
        --ncov_negative ${ncov_negative} \
        --vcf ${vcf} \
        --revision ${rev} \
        --scheme ${params.schemeVersion} \
        --scheme_bed ${scheme_bed} \
        --script_name 'nml-ncov2019-artic-nf' \
        --sequencing_technology ${params.sequencingTechnology} \
        --snpeff_tsv ${snp_eff_path}/${sampleName}_aa_table.tsv \
        --pcr_bed ${pcr_bed}
        """
}

process writeQCSummaryCSV {
    // Concatenate individual QC files together
    tag { params.prefix }

    input:
    val lines

    output:
    path "${params.prefix}.initial.qc.csv"

    exec:
    task.workDir.resolve("${params.prefix}.initial.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

process correctQCSummaryCSV {
    // Reorder and correct columns in the concatenated QC file. Focused on neg controls and tracking failed samples
    publishDir "${params.outdir}", pattern: "${params.prefix}.qc.csv", mode: 'copy'
    tag { params.prefix }

    input:
    path initial_qc_csv
    path read_count_failure_tsv
    path read_filter_failure_tsv

    output:
    path "${params.prefix}.qc.csv"

    script:
    def failCountArg = read_count_failure_tsv ? "--count_failure_tsv $read_count_failure_tsv" : ""
    def filterCountArg = read_filter_failure_tsv ? "--filter_failure_tsv $read_filter_failure_tsv" : ""
    """
    negative_control_fixes.py \\
        $failCountArg \\
        $filterCountArg \\
        --qc_csv ${initial_qc_csv} \\
        --output_prefix ${params.prefix}
    """
}
