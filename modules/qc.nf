process makeQCCSV {
    tag { sampleName }

    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(fasta), path(vcf), path(ref), path(lineage), path(ncov_summary), path(ncov_negative), path(sample_sheet), path(snp_eff_path)
    path(pcr_bed)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

    script:
    if ( params.illumina ) {
       qcSetting = "--illumina"
    } else {
       qcSetting = "--nanopore"
    }

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
        --sequencing_technology ${params.sequencingTechnology} \
        --snpeff_tsv ${snp_eff_path}/${sampleName}_aa_table.tsv \
        --pcr_bed ${pcr_bed}
        """
}


process writeQCSummaryCSV {

    tag { params.prefix }

    input:
    val lines

    output:
    file("${params.prefix}.initial.qc.csv")

    exec:
    task.workDir.resolve("${params.prefix}.initial.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

process correctQCSummaryCSV {

    publishDir "${params.outdir}", pattern: "${params.prefix}.qc.csv", mode: 'copy'

    tag { params.prefix }

    input:
    file(initial_qc_csv)

    output:
    file("${params.prefix}.qc.csv")

    script:
    """
    negative_control_fixes.py --qc_csv ${initial_qc_csv} --output_prefix ${params.prefix}
    """
}
