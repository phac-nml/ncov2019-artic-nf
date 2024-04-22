process makeQCCSV {
    // Make singular sample qc.csv files
    tag { sampleName }
    label 'conda_extra'
    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple val(sampleName), path(bam), path(fasta), path(vcf), path(nextclade_tsv)
    path ref
    path lineage_csv
    path ncov_summary
    path ncov_negative
    path snp_eff_path
    path scheme_bed
    path samplesheet
    path pcr_bed
    val seq_tech

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"
    path "versions.yml", emit: versions

    script:
    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    def samplesheetArg = samplesheet ? "--sample_sheet $samplesheet" : ""
    def pcrBedArg = pcr_bed ? "--pcr_bed ${pcr_bed}" : ""
    """
    qc.py \\
        --nanopore \\
        --outfile ${params.prefix}.${sampleName}.qc.csv \\
        --sample ${sampleName} \\
        --ref $ref \\
        --bam $bam \\
        --fasta $fasta \\
        --pangolin $lineage_csv \\
        --ncov_summary $ncov_summary \\
        --ncov_negative $ncov_negative \\
        --vcf $vcf \\
        --revision ${rev} \\
        --scheme ${params.schemeVersion} \\
        --scheme_bed $scheme_bed \\
        --script_name 'nml-ncov2019-artic-nf' \\
        --sequencing_technology ${seq_tech} \\
        --snpeff_tsv ${snp_eff_path}/${sampleName}_aa_table.tsv \\
        --nextclade_tsv $nextclade_tsv \\
        $samplesheetArg \\
        $pcrBedArg

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
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
    path "${params.prefix}.qc.csv", emit: final_csv
    path "versions.yml", emit: versions

    script:
    def failCountArg = read_count_failure_tsv ? "--count_failure_tsv $read_count_failure_tsv" : ""
    def filterCountArg = read_filter_failure_tsv ? "--filter_failure_tsv $read_filter_failure_tsv" : ""
    """
    negative_control_fixes.py \\
        $failCountArg \\
        $filterCountArg \\
        --qc_csv $initial_qc_csv \\
        --output_prefix ${params.prefix}

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
    """
}
