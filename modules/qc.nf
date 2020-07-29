process makeQCCSV {
    tag { sampleName }

    //conda 'environments/extras.txt'

    publishDir "${params.outdir}/qc_plots", pattern: "${sampleName}.depth.png", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(fasta), path(ref), path(lineage), path(sample_sheet)

    output:
    path "${params.prefix}.${sampleName}.qc.csv", emit: csv
    path "${sampleName}.depth.png"

    script:
    if ( params.illumina ) {
       qcSetting = "--illumina"
    } else {
       qcSetting = "--nanopore"
    }

    def rev = workflow.revision ?: workflow.commitId ?: workflow.scriptId

    if ( params.irida )

        """
        qc.py ${qcSetting} --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --pangolin ${lineage} --run_name ${sample_sheet} --revision ${rev}
        """
    
    else

        """
        qc.py ${qcSetting} --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --pangolin ${lineage} --revision ${rev}
        """
}


process writeQCSummaryCSV {

    publishDir "${params.outdir}", pattern: "${params.prefix}.qc.csv", mode: 'copy'

    tag { params.prefix }

    input:
    val lines

    output:
    file("${params.prefix}.qc.csv")

    exec:
    task.workDir.resolve("${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}
