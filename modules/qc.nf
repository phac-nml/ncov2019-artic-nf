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

    if ( params.irida )

        """
        qc.py ${qcSetting} --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --pangolin ${lineage} --run_name ${sample_sheet} --revision $workflow.revision
        """
    
    else

        """
        qc.py ${qcSetting} --outfile ${params.prefix}.${sampleName}.qc.csv --sample ${sampleName} --ref ${ref} --bam ${bam} --fasta ${fasta} --pangolin ${lineage} --revision $workflow.revision
        """
}


process writeQCSummaryCSV {
    tag { params.prefix }

    input:
    val lines

    exec:
    file("${params.outdir}/${params.prefix}.qc.csv").withWriter { writer ->
        for ( line in lines ) {
            writer.writeLine(line.join(','))
         }   
    }
}

