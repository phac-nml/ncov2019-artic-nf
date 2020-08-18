process renameSamples {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.fastq", mode: "copy"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    tuple file(fastq), file(samplecsv)

    output:
    file('*.fastq')

    script:
    """
    irida_fastq.py --fastq ${fastq} --prefix ${params.prefix} --sample_info ${samplecsv}
    """
}

process accountNoReadsInput {

    label 'smallmem'

    publishDir "${params.outdir}/", pattern: "samples_failing_no_input_reads.txt", mode: "copy"

    input:
    path(fastq_dir)

    output:
    file('samples_failing_no_input_reads.txt')

    script:
    """
    ls -d */ > samples_failing_no_input_reads.txt
    """
}

process accountReadFilterFailures {

    label 'smallmem'

    publishDir "${params.outdir}/", pattern: "samples_failing_read_filter.txt", mode: "copy"

    input:
    path(fastq)

    output:
    file('samples_failing_read_filter.txt')

    script:
    """
    ls *.fastq > samples_failing_read_filter.txt
    """
}

process generateFastqIridaReport {

    publishDir "${params.outdir}", pattern: "irida_fastq", mode: "copy"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    file(fastqs)
    file(sampletsv)

    output:
    path("irida_fastq")

    script:
    """
    mkdir irida_fastq
    mv ${fastqs} irida_fastq
    irida_fastq.py --sample_info ${sampletsv} --prefix ${params.prefix} --sample_dir irida_fastq
    """
}

process generateFastaIridaReport {

    publishDir "${params.outdir}", pattern: "irida_consensus", mode: "copy"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    file(fastas)
    file(sampletsv)

    output:
    path("irida_consensus")

    script:
    """
    mkdir irida_consensus
    mv *.consensus.fasta irida_consensus
    irida_fasta.py --sample_info ${sampletsv} --sample_dir irida_consensus
    """
}

process runNcovTools {

    publishDir "${params.outdir}/qc_plots", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*.tsv", mode: "copy"

    //conda 'environments/ncovtools.yml'

    label 'ncovtools'

    input:
    file(config)
    file(reference)
    file(amplicon)
    file(nanopolishresults)
    file(metadata)

    output:
    file("*.pdf")
    path "*.tsv"

    path "ncov-tools/lineages/*.csv" , emit: lineage
    path "*_summary_qc.tsv" , emit: ncovtools_qc

    script:
    
    if ( params.irida )
    
        """
        git clone https://github.com/jts/ncov-tools.git
        sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta
        mv ${config} ${reference} ${amplicon} ./ncov-tools
        mv ${metadata} ./ncov-tools/metadata.tsv
        mkdir ./ncov-tools/run
        mv *.* ./ncov-tools/run
        cd ncov-tools
        samtools faidx ${reference}
        snakemake -s qc/Snakefile all_qc_sequencing --cores 8
        snakemake -s qc/Snakefile all_qc_analysis --cores 8
        snakemake -s qc/Snakefile all_qc_reports --cores 4
        mv ./plots/*.pdf ../
        mv ./qc_reports/*.tsv ../
        """
    
    else
        """
        git clone https://github.com/jts/ncov-tools.git
        sed -i -e 's/^metadata/#metadata/' ${config}
        sed -i -e 's/^negative_control_samples/#negative_control_samples/' ${config}
        sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta
        mv ${config} ${reference} ${amplicon} ./ncov-tools
        mkdir ./ncov-tools/run
        mv *.* ./ncov-tools/run
        cd ncov-tools
        samtools faidx ${reference}
        snakemake -s qc/Snakefile all_qc_sequencing --cores 8
        snakemake -s qc/Snakefile all_qc_analysis --cores 8
        snakemake -s qc/Snakefile all_qc_reports --cores 4
        mv ./plots/*.pdf ../
        mv ./qc_reports/*.tsv ../
        """
}
