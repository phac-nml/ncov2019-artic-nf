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

process generateFast5IridaReport {

    publishDir "${params.outdir}", pattern: "irida_fast5", mode: "symlink"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    path(fast5_dirs)
    file(sampletsv)

    output:
    path("irida_fast5")

    script:
    """
    irida_fast5.py --sample_info ${sampletsv} --sample_dir ${fast5_dirs} --output_dir irida_fast5
    """
}

process runNcovTools {

    publishDir "${params.outdir}/qc_plots", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*.tsv", mode: "copy"

    //conda 'environments/ncovtools.yml'
    // Make conda env with mamba or it will error (takes 3+ hours without)

    label 'ncovtools'

    input:
    file(config)
    file(reference)
    file(amplicon)
    file(nanopolishresults)
    file(bed)
    file(metadata)

    // Currently have the nml_* outputs hardcoded as the config has the run name as nml
    // If you change the ncov-tools config change them as well in all instances below
    output:
    file("*.pdf")
    path "*.tsv"

    path "ncov-tools/lineages/*.csv" , emit: lineage
    path "nml_summary_qc.tsv" , emit: ncovtools_qc
    path "nml_negative_control_report.tsv" , emit: ncovtools_negative

    script:
    
    // Different with IRIDA param due to files being renamed and the addition of metadata
    // Touch nml_negative_control as it isn't always made and we need it even if its blank
    if ( params.irida )

        """
        bash run_ncovtools.sh ${params.negative_control} ${config} ${amplicon} ${reference} ${bed} ${metadata}
        """

    else
        """
        git clone https://github.com/jts/ncov-tools.git
        sed -i -e 's/^metadata/#metadata/' ${config}
        sed -i -e 's/^negative_control_samples/#negative_control_samples/' ${config}
        sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta
        mv ${amplicon} ./ncov-tools/input_amplicon.bed
        mv ${config} ${reference} ${bed} ./ncov-tools
        mkdir ./ncov-tools/run
        mv *.* ./ncov-tools/run
        cd ncov-tools
        samtools faidx ${reference}
        snakemake -s qc/Snakefile all_qc_sequencing --cores 8
        snakemake -s qc/Snakefile all_qc_analysis --cores 8
        snakemake -s qc/Snakefile all_qc_reports --cores 4
        mv ./plots/*.pdf ../
        mv ./qc_reports/*.tsv ../
        cd ..
        touch nml_negative_control_report.tsv
        """
}

process uploadIrida {

    //conda 'environments/irida_uploader.yml'

    label 'Upload'

    input:
    path(fastq_folder)
    path(consensus_folder)
    path(fast5_folder)
    file(irida_config)
    file(metadata_csv)

    script:
    """
    irida-uploader --config ${irida_config} -d ${fastq_folder}
    irida-uploader --config ${irida_config} -d ${consensus_folder} --upload_mode=assemblies
    irida-uploader --config ${irida_config} -d ${fast5_folder} --upload_mode=fast5
    upload.py --config ${irida_config} --metadata ${metadata_csv}
    """
}
