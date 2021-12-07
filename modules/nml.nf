process renameSamples {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.fastq", mode: "copy"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    tuple file(fastq), file(sampletsv)

    output:
    file('*.fastq')

    script:
    """
    rename_fastq.py --fastq ${fastq} --sample_info ${sampletsv}
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
    irida_upload_csv_generator.py --sample_info ${sampletsv} --sample_dir irida_fastq --fastq
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
    irida_upload_csv_generator.py --sample_info ${sampletsv} --sample_dir irida_consensus --fasta
    """
}

process generateFast5IridaReport {

    publishDir "${params.outdir}", pattern: "irida_fast5", mode: "symlink"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'fast5compress'

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

process correctFailNs {

    publishDir "${params.outdir}/corrected_consensus", pattern: "*.corrected.consensus.fasta", mode: "copy"
    publishDir "${params.outdir}/corrected_consensus", pattern: "logs/*.log", mode: "copy"

    //conda 'environments/extras.yml'
    // Experimental change to re-check N calls

    label 'smallmem'

    input:
    tuple(sampleName, path(bamfile), path(bambai), path(consensus), path(fail_vcf))
    file(reference)

    output:
    path "*.corrected.consensus.fasta", optional: true, emit: corrected_consensus
    path "logs/*.log", optional: true, emit: logs

    script:
    """
    gzip -f $fail_vcf
    correct_n.py --bam $bamfile --consensus $consensus --reference $reference --fail_vcf ${fail_vcf}.gz --force
    mkdir -p logs
    mv *.log logs/
    """
}

process runNcovTools {

    publishDir "${params.outdir}/qc_plots", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*.tsv", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*aligned.fasta", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "ncov-tools/qc_annotation/*", mode: "copy"

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
    path(corrected_fastas)

    // Currently have the nml_* outputs hardcoded as the config has the run name as nml
    // If you change the ncov-tools config change them as well in all instances below
    output:
    file("*.pdf")
    path("*.tsv")

    path "ncov-tools/lineages/*.csv" , emit: lineage
    path "nml_summary_qc.tsv" , emit: ncovtools_qc
    path "nml_negative_control_report.tsv" , emit: ncovtools_negative
    path "nml_aligned.fasta" , emit: aligned
    path "ncov-tools/qc_annotation/", emit: snpeff_path

    script:
    """
    bash run_ncovtools.sh ${config} ${amplicon} ${reference} ${bed} ${metadata}
    """
}

process snpDists {

    publishDir "${params.outdir}/", pattern: "matrix.tsv", mode: "copy"

    //conda 'environments/snpdist.yml'

    label 'smallmem'

    input:
    file(aligned_fasta)

    output:
    path("matrix.tsv")

    script:
    """
    snp-dists ${aligned_fasta} > matrix.tsv
    """

}

process uploadIrida {

    //conda 'environments/irida_uploader.yml'

    publishDir "${params.outdir}", pattern: "metadata_upload_status.csv", mode: "copy"

    label 'Upload'
    errorStrategy 'terminate'

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
    upload.py --config ${irida_config} --metadata_csv  ${metadata_csv}
    irida-uploader --config ${irida_config} -d ${fast5_folder} --upload_mode=fast5
    """
}

process uploadCorrectN{

    //conda 'environments/irida_uploader.yml'

    label 'Upload'
    errorStrategy 'terminate'

    input:
    path(fastas)
    file(irida_config)
    file(sampletsv)

    script:
    """
    mkdir -p corrected_consensus
    mv *.corrected.consensus.fasta corrected_consensus
    irida_upload_csv_generator.py --sample_info ${sampletsv} --sample_dir corrected_consensus --fasta
    irida-uploader --config ${irida_config} -d corrected_consensus --upload_mode=assemblies
    """
}
