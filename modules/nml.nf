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

    publishDir "${params.outdir}/", pattern: "samples_failing_no_input_reads.csv", mode: "copy"

    input:
    path(fastq_dir)
    file(sampletsv)

    output:
    file('samples_failing_no_input_reads.csv')

    // Use shell block so that we don't have to escape bash variables
    shell:
    '''
    echo "sample,qc_pass" > samples_failing_no_input_reads.csv

    if [ -f !{sampletsv} ]; then 
        for barcode in barcode*; do
            barcode_n="${barcode//[!0-9]/}"
            sample_name=$(awk -F'\t' -v barcode_n="$barcode_n"  '$3 == barcode_n' !{sampletsv} | cut -f 1)
            echo "$sample_name,TOO_FEW_INPUT_READS" >> samples_failing_no_input_reads.csv
        done
    else
        for barcode in barcode*; do
            echo "!{params.prefix}_${barcode},TOO_FEW_INPUT_READS" >> samples_failing_no_input_reads.csv
        done
    fi
    '''
}

process accountReadFilterFailures {

    label 'smallmem'

    publishDir "${params.outdir}/", pattern: "samples_failing_read_filter.csv", mode: "copy"

    input:
    path(fastq)

    output:
    file('samples_failing_read_filter.csv')

    shell:
    '''
    echo "sample,qc_pass" > samples_failing_read_filter.csv

    for fastq in *.fastq; do
        # Removes all extensions to get the sample name, . are not allowed in IRIDA sample names anyway
        filename="${fastq%%.*}"
        echo "$filename,TOO_FEW_MAPPING_READS" >> samples_failing_read_filter.csv
    done
    '''
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
    path "*.process.yml" , emit: versions

    script:
    """
    gzip -f $fail_vcf
    correct_n.py --bam $bamfile --consensus $consensus --reference $reference --fail_vcf ${fail_vcf}.gz --force
    mkdir -p logs
    mv *.log logs/

    # Versions #
    cat <<-END_VERSIONS > ncorrection.process.yml
        "${task.process}":
            correct_n.py: 0.1.0
    END_VERSIONS
    """
}

process runNcovTools {

    publishDir "${params.outdir}/qc_plots", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*.tsv", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*aligned.fasta", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "ncov-tools/qc_annotation", mode: "copy"

    //conda 'environments/ncovtools.yml'
    // Make conda env with mamba or it will error (takes 3+ hours without)

    label 'mediumcpu'

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
    path "*.process.yml" , emit: versions

    script:
    """
    bash run_ncovtools.sh ${config} ${amplicon} ${reference} ${bed} ${metadata} ${task.cpus}

    # Versions #
    cat <<-END_VERSIONS > ncovtools.process.yml
        "${task.process}":
            augur: \$(echo \$(augur --version | sed 's/augur //'))
            bedtools: \$(echo \$(bedtools --version | sed 's/bedtools v//'))
            minimap2: \$(echo \$(minimap2 --version))
            pangolin: \$(echo \$(pangolin --version | sed 's/pangolin //'))
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
            snakemake: \$(echo \$(snakemake --version))
            snpEff: \$(echo \$(snpEff -version))
    END_VERSIONS
    """
}

process snpDists {

    publishDir "${params.outdir}/", pattern: "matrix.tsv", mode: "copy"

    //conda 'environments/snpdist.yml'

    label 'smallmem'

    input:
    file(aligned_fasta)

    output:
    path "matrix.tsv", emit: matrix
    path "*.process.yml" , emit: versions

    script:
    """
    snp-dists ${aligned_fasta} > matrix.tsv

    # Versions #
    cat <<-END_VERSIONS > snpdists.process.yml
        "${task.process}":
            snp-dists: \$(echo \$(snp-dists -v | sed 's/snp-dists //'))
    END_VERSIONS
    """
}

process uploadIridaNanopolish {

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

    output:
    path "metadata_upload_status.csv", emit: upload_metadata
    path "*.process.yml" , emit: versions

    script:
    """
    irida-uploader --config ${irida_config} -d ${fastq_folder}
    irida-uploader --config ${irida_config} -d ${consensus_folder} --upload_mode=assemblies
    upload.py --config ${irida_config} --metadata_csv  ${metadata_csv}
    irida-uploader --config ${irida_config} -d ${fast5_folder} --upload_mode=fast5

    # Versions #
    cat <<-END_VERSIONS > upload_all_nanopolish.process.yml
        "${task.process}":
            irida-uploader: \$(echo \$(irida-uploader --version | sed 's/IRIDA Uploader //'))
    END_VERSIONS
    """
}

process uploadIridaMedaka {

    //conda 'environments/irida_uploader.yml'

    publishDir "${params.outdir}", pattern: "metadata_upload_status.csv", mode: "copy"

    label 'Upload'
    errorStrategy 'terminate'

    input:
    path(fastq_folder)
    path(consensus_folder)
    file(irida_config)
    file(metadata_csv)

    output:
    path "metadata_upload_status.csv", emit: upload_metadata
    path "*.process.yml" , emit: versions

    script:
    """
    irida-uploader --config ${irida_config} -d ${fastq_folder}
    irida-uploader --config ${irida_config} -d ${consensus_folder} --upload_mode=assemblies
    upload.py --config ${irida_config} --metadata_csv  ${metadata_csv}

    # Versions #
    cat <<-END_VERSIONS > upload_all_medaka.process.yml
        "${task.process}":
            irida-uploader: \$(echo \$(irida-uploader --version | sed 's/IRIDA Uploader //'))
    END_VERSIONS
    """
}

process uploadCorrectN {

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

process outputVersions {

    label 'smallmem'

    publishDir "${params.outdir}", pattern: "process_versions.yml", mode: "copy"

    input:
    path versions

    output:
    path "process_versions.yml"

    script:
    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId
    """
    echo '"articNcovNanopore-${rev}":' > process_versions.yml
    cat *.process.yml >> process_versions.yml
    """
}
