process renameSamples {
    // Rename barcoded fastq samples to their name in the --irida {input tsv} param
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.fastq", mode: "copy"
    label 'conda_extra'

    input:
    tuple val(sampleName), path(fastq)
    path samplesheet_tsv

    output:
    path "*.fastq", emit: fastq

    script:
    """
    rename_fastq.py --fastq ${fastq} --sample_info ${samplesheet_tsv}
    """
}

process accountNoReadsInput {
    // Account for no reads in the input barcode folder and rename it --irida given
    publishDir "${params.outdir}/", pattern: "samples_failing_no_input_reads.tsv", mode: "copy"

    input:
    path fastqs
    path samplesheet_tsv

    output:
    path "samples_failing_no_input_reads.tsv", optional: true, emit: count_filter

    // Use shell block so that we don't have to escape bash variables
    shell:
    '''
    ## If we have a samplesheet, we want all info along with a note in the QC Pass column carried along
    if [ -f !{samplesheet_tsv} ]; then 
        ## Make sure we have barcode column
        barcode_col=$(awk -v RS='\t' '/barcode/{print NR; exit}' !{samplesheet_tsv})
        if [ "$barcode_col" == "" ]; then
            echo "ERROR: Column 'barcode' does not exist and is required"
            exit 1
        fi

        ## Set header
        header=$(head -n 1 !{samplesheet_tsv})
        echo "$header	qc_pass	nextflow_qc_pass" > samples_failing_no_input_reads.tsv

        ## Populate
        for barcode in barcode*; do
            barcode_n="${barcode//[!0-9]/}"
            ## Awk line uses the column number of the barcode column found and checks that it matches to the barcode number
            fileline=$(awk -F'\t' -v col="$barcode_col" -v barcode_n="$barcode_n"  '$col == barcode_n' !{samplesheet_tsv})
            ## No matches, skip line
            if [ "$fileline" != "" ]; then
                echo "$fileline	TOO_FEW_INPUT_READS	FALSE" >> samples_failing_no_input_reads.tsv
            fi
        done
    ## No samplesheet, just track that the barcode failed
    else
        echo "sample	qc_pass	nextflow_qc_pass" > samples_failing_no_input_reads.tsv
        for barcode in barcode*; do
            echo "!{params.prefix}_${barcode}	TOO_FEW_INPUT_READS	FALSE" >> samples_failing_no_input_reads.tsv
        done
    fi
    '''
}

process accountReadFilterFailures {
    // Account for samples that fail after the read filtering step
    label 'smallmem'
    publishDir "${params.outdir}/", pattern: "samples_failing_read_size_filter.tsv", mode: "copy"

    input:
    path fastqs
    path samplesheet_tsv

    output:
    path "samples_failing_read_size_filter.tsv", optional: true, emit: size_filter

    shell:
    '''
    ## If we have a samplesheet, use it to keep all of the values
    if [ -f !{samplesheet_tsv} ]; then 
        ## Make sure we have sample column
        sample_col=$(awk -v RS='\t' '/sample/{print NR; exit}' !{samplesheet_tsv})
        if [ "$sample_col" == "" ]; then
            echo "ERROR: Column 'sample' does not exist and is required"
            exit 1
        fi

        ## Set header
        header=$(head -n 1 !{samplesheet_tsv})
        echo "$header	qc_pass	nextflow_qc_pass" > samples_failing_read_size_filter.tsv

        for fastq in *.fastq; do
            ## Removes all extensions to get the sample name, "." are not allowed in IRIDA sample names anyway
            filename="${fastq%%.*}"
            ## Use AWK to get the column as grep is having issues with extra data (barcode not in samplesheet, 'extra_nml_barcode' samples)
            fileline=$(awk -F'\t' -v col="$sample_col" -v filename="$filename"  '$col == filename' !{samplesheet_tsv})
            ## No matches, skip line
            if [ "$fileline" != "" ]; then
                echo "$fileline	TOO_FEW_SIZE_SELECTED_READS	FALSE" >> samples_failing_read_size_filter.tsv
            fi
        done
    ## No samplesheet, just track that the barcode failed
    else
        echo "sample	qc_pass	nextflow_qc_pass" > samples_failing_read_size_filter.tsv
        for fastq in *.fastq; do
            ## Removes all extensions to get the sample name, "." are not allowed in IRIDA sample names anyway
            filename="${fastq%%.*}"
            echo "$filename	TOO_FEW_SIZE_SELECTED_READS	FALSE" >> samples_failing_read_size_filter.tsv
        done
    fi
    '''
}

process generateFastqIridaReport {
    // Create directory for fastq files that can be uploaded to IRIDA if needed
    publishDir "${params.outdir}", pattern: "irida_fastq", mode: "copy"
    label 'smallmem'
    label 'conda_extra'

    input:
    path fastqs
    path samplesheet_tsv

    output:
    path "irida_fastq", emit: fastq_dir

    script:
    """
    mkdir irida_fastq
    mv ${fastqs} irida_fastq
    irida_upload_csv_generator.py --sample_info ${samplesheet_tsv} --sample_dir irida_fastq --fastq
    """
}

process generateFastaIridaReport {
    // Create directory for fasta files that can be uploaded to IRIDA if needed
    publishDir "${params.outdir}", pattern: "irida_consensus", mode: "copy"
    label 'smallmem'
    label 'conda_extra'

    input:
    path fastas
    path samplesheet_tsv

    output:
    path "irida_consensus", emit: fasta_dir

    script:
    """
    mkdir irida_consensus
    mv *.consensus.fasta irida_consensus
    irida_upload_csv_generator.py --sample_info ${samplesheet_tsv} --sample_dir irida_consensus --fasta
    """
}

process generateFast5IridaReport {
    // Create directory for fast5 files that can be uploaded to IRIDA.
    //  Only ran with --upload_irida as it is an intensive process
    publishDir "${params.outdir}", pattern: "irida_fast5", mode: "symlink"
    label 'fast5compress'
    label 'conda_extra'

    input:
    path fast5_dirs
    path samplesheet_tsv

    output:
    path "irida_fast5", emit: fast5_dir

    script:
    """
    irida_fast5.py --sample_info ${samplesheet_tsv} --sample_dir ${fast5_dirs} --output_dir irida_fast5
    """
}

process correctFailNs {
    // Nanopore - Correct positions that are designated as an N but can be called a reference base based on bcftools
    tag { sampleName }
    publishDir "${params.outdir}/corrected_consensus", pattern: "*.corrected.consensus.fasta", mode: "copy"
    publishDir "${params.outdir}/corrected_consensus", pattern: "logs/*.log", mode: "copy"
    label 'smallmem'
    label 'conda_extra'

    input:
    tuple val(sampleName), path(bamfile), path(bambai), path(consensus), path(fail_vcf)
    path reference

    output:
    tuple val(sampleName), path("*.corrected.consensus.fasta"), optional: true, emit: corrected_consensus
    path "logs/*.log", optional: true, emit: logs
    path "*.process.yml", emit: versions

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
    // Run ncov-tools with the bash script in the bin
    //  Script written to not have to escape a lot of variables
    publishDir "${params.outdir}/qc_plots", pattern: "*.pdf", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*.tsv", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "*aligned.fasta", mode: "copy"
    publishDir "${params.outdir}/ncov-tools_qc", pattern: "ncov-tools/qc_annotation", mode: "copy"
    label 'ncovtools'

    input:
    path config
    path reference
    path amplicon_bed
    path nanopolishresults
    path scheme_bed
    path samplesheet_tsv
    path corrected_fastas

    // Currently have the nml_* outputs hardcoded as the config has the run name as nml
    // If you change the ncov-tools config change them as well in all instances below
    output:
    path("*.pdf")
    path("*.tsv")

    path "ncov-tools/lineages/*.csv", emit: lineage
    path "nml_summary_qc.tsv", emit: ncovtools_qc
    path "nml_negative_control_report.tsv", emit: ncovtools_negative
    path "nml_aligned.fasta", emit: aligned
    path "ncov-tools/qc_annotation/", emit: snpeff_path
    path "*.process.yml", emit: versions

    script:
    """
    bash run_ncovtools.sh ${config} ${amplicon} ${reference} ${amplicon_bed} ${samplesheet_tsv} ${task.cpus}

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
    // Run snpDist check for samples in the analysis
    publishDir "${params.outdir}/", pattern: "matrix.tsv", mode: "copy"
    label 'smallmem'

    input:
    path aligned_fasta

    output:
    path "matrix.tsv", emit: matrix
    path "*.process.yml", emit: versions

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

process uploadIridaFiles {
    // Upload all data to IRIDA
    //  Includes: Fastq, Fasta, Fast5 (if nanopolish), and metadata 
    publishDir "${params.outdir}", pattern: "metadata_upload_status.csv", mode: "copy"
    label 'upload'
    label 'conda_iridaupload'
    errorStrategy 'terminate' // Don't want to duplicate samples if there is an issue

    input:
    path fastq_folder
    path consensus_folder
    path fast5_folder
    path irida_config
    path pipeline_data_csv

    output:
    path "metadata_upload_status.csv", emit: upload_metadata
    path "*.process.yml", emit: versions

    script:
    """
    # Always will be uploaded
    # -----------------------
    irida-uploader --config ${irida_config} -d ${fastq_folder}
    irida-uploader --config ${irida_config} -d ${consensus_folder} --upload_mode=assemblies
    upload.py --config ${irida_config} --metadata_csv  ${pipeline_data_csv}

    # Only if Nanopolish is run
    # -----------------------
    if [ -d ${fast5_folder} ]; then
        irida-uploader --config ${irida_config} -d ${fast5_folder} --upload_mode=fast5
    fi

    # Versions #
    cat <<-END_VERSIONS > upload_all.process.yml
        "${task.process}":
            irida-uploader: \$(echo \$(irida-uploader --version | sed 's/IRIDA Uploader //'))
    END_VERSIONS
    """
}

process uploadCorrectN {
    // Upload the N corrected consensus sequences to IRIDA
    //  If both --irida and --upload_irida params given
    label 'upload'
    label 'conda_iridaupload'
    errorStrategy 'terminate'

    input:
    path corrected_fastas
    path irida_config
    path samplesheet_tsv

    script:
    """
    mkdir -p corrected_consensus
    mv *.corrected.consensus.fasta corrected_consensus
    irida_upload_csv_generator.py --sample_info ${samplesheet_tsv} --sample_dir corrected_consensus --fasta
    irida-uploader --config ${irida_config} -d corrected_consensus --upload_mode=assemblies
    """
}

process outputVersions {
    // Output versions for each process
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
