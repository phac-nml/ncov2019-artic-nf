process renameBarcodeSamples {
    // Rename barcoded fastq samples to their name in the samplesheet
    //  Non-barcoded samples will have this step run but it'll just pass back the input
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.fastq", mode: "copy"
    label 'conda_extra'

    input:
    tuple val(sampleName), path(fastq)
    path samplesheet_tsv

    output:
    path "*.fastq", includeInputs: true, emit: fastq
    path "versions.yml", emit: versions

    script:
    """
    rename_fastq.py \\
        --fastq $fastq \\
        --metadata $samplesheet_tsv \\
        --barcode ${sampleName}

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

process accountNoReadsInput {
    // Account for no reads in the input barcode folder and rename it --irida given
    publishDir "${params.outdir}/", pattern: "samples_failing_no_input_reads.tsv", mode: "copy"

    input:
    path reads
    path samplesheet_tsv

    output:
    path "samples_failing_no_input_reads.tsv", optional: true, emit: count_filter

    // Use shell block so that we don't have to escape bash variables
    shell:
    '''
    # If we have a samplesheet get the correct columns we need for required info
    if [ -f "!{samplesheet_tsv}" ]; then 
        barcode_col=$(awk -v RS='\t' '/barcode/{print NR; exit}' "!{samplesheet_tsv}")
        sample_col=$(awk -v RS='\t' '/sample/{print NR; exit}' "!{samplesheet_tsv}")

        ## Set header
        header=$(head -n 1 "!{samplesheet_tsv}")
        echo "$header	qc_pass	nextflow_qc_pass" > samples_failing_no_input_reads.tsv
    
        ## Inputs could be either a directory called barcode## or a fastq file called *.fastq* thus have to have info for both
        for read_input in !{reads}; do
            ## Barcode## dir
            if [ -d $read_input ]; then
                barcode_n="${read_input//[!0-9]/}"
                ## Awk line uses the column number of the barcode column found and checks that it matches to the barcode number
                fileline=$(awk -F'\t' -v col="$barcode_col" -v barcode_n="$barcode_n"  '$col == barcode_n' "!{samplesheet_tsv}")
                ## No matches, skip line
                if [ "$fileline" != "" ]; then
                    echo "$fileline	TOO_FEW_INPUT_READS	FALSE" >> samples_failing_no_input_reads.tsv
                fi
            ## Fastq file
            else
                ## Removes all extensions to get the sample name, "." are not allowed in IRIDA sample names anyway
                filename="${read_input%%.*}"
                ## Use AWK to get the column as grep is having issues with extra data (barcode not in samplesheet, 'extra_nml_barcode' samples)
                fileline=$(awk -F'\t' -v col="$sample_col" -v filename="$filename"  '$col == filename' "!{samplesheet_tsv}")
                ## No matches, skip line
                if [ "$fileline" != "" ]; then
                    echo "$fileline	TOO_FEW_INPUT_READS	FALSE" >> samples_failing_no_input_reads.tsv
                fi
            fi
        done
    ## No samplesheet, just track that the input failed
    else
        echo "sample	qc_pass	nextflow_qc_pass" > samples_failing_no_input_reads.tsv
        for read_input in !{reads}; do
            if [ -d $read_input ]; then
                echo "!{params.prefix}_${read_input}	TOO_FEW_INPUT_READS	FALSE" >> samples_failing_no_input_reads.tsv
            else
                filename="${read_input%%.*}"
                echo "$filename	TOO_FEW_INPUT_READS	FALSE" >> samples_failing_no_input_reads.tsv
            fi
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
    ## If we have a samplesheet, use it to keep all of metadata the values
    if [ -f "!{samplesheet_tsv}" ]; then 
        ## Make sure we have sample column
        sample_col=$(awk -v RS='\t' '/sample/{print NR; exit}' "!{samplesheet_tsv}")
        if [ "$sample_col" == "" ]; then
            echo "ERROR: Column 'sample' does not exist and is required"
            exit 1
        fi

        ## Set header
        header=$(head -n 1 "!{samplesheet_tsv}")
        echo "$header	qc_pass	nextflow_qc_pass" > samples_failing_read_size_filter.tsv

        for fastq in *.fastq; do
            ## Removes all extensions to get the sample name, "." are not allowed in IRIDA sample names anyway
            filename="${fastq%%.*}"
            ## Use AWK to get the column as grep is having issues with extra data (barcode not in samplesheet, 'extra_nml_barcode' samples)
            fileline=$(awk -F'\t' -v col="$sample_col" -v filename="$filename"  '$col == filename' "!{samplesheet_tsv}")
            ## No matches, skip line
            if [ "$fileline" != "" ]; then
                echo "$fileline	TOO_FEW_SIZE_SELECTED_READS	FALSE" >> samples_failing_read_size_filter.tsv
            fi
        done
    ## No samplesheet, just track that the barcode failed
    ##  No prefix here as artic guppyplex size selection steps adds it
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
    path "versions.yml", emit: versions

    script:
    """
    # Create folder for upload
    mkdir irida_fastq
    mv $fastqs irida_fastq

    # Create SampleList.csv
    irida_upload_csv_generator.py \\
        --sample_info $samplesheet_tsv \\
        --sample_dir irida_fastq \\
        --fastq

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
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
    path "versions.yml", emit: versions

    script:
    """
    # Create folder for upload
    mkdir irida_consensus
    mv $fastas irida_consensus

    # Create SampleList.csv
    irida_upload_csv_generator.py \\
        --sample_info $samplesheet_tsv \\
        --sample_dir irida_consensus \\
        --fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
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
    path "versions.yml", emit: versions

    script:
    """
    irida_fast5.py \\
        --sample_info $samplesheet_tsv \\
        --sample_dir $fast5_dirs \\
        --output_dir irida_fast5

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas as pd; print(pd.__version__)")
        pigz: \$(pigz --version 2>&1 | tee | sed 's/pigz //g')
    END_VERSIONS
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
    path "versions.yml", emit: versions

    script:
    """
    # Zip fail vcf so it can be read by script
    gzip -f $fail_vcf

    # Run correct N
    correct_n.py \\
        --bam $bamfile \\
        --consensus $consensus \\
        --reference $reference \\
        --fail_vcf ${fail_vcf}.gz \\
        --force
    
    # Log output for tracking
    mkdir -p logs
    mv *.log logs/

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        correct_n.py: 0.1.0
        python: \$(python --version | sed 's/Python //g')
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
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
    path primer_bed
    path samplesheet_tsv
    path corrected_fastas
    val primer_prefix

    output:
    path("*.pdf")
    path("*.tsv")

    path "ncov-tools/lineages/*.csv", emit: lineage
    path "${params.prefix}_summary_qc.tsv", emit: ncovtools_qc
    path "${params.prefix}_negative_control_report.tsv", emit: ncovtools_negative
    path "ncov-tools/qc_annotation/", emit: snpeff_path
    path "${params.prefix}_aligned.fasta", optional: true, emit: aligned // Optional as <2 samples will not always create
    path "versions.yml", emit: versions

    script:
    """
    bash run_ncovtools.sh \\
        $config \\
        $amplicon_bed \\
        $reference \\
        $primer_bed \\
        ${primer_prefix} \\
        ${params.prefix} \\
        ${task.cpus} \\
        $samplesheet_tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augur: \$(echo \$(augur --version | sed 's/augur //'))
        bedtools: \$(echo \$(bedtools --version | sed 's/bedtools v//'))
        minimap2: \$(echo \$(minimap2 --version))
        pangolin: \$(echo \$(pangolin --version | sed 's/pangolin //'))
        samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
        snakemake: \$(echo \$(snakemake --version))
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
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
    path "versions.yml", emit: versions

    script:
    """
    snp-dists $aligned_fasta > matrix.tsv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
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
    path irida_config
    path pipeline_data_csv

    output:
    path "metadata_upload_status.csv", emit: upload_metadata
    path "versions.yml", emit: versions

    script:
    """
    irida-uploader --config $irida_config -d $fastq_folder
    irida-uploader --config $irida_config -d $consensus_folder --upload_mode=assemblies
    upload.py --config $irida_config --metadata_csv  $pipeline_data_csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
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

    output:
    path "versions.yml", emit: versions

    script:
    """
    # Create folder for upload
    mkdir -p corrected_consensus
    mv $corrected_fastas corrected_consensus

    # Create SampleList.csv
    irida_upload_csv_generator.py \\
        --sample_info $samplesheet_tsv \\
        --sample_dir corrected_consensus \\
        --fasta

    # Upload
    irida-uploader --config $irida_config -d corrected_consensus --upload_mode=assemblies

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irida-uploader: \$(echo \$(irida-uploader --version | sed 's/IRIDA Uploader //'))
    END_VERSIONS
    """
}

process outputVersions {
    // Output final version information into one file
    label 'smallmem'
    publishDir "${params.outdir}", pattern: "process_versions.yml", mode: "copy"

    input:
    path versions

    output:
    path "process_versions.yml"

    script:
    """
    cat <<-END_VERSIONS > process_versions.yml
    "Workflow":
        ${workflow.manifest.name}: ${workflow.manifest.version}
        Nextflow: ${nextflow.version}
    END_VERSIONS

    cat $versions >> process_versions.yml
    """
}
