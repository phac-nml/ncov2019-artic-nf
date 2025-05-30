// ARTIC processes
process checkFastqForModel {
    // Check if artic will be able to auto-detect the model from the fastq headers OR model parameter has been set
    label 'smallmem'
    tag { sampleName }
    errorStrategy = 'terminate'

    input:
    tuple val(sampleName), path(fastq)

    output:
    tuple val(sampleName), val(true), emit: check_done 

    script:
    """
    # Have to check if the input is a directory first
    test_fastq_file="$fastq"
    if [ -d $fastq ]; then
        test_fastq_file="\$(find $fastq/ -type f -name *.fastq* | head -n 1)"
    fi

    if [[ "\$test_fastq_file" == *.gz ]]; then
        header=\$(zgrep -m 1 '^@' "\$test_fastq_file" || echo "")
    else
        header=\$(grep -m 1 '^@' "\$test_fastq_file" || echo "")
    fi

    if [[ "\$header" == *"basecall_model_version_id"* ]] && [[ -z "${params.clair3_model}" || "${params.clair3_model}" == "null" ]]; then
        echo "FastQ header for $sampleName contains basecall model information for Clair3 model selection, artic will choose clair3 model automatically."
    elif [[ -z "${params.clair3_model}" || "${params.clair3_model}" == "null" ]]; then
        echo "ERROR: No Clair3 model provided and no basecall model found in $sampleName FastQ header!" >&2
        echo "Please make sure your input files have basecall model information or specify which Clair3 model to use with --clair3_model." >&2
        exit 1
    else
        echo "Using Clair3 model: ${params.clair3_model}"
    fi
    """
}

process articDownloadModels {
    // Pulls r10 models for clair3, models are saved here by default: $CONDA_PREFIX/bin/models
    label 'smallmem'
    script:
    """
    artic_get_models
    """
}

process articGuppyPlex {
    // Filter reads based on given length
    //  Length should be based on amplicon size
    tag { sampleName }
    label 'mediumcpu'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${newSampleName}*.fastq", mode: "copy"

    input:
    tuple val(sampleName), val(check_done), path(fastq)

    output:
    tuple val(newSampleName), path("${newSampleName}.fastq"), emit: fastq
    path "versions.yml", emit: versions

    script:
    // Fastq input can either be a directory or a set of fastq files
    //  Outputs are the same then after allowing a streamlined pipeline
    if ( fastq.isDirectory() ) {
        // For directories rename with the prefix
        newSampleName = params.prefix + '_' + sampleName
        """
        artic guppyplex \\
            --min-length ${params.min_length} \\
            --max-length ${params.max_length} \\
            --output ${newSampleName}.fastq \\
            --directory $fastq

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    } else {
        // For flat files just keep the name
        //  Still need to set the newSampleName to be added as part of the output tuple
        newSampleName = sampleName
        """
        mkdir -p input_fastq
        mv $fastq input_fastq/
        artic guppyplex \\
            --min-length ${params.min_length} \\
            --max-length ${params.max_length} \\
            --output ${newSampleName}.fastq \\
            --directory input_fastq

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    }
}

process articMinION {
    tag { sampleName }
    label 'mediumcpu'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

    input:
    tuple val(sampleName), path(fastq)
    path reference
    path primer_bed

    output:
    path "${sampleName}*", emit: all
    
    // Individual files
    tuple val(sampleName), path("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple val(sampleName), path("${sampleName}.primertrimmed.rg.sorted.bam.bai"), emit: ptrimbai
    tuple val(sampleName), path("${sampleName}.sorted.bam"), emit: mapped
    tuple val(sampleName), path("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple val(sampleName), path("${sampleName}.pass.vcf.gz"), emit: vcf
    tuple val(sampleName), path("${sampleName}.fail.vcf"), emit: fail_vcf
    path "versions.yml", emit: versions

    script:
    // Clair3 model is added conditonally if it's been set
    // Nextflow parameters to minion args
    def argsList = []
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    } else {
        argsList.add("--normalise 0")
    }
    if ( params.no_frameshift ) {
        argsList.add("--no-frameshifts")
    }
    // If no model, then will detect and use builtin one if it can and fail if not
    //  Otherwise if we have a directory as the model, we split the path
    //  and if we don't we use it as a model
    if ( params.clair3_model && params.clair3_model != 'null') {
        def modelFile = new File(params.clair3_model)
        if ( modelFile.isDirectory() ) {
            def modelPath = modelFile.getParent()
            def modelDir  = modelFile.getName()
            argsList.add("--model-dir ${modelPath}")
            argsList.add("--model ${modelDir}")
        } else {
            argsList.add("--model ${params.clair3_model}")
        }
    }
    def finalArgsConfiguration = argsList.join(" ")
    """
    artic minion \\
        ${finalArgsConfiguration} \\
        --threads ${task.cpus} \\
        --ref $reference \\
        --bed $primer_bed \\
        --read-file $fastq \\
        $sampleName

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        bcftools: \$(echo \$(bcftools --version | grep bcftools | sed 's/bcftools //'))
        minimap2: \$(echo \$(minimap2 --version))
        samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
        clair3: \$(echo \$(run_clair3.sh --version | sed 's/Clair3 v//g'))
    END_VERSIONS
    """
}
