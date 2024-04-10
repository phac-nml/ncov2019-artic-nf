// ARTIC processes
process articGuppyPlex {
    // Filter reads based on given length
    //  Length should be based on amplicon size
    tag { sampleName }
    label 'mediumcpu'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${newSampleName}*.fastq", mode: "copy"

    input:
    tuple val(sampleName), path(fastq)

    output:
    tuple val(newSampleName), path("${newSampleName}.fastq"), emit: fastq
    path "*.process.yml", emit: versions

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
        cat <<-END_VERSIONS > artic_guppyplex.process.yml
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
        cat <<-END_VERSIONS > artic_guppyplex.process.yml
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
    path fast5_dir
    path sequencing_summary
    tuple val(scheme_version), path(scheme)

    output:
    path "${sampleName}*", emit: all
    
    // Individual files
    tuple val(sampleName), path("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple val(sampleName), path("${sampleName}.primertrimmed.rg.sorted.bam.bai"), emit: ptrimbai
    tuple val(sampleName), path("${sampleName}.sorted.bam"), emit: mapped
    tuple val(sampleName), path("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple val(sampleName), path("${sampleName}.pass.vcf.gz"), emit: vcf
    tuple val(sampleName), path("${sampleName}.fail.vcf"), emit: fail_vcf
    path "*.process.yml", emit: versions

    script:
    // Setup args for medaka vs nanopolish
    def argsList = []
    if ( params.medaka ) {
        argsList.add("--medaka")
        argsList.add("--medaka-model ${params.medaka_model}")
        // Medaka only no longshot
        if ( params.no_longshot ) {
            argsList.add("--no-longshot")
        }
    } else {
        argsList.add("--fast5-directory $fast5_dir")
        argsList.add("--sequencing-summary $sequencing_summary")
    }
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    } else {
        argsList.add("--normalise 0")
    }
    if ( params.no_frameshift ) {
        argsList.add("--no-frameshifts")
    }
    def finalArgsConfiguration = argsList.join(" ")

    // Aligner
    def alignerArg = "--minimap2"
    if ( params.use_bwa ) {
        alignerArg = "--bwa"
    }

    // Cmd
    """
    artic minion \\
        ${finalArgsConfiguration} \\
        ${alignerArg} \\
        --threads ${task.cpus} \\
        --read-file $fastq \\
        --scheme-version ${scheme_version} \\
        --scheme-directory $scheme \\
        ${params.scheme} \\
        $sampleName

    # Versions #
    cat <<-END_VERSIONS > artic_minion.process.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
