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
    // --model or --model-dir based on if input is a model string or a path or nothing
    
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
    def finalArgsConfiguration = argsList.join(" ")
    """
    artic minion \\
        ${finalArgsConfiguration} \\
        --model ${params.medaka_model} \\
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
