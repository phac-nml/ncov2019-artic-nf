// ARTIC processes
process articDownloadScheme{
    // Download primer scheme from given repo URL
    tag { params.schemeRepoURL }
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.schemeDir}", mode: "copy"

    output:
    path("${params.schemeDir}/${params.scheme}/${params.schemeVersion}/${params.scheme}.reference.fasta") , emit: reffasta
    path("${params.schemeDir}/${params.scheme}/${params.schemeVersion}/${params.scheme}.bed") , emit: bed
    path("${params.schemeDir}/${params.scheme}/${params.schemeVersion}/ncov-qc_*.scheme.bed") , emit: ncov_amplicon
    path("${params.schemeDir}") , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} ${params.schemeDir}
    """
}

process articGuppyPlex {
    // Filter reads based on given length
    //  Length should be based on amplicon size
    tag { sampleName }
    label 'mediumcpu'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.fastq", mode: "copy"

    input:
    tuple val(sampleName), path(fastq)

    output:
    tuple val(sampleName), path("${sampleName}.fastq"), emit: fastq
    path "*.process.yml", emit: versions

    script:
    // Fastq input can either be a directory or a set of fastq files
    //  Outputs are the same then after allowing a streamlined pipeline
    if ( fastq.isDirectory() ) {
        """
        artic guppyplex \\
            --min-length ${params.min_length} \\
            --max-length ${params.max_length} \\
            --output ${sampleName}.fastq \\
            --directory $fastq

        # Versions #
        cat <<-END_VERSIONS > artic_guppyplex.process.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    } else {
        """
        mkdir -p input_fastq
        mv $fastq input_fastq/
        artic guppyplex \\
            --min-length ${params.min_length} \\
            --max-length ${params.max_length} \\
            --output ${sampleName}.fastq \\
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
    path scheme

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
        --scheme-version ${params.schemeVersion} \\
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

process articRemoveUnmappedReads {
    // Remove reads that were not mapped to reference sequence
    tag { sampleName }

    input:
    tuple val(sampleName), path(bamfile)

    output:
    tuple val(sampleName), path("${sampleName}.mapped.sorted.bam"), emit: mapped_bam
    path "*.process.yml", emit: versions

    script:
    """
    samtools view -F4 -o ${sampleName}.mapped.sorted.bam ${bamfile}

    # Versions #
    cat <<-END_VERSIONS > artic_remove_unmapped.process.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}
