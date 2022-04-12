// ARTIC processes

process articDownloadScheme{
    tag params.schemeRepoURL

    label 'internet'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.schemeDir}", mode: "copy"

    output:
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/${params.scheme}.reference.fasta" , emit: reffasta
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/${params.scheme}.bed" , emit: bed
    path "${params.schemeDir}/${params.scheme}/${params.schemeVersion}/ncov-qc_*.scheme.bed" , emit: ncov_amplicon
    path "${params.schemeDir}" , emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} ${params.schemeDir}
    """
}

process articGuppyPlex {
    tag { params.prefix + "-" + fastqDir }

    label 'mediumcpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${params.prefix}*.fastq", mode: "copy"

    input:
    path(fastqDir)

    output:
    path "${params.prefix}*.fastq", emit: fastq
    path "*.process.yml" , emit: versions

    script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.prefix} \
    --directory ${fastqDir}

    # Versions #
    cat <<-END_VERSIONS > artic_guppyplex.process.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process articGuppyPlexFlat {
    tag { params.prefix + "-" + fastq }

    label 'mediumcpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "out/${sampleName}.fastq", mode: "copy"

    input:
    path(fastq)

    output:
    path "out/${sampleName}.fastq", emit: fastq
    path "*.process.yml" , emit: versions

    // Utilize the sampleName to keep the name consistent
    //  Have to use out dir for the end name as otherwise the file is the same as the input
    script:
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')
    """
    mkdir -p input_fastq
    mkdir -p out
    mv $fastq input_fastq
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --prefix ${params.prefix} \
    --directory input_fastq

    mv ${params.prefix}_input_fastq.fastq out/${sampleName}.fastq

    # Versions #
    cat <<-END_VERSIONS > artic_guppyplex_flat.process.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process articMinIONMedaka {
    tag { sampleName }

    label 'mediumcpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

    input:
    tuple file(fastq), file(schemeRepo)

    output:
    file("${sampleName}*")
    
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam.bai"), emit: ptrimbai
    tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
    tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple sampleName, file("${sampleName}.pass.vcf.gz"), emit: vcf
    tuple sampleName, file("${sampleName}.fail.vcf"), emit: fail_vcf
    path "*.process.yml" , emit: versions

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion --medaka \
    ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo} \
    --read-file ${fastq} \
    --medaka-model ${params.medakaModel} \
    --scheme-version ${params.schemeVersion} \
    ${params.scheme} \
    ${sampleName}

    # Versions #
    cat <<-END_VERSIONS > artic_minion_medaka.process.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
            artic-tools: \$(artic-tools --version)
            bcftools: \$(echo \$(bcftools --version | grep bcftools | sed 's/bcftools //'))
            medaka: \$(echo \$(medaka --version | sed 's/medaka //'))
            minimap2: \$(echo \$(minimap2 --version))
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process articMinIONNanopolish {
    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

    input:
    tuple file(fastq), file(schemeRepo), file(fast5Pass), file(seqSummary)

    output:
    file("${sampleName}*")
    
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
    tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam.bai"), emit: ptrimbai
    tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
    tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta
    tuple sampleName, file("${sampleName}.pass.vcf.gz"), emit: vcf
    tuple sampleName, file("${sampleName}.fail.vcf"), emit: fail_vcf
    path "*.process.yml" , emit: versions

    script:
    // Make an identifier from the fastq filename
    sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

    // Configure artic minion pipeline
    minionRunConfigBuilder = []

    if ( params.normalise ) {
    minionRunConfigBuilder.add("--normalise ${params.normalise}")
    }
    
    if ( params.bwa ) {
    minionRunConfigBuilder.add("--bwa")
    } else {
    minionRunConfigBuilder.add("--minimap2")
    }

    minionFinalConfig = minionRunConfigBuilder.join(" ")

    """
    artic minion ${minionFinalConfig} \
    --threads ${task.cpus} \
    --scheme-directory ${schemeRepo} \
    --read-file ${fastq} \
    --fast5-directory ${fast5Pass} \
    --sequencing-summary ${seqSummary} \
    --scheme-version ${params.schemeVersion} \
    ${params.scheme} \
    ${sampleName}

    # Versions #
    cat <<-END_VERSIONS > artic_minion_nanopolish.process.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
            artic-tools: \$(artic-tools --version)
            bcftools: \$(echo \$(bcftools --version | grep bcftools | sed 's/bcftools //'))
            nanopolish: \$(echo \$(nanopolish --version | grep nanopolish | sed 's/nanopolish version //'))
            minimap2: \$(echo \$(minimap2 --version))
            samtools: \$(echo \$(samtools --version | head -n 1 | grep samtools | sed 's/samtools //'))
    END_VERSIONS
    """
}

process articRemoveUnmappedReads {
    tag { sampleName }

    cpus 1

    input:
    tuple(sampleName, path(bamfile))

    output:
    tuple sampleName, file("${sampleName}.mapped.sorted.bam"), emit: mapped_bam
    path "*.process.yml" , emit: versions

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
