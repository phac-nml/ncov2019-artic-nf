process nextcladeDatasetGet {
    tag { version_tag }
    label 'smallmem'
    label 'nextclade'
    publishDir "${params.outdir}/nextclade", pattern: "$dataset", mode: "copy"

    input:
    val dataset
    val version_tag

    output:
    path "$dataset", emit: dataset
    path "versions.yml", emit: versions

    script:
    def tagVersion = version_tag ? "--tag ${version_tag}" : ""
    """
    nextclade dataset get \\
        --name $dataset \\
        ${tagVersion} \\
        --output-dir $dataset

    # Versions #
    # From Katherine's tracking
    name=\$(grep -m 1 "name" ${dataset}/pathogen.json | sed 's/ \\|,\\|"//g' | cut -d ":" -f 2)
    tag=\$(grep -m 1 "tag" ${dataset}/pathogen.json | sed 's/ \\|,\\|"//g' | cut -d ":" -f 2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
        nextclade_db: \$name \$tag
    END_VERSIONS
    """
}

process nextcladeRun {
    tag { sampleName }
    label 'smallmem'
    label 'nextclade'
    publishDir "${params.outdir}/nextclade", pattern: "${sampleName}_nextclade.tsv", mode: "copy"

    input:
    tuple val(sampleName), path(fasta)
    path dataset

    output:
    tuple val(sampleName), path("${sampleName}_nextclade.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    # Correct header for nc output line
    sed -iE 's|/ARTIC/.*||g' $fasta

    # Run nc
    nextclade run \\
        -D $dataset \\
        -t ${sampleName}_nextclade.tsv \\
        -j ${task.cpus} \\
        $fasta

    # Versions #
    # From Katherine's tracking
    name=\$(grep -m 1 "name" ${dataset}/pathogen.json | sed 's/ \\|,\\|"//g' | cut -d ":" -f 2)
    tag=\$(grep -m 1 "tag" ${dataset}/pathogen.json | sed 's/ \\|,\\|"//g' | cut -d ":" -f 2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
        nextclade_db: \$name \$tag
    END_VERSIONS
    """
}
