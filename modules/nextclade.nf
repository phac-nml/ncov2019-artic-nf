process nextcladeDatasetGet {
    label 'smallmem'
    label 'nextclade'
    publishDir "${params.outdir}/nextclade", pattern: "$dataset", mode: "copy"

    input:
    val dataset
    val tag

    output:
    path "$dataset", emit: dataset

    script:
    def tag_version = tag ? "--tag ${tag}" : ""
    """
    nextclade dataset get \\
        --name $dataset \\
        $tag_version \\
        --output-dir $dataset
    """
}
process nextcladeRun {
    label 'smallmem'
    label 'nextclade'
    publishDir "${params.outdir}/nextclade", pattern: "${sampleName}_nextclade.tsv", mode: "copy"

    input:
    tuple val(sampleName), path(fasta)
    path dataset

    output:
    tuple val(sampleName), path("${sampleName}_nextclade.tsv"), emit: tsv
    path "*.process.yml", emit: versions

    script:
    """
    nextclade run \\
        -D $dataset \\
        -t ${sampleName}_nextclade.tsv \\
        -j ${task.cpus} \\
        $fasta

    # From Katherine's tracking
    name=\$(grep -m 1 "name" ${dataset}/pathogen.json | sed 's/ \\|,\\|"//g' | cut -d ":" -f 2)
    tag=\$(grep -m 1 "tag" ${dataset}/pathogen.json | sed 's/ \\|,\\|"//g' | cut -d ":" -f 2)

    cat <<-END_VERSIONS > nextclade.process.yml
    "${task.process}":
        nextclade: \$(echo \$(nextclade --version 2>&1) | sed 's/^.*nextclade //; s/ .*\$//')
        nextclade_db: \$name \$tag
    END_VERSIONS
    """
}
