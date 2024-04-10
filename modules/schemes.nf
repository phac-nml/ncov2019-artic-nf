// Scheme related processes to go with scheme subworkflow
process downloadScheme {
    // Download primer scheme from given repo URL to primer-schemes directory
    tag { params.schemeRepoURL }

    output:
    path "primer-schemes", emit: scheme

    script:
    """
    git clone ${params.schemeRepoURL} primer-schemes
    """
}
process validateScheme {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "primer-schemes", mode: "copy"

    input:
    path scheme

    output:
    path "primer-schemes/${params.scheme}/${schemeVersion}/*.reference.fasta", emit: reference
    path "primer-schemes/${params.scheme}/${schemeVersion}/*.primer.bed", emit: primer_bed
    tuple val(schemeVersion), path("primer-schemes"), emit: scheme

    script:
    // Check for starting with a V as that is required by artic
    schemeVersion = params.schemeVersion
    def adjust_scheme_version = false
    if ( ! schemeVersion.startsWith("V") ) {
        schemeVersion = "V" + schemeVersion
        adjust_scheme_version = true
    }
    """
    # Check for local scheme name to match everything up
    if [[ "$scheme" != "primer-schemes" ]]; then
        mv $scheme primer-schemes
    fi

    # Adjust version if required
    if $adjust_scheme_version; then
        mv primer-schemes/${params.scheme}/${params.schemeVersion} primer-schemes/${params.scheme}/${schemeVersion}
    fi

    # File checks
    if [ ! -f primer-schemes/${params.scheme}/${schemeVersion}/*reference.fasta ]; then
        echo "ERROR: Reference Fasta not found in 'primer-schemes/${params.scheme}/${schemeVersion}/*reference.fasta'"
        exit 1
    elif [ ! -f primer-schemes/${params.scheme}/${schemeVersion}/*scheme.bed ]; then
        echo "ERROR: Scheme bed file not found in 'primer-schemes/${params.scheme}/${schemeVersion}/*primer.bed'"
        exit 1
    fi
    """
}
process generateAmpliconBed {
    publishDir "${params.outdir}/amplicon_info", pattern: "*.bed", mode: "copy"
    publishDir "${params.outdir}/amplicon_info", pattern: "primer_prefix.txt", mode: "copy"

    input:
    path primer_bed

    output:
    path "amplicon.bed", emit: amplicon_bed
    path "tiling_region.bed", emit: tiling_bed
    env PREFIX, emit: primer_prefix

    script:
    """
    primers_to_amplicons.py \\
        --bed $primer_bed

    PREFIX=\$(cat primer_prefix.txt)
    """
}
