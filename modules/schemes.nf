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
    path "primer-schemes/${params.scheme}/${schemeFolder}/*.reference.fasta", emit: reference
    path "primer-schemes/${params.scheme}/${schemeFolder}/*.primer.bed", emit: primer_bed // primer.bed for now as the repos have some conflicting scheme.bed
    tuple val(schemeVersion), path("primer-schemes"), emit: scheme

    script:
    /*
        ARTIC Requires the following:
            - Scheme folder starting with a V (ex. Vfreed/)
            - Scheme input not starting with a V (ex. freed)

        So a few checks to attempt to meet those requirements
        Note that no def for the params as we need them in the outputs
    */
    schemeVersion = params.schemeVersion
    schemeFolder = params.schemeVersion
    if ( ! schemeVersion.startsWith("V") ) {
        schemeFolder = "V" + schemeFolder
    } else {
        schemeVersion = schemeVersion.substring(1)
    }
    """
    # Check for local scheme name to match everything up
    if [[ "$scheme" != "primer-schemes" ]]; then
        mv $scheme primer-schemes
    fi

    # Adjust folder if required
    if [ ! -d primer-schemes/${params.scheme}/${schemeFolder} ]; then
        if [ -d primer-schemes/${params.scheme}/${schemeVersion} ]; then
            mv primer-schemes/${params.scheme}/${schemeVersion} primer-schemes/${params.scheme}/${schemeFolder}
        else
            echo "ERROR: Cannot find input scheme version ${schemeVersion}"
            exit 1
        fi
    fi

    # File checks
    if [ ! -f primer-schemes/${params.scheme}/${schemeFolder}/*reference.fasta ]; then
        echo "ERROR: Reference Fasta not found in 'primer-schemes/${params.scheme}/${schemeFolder}/*reference.fasta'"
        exit 1
    elif [ ! -f primer-schemes/${params.scheme}/${schemeFolder}/*scheme.bed ]; then
        echo "ERROR: Scheme bed file not found in 'primer-schemes/${params.scheme}/${schemeFolder}/*primer.bed'"
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
