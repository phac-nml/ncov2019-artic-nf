/*
    Workflow to check, validate, and fix input primer schemes
*/
// Modules to include
include {
    downloadScheme ;
    validateScheme ;
    generateAmpliconBed
} from '../modules/schemes.nf' 

/*
    Initialize channels from params
*/
ch_scheme = params.local_scheme ? file(params.local_scheme, type: 'dir', checkIfExists: true) : []

// Workflow
workflow schemeValidate {
    main:
    // Version tracking
    ch_versions = Channel.empty()

    // =============================== //
    // Pull scheme
    // =============================== //
    if ( ! ch_scheme ) {
        downloadScheme()
        ch_scheme = downloadScheme.out.scheme
    }

    // Validate/correct scheme to fit artic minion requirements
    validateScheme(
        ch_scheme
    )

    // Generate amplicon bed file and other scheme specific info needed
    generateAmpliconBed(
        validateScheme.out.primer_bed
    )
    ch_versions = ch_versions.mix(generateAmpliconBed.out.versions)

    emit:
    scheme = validateScheme.out.scheme
    reference = validateScheme.out.reference
    primer_bed = validateScheme.out.primer_bed
    amplicon_bed = generateAmpliconBed.out.amplicon_bed
    primer_prefix = generateAmpliconBed.out.primer_prefix
    versions = ch_versions
}
