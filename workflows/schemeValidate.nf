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
    // =============================== //
    // Pull scheme from online
    // =============================== //
    if ( ! ch_scheme ) {
        downloadScheme()
        ch_scheme = downloadScheme.out.scheme
    }

    validateScheme(
        ch_scheme
    )

    generateAmpliconBed(
        validateScheme.out.primer_bed
    )

    emit:
    scheme = validateScheme.out.scheme
    reference = validateScheme.out.reference
    primer_bed = validateScheme.out.primer_bed
    amplicon_bed = generateAmpliconBed.out.amplicon_bed
    primer_prefix = generateAmpliconBed.out.primer_prefix
}
