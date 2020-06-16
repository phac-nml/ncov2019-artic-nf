process generateIridaReport {

    publishDir "${params.outdir}", pattern: "irida_upload", mode: "copy"

    //conda 'environments/extras.txt'

    input:
    file(fastq)
    file(samplecsv)

    output:
    path("irida_upload")

    script:
    """
    mkdir irida_upload
    mv ${fastq} irida_upload
    irida_samples.py --sample_info ${samplecsv} --prefix ${params.prefix} --sample_dir irida_upload
    """
}

process runNcovTools {

    publishDir "${params.outdir}/irida_upload", pattern: "*.pdf", mode: "copy"

    //conda 'environments/ncovtools.yml'

    input:
    file(config)
    file(reference)
    file(amplicon)
    file(nanopolishresults)

    output:
    file("*.pdf")

    script:
    """
    git clone https://github.com/jts/ncov-tools.git
    mv ${config} ${reference} ${amplicon} ./ncov-tools
    mkdir ./ncov-tools/run
    mv ${params.prefix}* ./ncov-tools/run
    cd ncov-tools
    snakemake -s qc/Snakefile all_qc_sequencing --cores 6
    mv ./plots/* ../
    """
}
