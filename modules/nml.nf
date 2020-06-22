process renameSamples {

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.fastq", mode: "copy"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    tuple file(fastq), file(samplecsv)

    output:
    file('*.fastq')

    script:
    """
    irida_samples.py --fastq ${fastq} --prefix ${params.prefix} --sample_info ${samplecsv}
    """
}

process generateIridaReport {

    publishDir "${params.outdir}", pattern: "irida_upload", mode: "copy"

    //conda 'environments/extras.txt'
    // Only with --irida flag

    label 'smallmem'

    input:
    file(fastqs)
    file(samplecsv)

    output:
    path("irida_upload")

    script:
    """
    mkdir irida_upload
    mv ${fastqs} irida_upload
    irida_samples.py --sample_info ${samplecsv} --prefix ${params.prefix} --sample_dir irida_upload
    """
}

process runNcovTools {

    publishDir "${params.outdir}/qc_plots", pattern: "*.pdf", mode: "copy"

    //conda 'environments/ncovtools.yml'

    label 'ncovtools'

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
    mv *.* ./ncov-tools/run
    cd ncov-tools
    snakemake -s qc/Snakefile all_qc_sequencing --cores 8
    snakemake -s qc/Snakefile all_qc_analysis --cores 8
    mv ./plots/* ../
    """
}
