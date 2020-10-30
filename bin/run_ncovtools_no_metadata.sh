#!/bin/bash

# Get the needed variables from the pipeline
# Names are the same to avoid confusion!
config=$1
amplicon=$2
reference=$3
bed=$4

# Moves corrected consensus data to the right header name and replaces the non-corrected ones
sed -i "s|-updated||" *.corrected.consensus.fasta

# Overwrites the original consensus files so that we use the corrected ones
for i in *.corrected.consensus.fasta
    do
        name="${i%%.*}"
        echo $name
        mv $i ${name}.consensus.fasta
    done

# Get ncov-tools
git clone https://github.com/jts/ncov-tools.git

# Remove negative controls and metadata from config as we don't have it without irida data passed
sed -i -e 's/^metadata/#metadata/' ${config}
sed -i -e 's/^negative_control_samples/#negative_control_samples/' ${config}
sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta

# Move files
mv ${amplicon} ./ncov-tools/input_amplicon.bed
mv ${config} ${reference} ${bed} ./ncov-tools
mkdir ./ncov-tools/run
mv *.* ./ncov-tools/run

# Run ncovtools
cd ncov-tools
samtools faidx ${reference}
snakemake -s workflow/Snakefile all --cores 8

# Move data to be detected
mv ./plots/*.pdf ../
mv ./qc_reports/*.tsv ../

# Make the negative control report so that we don't error out
cd ..
touch nml_negative_control_report.tsv