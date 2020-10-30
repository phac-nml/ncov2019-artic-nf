#!/bin/bash

# Get the needed variables from the pipeline
# Names are the same to avoid confusion!
negativeControlGrep=$1
config=$2
amplicon=$3
reference=$4
bed=$5
metadata=$6

# Moves corrected consensus data to the right header name and replaces the non-corrected ones
sed -i "s|-updated||" *.corrected.consensus.fasta

# Overwrites the original consensus files so that we use the corrected ones
for i in *.corrected.consensus.fasta
    do
        name="${i%%.*}"
        echo $name
        mv $i ${name}.consensus.fasta
    done

# Clone in ncov-tools
git clone https://github.com/jts/ncov-tools.git

# Remove the /ARTIC/nanopolish from the files or it seems to fail
sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta

# If we have a matching negative control, we modify the config to make sure its gotten
# Else, we comment out the negative controls if its not found
if $(ls | grep -q "${negativeControlGrep}")
then
    run=$(awk '{if(NR==2){ print $2; }}' < $metadata)
    sed -i -e "s/PLACEHOLDER/$run/g" $config
else
    sed -i -e 's/^negative_control_samples/#negative_control_samples/' ${config}
fi

# Move files into the correct spots
mv ${amplicon} ./ncov-tools/input_amplicon.bed
mv ${config} ${reference} ${bed} ./ncov-tools
mv ${metadata} ./ncov-tools/metadata.tsv
mkdir ./ncov-tools/run
mv *.* ./ncov-tools/run

# Go in, run the commands and generate the indexed reference sequence
cd ncov-tools
samtools faidx ${reference}
snakemake -s workflow/Snakefile all --cores 8

# Move files out so that they can be easily detected by nextflow
mv ./plots/*.pdf ../
mv ./qc_reports/*.tsv ../
cd ..

# Touching a negative control so that there always is one (even if we remove the check)
# To allow us to always smash together all qc outputs into one file!
touch nml_negative_control_report.tsv
