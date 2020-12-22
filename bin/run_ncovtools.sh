#!/bin/bash

# Get the needed variables from the pipeline
# Names are the same to avoid confusion!
config=$1
amplicon=$2
reference=$3
bed=$4
metadata=$5

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

#### Moving files and modifying the config to get what we expect based on inputs ####
# cp config to allow us to mess with its values
mv ${amplicon} ./ncov-tools/input_amplicon.bed
cp ${config} ./ncov-tools
mv ${reference} ${bed} ./ncov-tools

# If we have a matching negative control, we modify the config to make sure its gotten
# If we don't find any, then no negative controls are added
if $(ls | grep -q -i "negative\|ntc\|water\|blank")
then
   negative_list=$(grep -i -e ntc -e negative -e water -e blank ${metadata} | cut -f 1 | sed 's/^/"/g' | sed 's/$/"/g' | tr "\n" ',' | sed 's/^/[/' | sed 's/$/]/')
   echo "negative_control_samples: ${negative_list}" >> ./ncov-tools/${config}
fi

# Check for metadata file
# If irida sample sheet is used, we will have some and will move it into ncov-tools folder
# If not, then we will have false passed and will comment out the metadata line to allow ncov-tools to run
if [ -f "$metadata" ];
then
    mv ${metadata} ./ncov-tools/metadata.tsv
else
    sed -i -e 's/^metadata/#metadata/' ./ncov-tools/${config}
fi

# mv all the files into a folder to run on
mkdir ./ncov-tools/run
mv *.* ./ncov-tools/run


# Go in, run the commands and generate the indexed reference sequence
cd ncov-tools
samtools faidx ${reference}
snakemake -s workflow/Snakefile all --cores 8
snakemake -s workflow/Snakefile --cores 1 build_snpeff_db
snakemake -s workflow/Snakefile --cores 2 all_qc_annotation

# Move files out so that they can be easily detected by nextflow
mv ./plots/*.pdf ../
mv ./qc_reports/*.tsv ../
mv ./qc_analysis/nml_aligned.fasta ../
cd ..

# Touching a negative control so that there always is one (even if we remove the check)
# This will allow us to always smash together all qc outputs into one file!
touch nml_negative_control_report.tsv
