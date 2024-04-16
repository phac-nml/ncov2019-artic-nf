#!/bin/bash
# Basic script to organize inputs and run ncov-tools

## Variables
# Get the needed variables from the pipeline
CONFIG=$1
AMPLICON_BED=$2
REFERENCE=$3
PRIMER_BED=$4
METADATA=$5
PRIMER_PREFIX=$6
RUN_PREFIX=$7
CORES=$8

## Corrected consensus (if any)
# Updates the corrected consensus header lines to match the normal output
sed -i "s|-updated||" *.corrected.consensus.fasta

# Overwrites the original consensus files so that we use the corrected ones
for corrected_consensus in *.corrected.consensus.fasta
do
    # Names match fully before the first `.` so use that
    name="${corrected_consensus%%.*}"
    echo $name
    mv $corrected_consensus ${name}.consensus.fasta
done

## Setup
# Clone in ncov-tools
git clone https://github.com/jts/ncov-tools.git

# Remove the /ARTIC/nanopolish or /ARTIC/medaka from fasta file headers (for matching later)
sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta
sed -i 's|/ARTIC/medaka||' *.consensus.fasta

# Move files to their correct spot
#  cp config to allow us to mess with its values
mv ${AMPLICON_BED} ./ncov-tools/input_amplicon.bed
cp ${CONFIG} ./ncov-tools
mv ${REFERENCE} ./ncov-tools/nCoV-2019.reference.fasta
mv ${PRIMER_BED} ./ncov-tools/nCoV-2019.bed

## Adjust config for what we find
# Add in negative control names if we find any matching our pattern
if $(ls | grep -q -i "negative\|neg\|ntc\|water\|blank")
then
   negative_list=$(grep -i -e ntc -e negative -e water -e blank -e neg ${METADATA} | cut -f 1 | sed 's/^/"/g' | sed 's/$/"/g' | tr "\n" ',' | sed 's/^/[/' | sed 's/$/]/')
   echo "negative_control_samples: ${negative_list}" >> ./ncov-tools/${CONFIG}
fi

# Add in our primer prefix and run prefix
echo "run_name: $RUN_PREFIX" >> ./ncov-tools/${CONFIG}
echo "primer_prefix: '$PRIMER_PREFIX'" >> ./ncov-tools/${CONFIG}

# Adjust for metadata file if given
#  If we pass one, we add it to the ncov-tools folder
#  If we don't, then we comment out that config line so it isn't run
if [ -f "$METADATA" ];
then
    mv ${METADATA} ./ncov-tools/metadata.tsv
else
    sed -i -e 's/^metadata/#metadata/' ./ncov-tools/${CONFIG}
fi

## Add all the files in
# Move all the artic minion files into a folder to run on
mkdir ./ncov-tools/run
mv *.* ./ncov-tools/run

## Run ncov-tools
# Go into folder, run the commands and generate the indexed reference sequence
cd ncov-tools
samtools faidx nCoV-2019.reference.fasta
snakemake -kp -s workflow/Snakefile --cores 1 build_snpeff_db
snakemake -kp -s workflow/Snakefile all --cores ${CORES}

## Postprocessing
# Move files out so that they can be more easily viewable in the output nextflow results folder
mv ./plots/*.pdf ../
mv ./qc_reports/*.tsv ../
mv ./qc_analysis/${RUN_PREFIX}_amplicon_coverage_table.tsv ../
# Less than 2 samples may not create MSA and will not create a tree (iqtree needs 3)
if [ -f ./qc_analysis/${RUN_PREFIX}_aligned.fasta ]; then
    mv ./qc_analysis/${RUN_PREFIX}_aligned.fasta ../
fi
cd ..

# Touching a negative control so that there always is one
#   This just allows an easier control of flow when combining files
touch ${RUN_PREFIX}_negative_control_report.tsv
