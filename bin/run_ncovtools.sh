#!/bin/bash

# Get the needed variables from the pipeline
# Names are the same to avoid confusion!
CONFIG=$1
AMPLICON_BED=$2
NCOV_REF=$3
PRIMER_BED=$4
METADATA=$5
PRIMER_PREFIX=$6
CORES=$7

# Moves corrected consensus data to the right header name and replaces the non-corrected ones
sed -i "s|-updated||" *.corrected.consensus.fasta

# Overwrites the original consensus files so that we use the corrected ones
for corrected_consensus in *.corrected.consensus.fasta
    do
        name="${corrected_consensus%%.*}"
        echo $name
        mv $corrected_consensus ${name}.consensus.fasta
    done

# Clone in ncov-tools
git clone https://github.com/jts/ncov-tools.git

# Remove the /ARTIC/nanopolish or /ARTIC/medaka from the files
sed -i 's|/ARTIC/nanopolish||' *.consensus.fasta
sed -i 's|/ARTIC/medaka||' *.consensus.fasta

#### Moving files and modifying the config to get what we expect based on inputs ####
# cp config to allow us to mess with its values
mv ${AMPLICON_BED} ./ncov-tools/input_amplicon.bed
cp ${CONFIG} ./ncov-tools
mv ${NCOV_REF} ./ncov-tools/nCoV-2019.reference.fasta
mv ${PRIMER_BED} ./ncov-tools/nCoV-2019.bed

# If we have a matching negative control, we modify the config to make sure its gotten
# If we don't find any, then no negative controls are added
if $(ls | grep -q -i "negative\|ntc\|water\|blank")
then
   negative_list=$(grep -i -e ntc -e negative -e water -e blank -e neg ${METADATA} | cut -f 1 | sed 's/^/"/g' | sed 's/$/"/g' | tr "\n" ',' | sed 's/^/[/' | sed 's/$/]/')
   echo "negative_control_samples: ${negative_list}" >> ./ncov-tools/${CONFIG}
fi

# Add in our primerprefix
echo "primer_prefix: '$PRIMER_PREFIX'" >> ./ncov-tools/${CONFIG}

# Check for metadata file
# If irida sample sheet is used, we will have some and will move it into ncov-tools folder
# If not, then we will have false passed and will comment out the metadata line to allow ncov-tools to run
if [ -f "$METADATA" ];
then
    mv ${METADATA} ./ncov-tools/metadata.tsv
else
    sed -i -e 's/^metadata/#metadata/' ./ncov-tools/${CONFIG}
fi

# mv all the files into a folder to run on
mkdir ./ncov-tools/run
mv *.* ./ncov-tools/run


# Go in, run the commands and generate the indexed reference sequence
cd ncov-tools
samtools faidx nCoV-2019.reference.fasta
snakemake -kp -s workflow/Snakefile --cores 1 build_snpeff_db
snakemake -kp -s workflow/Snakefile all --cores ${CORES}

# Move files out so that they can be easily detected by nextflow
mv ./plots/*.pdf ../
mv ./qc_reports/*.tsv ../
mv ./qc_analysis/nml_aligned.fasta ../
mv ./qc_analysis/nml_amplicon_coverage_table.tsv ../
cd ..

# Touching a negative control so that there always is one (even if we remove the check)
# This will allow us to always smash together all qc outputs into one file!
touch nml_negative_control_report.tsv
