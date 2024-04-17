#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p ../conda_cache

### Run Medaka Pipeline ###
nextflow run ./main.nf \
    -profile mamba,test \
    --cache ../conda_cache \
    --medaka \
    --medaka_model r941_min_hac_g507 \
    --prefix 'nml' \
    --basecalled_fastq $PWD/.github/data/nanopore/fastq_pass/ \
    --schemeVersion freed_V2_nml \
    --min_length 800 \
    --max_length 1600 \
    --irida $PWD/.github/data/metadata.tsv

### Check Outputs ###
# 1. Num Reads
READS=`awk -F, '$1 == "TestSample1" {print $5}' ./results/nml.qc.csv`
if [[ "$READS" != "10130" ]]; then 
    echo "Incorrect output: Number of reads mapped"
    echo "  Expected: 10130, Got: $READS"
    exit 1
fi
# 2. Number Consensus Ns
N_COUNT=`awk -F, '$1 == "TestSample1" {print $6}' ./results/nml.qc.csv`
if [[ "$N_COUNT" != "189" ]]; then 
    echo "Incorrect output: Number Consensus Ns"
    echo "  Expected: 189, Got: $N_COUNT"
    exit 1
fi

# Reset
rm -rf ./results ./work/ .nextflow*
echo "Passed Test"
