#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p ../conda_cache

### Run Medaka Pipeline Flat ###
nextflow run ./main.nf \
    -profile mamba,test \
    --cache ../conda_cache \
    --medaka \
    --medaka_model r941_min_hac_g507 \
    --prefix 'nml' \
    --basecalled_fastq $PWD/.github/data/nanopore/fastq_pass/barcode78/ \
    --scheme_version freed_V2_nml \
    --min_length 800 \
    --max_length 1600 \
    --sequencing_technology GridION \
    --irida $PWD/.github/data/metadata.tsv

### Check Outputs ###
# 1. Num Reads
READS=`awk -F, '$1 == "fastq_runid_20e5f24723fa8df117afef2fbcd04f61f053acdb_0_0" {print $5}' ./results/nml.qc.csv`
if [[ "$READS" != "3350" ]]; then 
    echo "Incorrect output: Number of reads mapped"
    echo "  Expected: 3350, Got: $READS"
    exit 1
fi
READS=`awk -F, '$1 == "fastq_runid_20e5f24723fa8df117afef2fbcd04f61f053acdb_1_0" {print $5}' ./results/nml.qc.csv`
if [[ "$READS" != "3419" ]]; then 
    echo "Incorrect output: Number of reads mapped"
    echo "  Expected: 3419, Got: $READS"
    exit 1
fi
# 2. Number Consensus Ns
N_COUNT=`awk -F, '$1 == "fastq_runid_20e5f24723fa8df117afef2fbcd04f61f053acdb_0_0" {print $6}' ./results/nml.qc.csv`
if [[ "$N_COUNT" != "1238" ]]; then 
    echo "Incorrect output: Number Consensus Ns"
    echo "  Expected: 1238, Got: $N_COUNT"
    exit 1
fi
N_COUNT=`awk -F, '$1 == "fastq_runid_20e5f24723fa8df117afef2fbcd04f61f053acdb_1_0" {print $6}' ./results/nml.qc.csv`
if [[ "$N_COUNT" != "1238" ]]; then 
    echo "Incorrect output: Number Consensus Ns"
    echo "  Expected: 1238, Got: $N_COUNT"
    exit 1
fi

# Reset
rm -rf ./results ./work/ .nextflow*
echo "Passed Test"
