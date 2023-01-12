#!/bin/env/usr bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# Set if metadata will be used or not for test
if [ "$1" = "--no_metadata" ]; then
    METADATA=""
else
    METADATA="--irida $PWD/.github/data/metadata.tsv"
fi

### Run Nanopolish Pipeline ###
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopolish \
    --prefix 'nml' \
    --basecalled_fastq $PWD/.github/data/nanopore/fastq_pass/ \
    --fast5_pass $PWD/.github/data/nanopore/fast5_pass/ \
    --sequencing_summary $PWD/.github/data/nanopore/sequencing_summary.txt \
    --schemeVersion freed_V2_nml \
    --min_length 800 \
    --max_length 1600 \
    --sequencingTechnology GridION \
    --schemeRepoURL 'https://github.com/DarianHole/primer-schemes.git' \
    $METADATA

### Check Outputs ###
# 1. Num Reads
READS=`awk -F, '$1 == "TestSample1" {print $5}' ./results/nml.qc.csv`
if [[ "$READS" != "8744" ]]; then 
    echo "Incorrect output: Number of reads mapped"
    echo "  Expected: 8744, Got: $READS"
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
