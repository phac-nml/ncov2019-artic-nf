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

# Run Nanopolish Pipeline
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
