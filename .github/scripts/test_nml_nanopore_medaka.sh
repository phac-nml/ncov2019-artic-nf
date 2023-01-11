#!/bin/env/usr bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# Run Medaka Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --medaka \
    --prefix 'nml' \
    --basecalled_fastq $PWD/.github/data/nanopore/fastq_pass/ \
    --medakaModel r941_min_hac_g507 \
    --schemeVersion freed_V2_nml \
    --min_length 800 \
    --max_length 1600 \
    --sequencingTechnology GridION \
    --schemeRepoURL 'https://github.com/DarianHole/primer-schemes.git' \
    --irida $PWD/.github/data/metadata.tsv

### Check Outputs ###
