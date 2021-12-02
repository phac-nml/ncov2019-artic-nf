#!/usr/bin/env python3

import argparse
import re
import subprocess
import pandas as pd


def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fastq',
        required=False,
        default=False,
        help='Path to single fastq input file'
    )
    parser.add_argument(
        '--sample_info',
        required=True,
        help='File path to tab separated sample sheet with the following columns ["sample", "run", "barcode", "project_id", "ct", "date"]'
    )
    return parser

def parse_sample_tsv(sample_tsv, fastq):
    '''
    Get name for fastq file to be renamed to
    '''
    BARCODE_REGEX = re.compile(r'barcode(\d+)')
    df = pd.read_csv(sample_tsv, sep='\t')

    # If we find the barcode we extract capture group 1 for the number and check the df for it
    barcode = re.search(BARCODE_REGEX, fastq)
    if barcode:
        barcode_number = int(barcode.group(1))
        df_slice = df[df['barcode'] == barcode_number]

        # If no matches found, exit it as an extra value
        if df_slice.empty:
            print('WARNING: No barcode matches {} in samplesheet {}'.format(fastq, sample_tsv))
            sample_name = 'extra_{}'.format(fastq)
        elif len(df_slice) > 1:
            print('ERROR: Multiple barcode {} in input samplesheet. Please double check them'.format(barcode_number))
            exit(1)
        else:
            # Fastq name made from the sample name and run_name
            sample_name = '{}_{}.fastq'.format(df_slice['sample'].values[0], df_slice['run'].values[0])
        return sample_name

    # If there is no barcode in the value something has gone wrong and we exit
    else:
        print('ERROR: No barcode identified in input fastq file {}'.format(fastq))
        exit(1)

def main():
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Create sample name from input TSV file using input fastq file
    sample_name = parse_sample_tsv(args.sample_info, args.fastq)
    
    # Rename input fastq file with subprocess
    subprocess.run('mv {} {}'.format(args.fastq, sample_name), shell=True)

if __name__ == "__main__":
    main()
