#!/usr/bin/env python3
'''Rename input fastq file if barcode is found in metadata'''
import argparse
import csv
import re

from pathlib import Path
from typing import Optional

BARCODE_REGEX = re.compile(r'barcode(\d+)')

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--fastq',
        required=True,
        type=str,
        help='Input fastq file to rename'
    )
    parser.add_argument(
        '-m',
        '--metadata',
        required=True,
        type=str,
        help='TSV metadata file containing the "barcode" and "sample" columns to be used to rename'
    )
    parser.add_argument(
        '-b',
        '--barcode',
        required=True,
        type=str,
        help='Barcode number. Ex. "barcode05"'
    )
    return parser

def match_sample_barcode(barcode: str, metadata: str) -> Optional[str]:
    """
    Match the barcode given to the correct sample in the input metadata file
        If match found, that is returned
        Otherwise None
    """
    with open(metadata) as tsv:
        tsv_reader = csv.DictReader(tsv, delimiter='\t')
        for row in tsv_reader:
            # If the barcode is only 1 char its between 1-9 so add a 0 infront to match nanopore format
            if len(row['barcode']) == 1:
                row['barcode'] = f"0{row['barcode']}"
            # Return only exact match (should be 1 or weird input, maybe need a validation)
            if barcode == row['barcode']:
                return row['sample']
    return None

def main() -> None:
    '''Run the program'''
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Match barcode input to tsv
    barcode_match = re.search(BARCODE_REGEX, args.barcode)
    if barcode_match:
        barcode = str(barcode_match.group(1))
        sample = match_sample_barcode(barcode, args.metadata)

        if sample:
            fastq = Path(args.fastq)
            fastq.rename(f'{sample}.fastq')
            print(sample)

if __name__ == '__main__':
    main()
