#!/usr/bin/env python3

import argparse
import re
import pandas as pd

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--barcode_list',
        required=True,
        help='Path to list of barcodes to determine their sample name'
    )

    parser.add_argument(
        '--samplesheet',
        required=True,
        help='Path so samplesheet linking barcode number to sample name'
    )
    return parser

def main():
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    df = pd.read_csv(args.samplesheet, sep="\t")

    # Add columns
    print("sample,qc_pass")

    with open(args.barcode_list) as input_handle:
        for line in input_handle:
            # Remove "barcode" from name to just get the number
            barcode_number = re.sub("\D", "", line.strip())

            # Use barcode number to find the sample
            try:
                sample = str(df[df['barcode'] == int(barcode_number)]['sample'].item())
            except ValueError:
                sample = "extra_nml_barcode{}".format(barcode_number)
            # Just print sample to send to output file
            print("{},TOO_FEW_READS".format(sample))

if __name__ == "__main__":
    main()
