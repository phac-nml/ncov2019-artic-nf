#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import numpy as np

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--qc_csv',
        required=True,
        help='Path to qc_csv to check for negative control status'
    )

    parser.add_argument(
        '--output_prefix',
        required=True,
        help='Output file prefix'
    )

    return parser


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Read in qc csv file
    df = pd.read_csv(args.qc_csv, sep=',')

    # Set the negative control headers which are static based on ncov-tools
    negative_columns = ['qc', 'genome_covered_bases', 'genome_total_bases', 'genome_covered_fraction', 'amplicons_detected']

    # Place data here as column: [statement(s)] to generate the fill for the columns
    replace_dict = {}

    for column in negative_columns:
        # Get integer locations of the non-null data in the negative columns
        replace_dict[column] = []
        negative_data_index = np.where(df[column].notnull())[0].tolist()
        
        if negative_data_index == []:
            continue

        for spot in negative_data_index:
            # Column is always called sample based on ncov-tools output, if it changes have to change this
            print(spot)
            sample = df['sample'][spot]
            replace_dict[column].append('{}:{}'.format(sample, df[column][spot]))

    for key in replace_dict.keys():
        if replace_dict[key] == []:
            continue
        df[key] = '-'.join(replace_dict[key])

    df.fillna('NA', inplace=True)
    df.to_csv('{}.qc.csv'.format(args.output_prefix), index=False)


if __name__ == "__main__":
    main()
