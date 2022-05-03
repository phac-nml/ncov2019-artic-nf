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
    parser.add_argument(
        '--read_tsv',
        required=False,
        help='Path to tsv containing samples failing read count filter'
    )
    parser.add_argument(
        '--mapping_tsv',
        required=False,
        help='Path to tsv containing samples failing read mapping count filter'
    )

    return parser


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Read in qc csv file
    df = pd.read_csv(args.qc_csv, sep=',', dtype=object)

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

    # Setup to concat DFs if we have any
    frames = [df]
    df_columns = df.columns.values
    # Rename to match what comes out of qc.py
    rename_columns = {
            'run': 'run_identifier',
            'ct': 'qpcr_ct',
            'date': 'collection_date'
        }

    # Append on the other failed samples for tracking
    if args.read_tsv:
        read_df = pd.read_csv(args.read_tsv, sep='\t', dtype=object)
        read_df.rename(columns={key: val for key, val in rename_columns.items() if val in df_columns}, inplace=True)
        frames.append(read_df)
    if args.mapping_tsv:
        mapping_df = pd.read_csv(args.mapping_tsv, sep='\t', dtype=object)
        mapping_df.rename(columns={key: val for key, val in rename_columns.items() if val in df_columns}, inplace=True)
        frames.append(mapping_df)

    # Create final concated df and then output
    final_df = pd.concat(frames)
    # Fill blank columns
    final_df.fillna('NA', inplace=True)
    # Rearrange final output columns
    cols = list(df.columns)
    key_cols = ['sample', 'run_identifier', 'barcode', 'project_id', 'num_aligned_reads', 'num_consensus_n', 'lineage', 'variants', 'protein_variants']
    extra_cols = [x for x in cols if x not in key_cols]
    final_df = final_df[key_cols+extra_cols]
    final_df.to_csv('{}.qc.csv'.format(args.output_prefix), index=False)

if __name__ == "__main__":
    main()
