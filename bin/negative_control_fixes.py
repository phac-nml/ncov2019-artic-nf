#!/usr/bin/env python3
"""Finalize pipeline CSV file by correcting blank values and evaluating controls"""

import argparse
import pandas as pd
import numpy as np

def init_parser() -> argparse.ArgumentParser:
    """
    Specify command line arguments
    Returns command line parser with inputs
    """
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
        '--count_failure_tsv',
        required=False,
        help='Path to tsv containing samples failing read count filter'
    )
    parser.add_argument(
        '--filter_failure_tsv',
        required=False,
        help='Path to tsv containing samples failing read length filter (artic Guppyplex too few reads)'
    )
    return parser

def main() -> None:
    '''
    One function to:
        - Fix up the negative control columns
        - Add in failing samples removed by pipeline for failing read count checks
        - Order final CSV columns that are available
    '''
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Read in qc csv file
    df = pd.read_csv(args.qc_csv, sep=',', dtype=object)

    # Set numeric columns so that they can be filled with 0 instead of NA. These columns ALWAYS part of output
    numeric_columns = [
        'num_aligned_reads',
        'num_consensus_n',
        'num_consensus_snvs',
        'num_consensus_iupac',
        'num_variants_snvs',
        'num_variants_indel',
        'num_variants_indel_triplet',
        'mean_sequencing_depth',
        'median_sequencing_depth',
        'genome_completeness'
    ]

    ### Negative Control Fixes ###
    # Set the negative control headers which are static based on ncov-tools
    negative_columns = ['qc', 'genome_covered_bases', 'genome_total_bases', 'genome_covered_fraction', 'amplicons_detected']
    # Place data here as column: [statement(s)] to generate the fill for the columns
    replace_dict = {}

    # Fix negative control columns to contain the negative control values
    for column in negative_columns:
        # Get integer locations of the non-null data in the negative columns
        replace_dict[column] = []
        negative_data_index = np.where(df[column].notnull())[0].tolist()
        if negative_data_index == []:
            continue
        # Combine all negative control column values
        for spot in negative_data_index:
            # Column is always called sample based on ncov-tools output, if it changes have to change this
            print(spot)
            sample = df['sample'][spot]
            replace_dict[column].append('{}:{}'.format(sample, df[column][spot]))

    # Add the combined negative columns to their column separated by a `-`
    for key in replace_dict.keys():
        if replace_dict[key] == []:
            continue
        df[key] = '-'.join(replace_dict[key])
    ### End Negative Control Fixes ###

    ### Adding in failed sample tsv files ###
    # Setup to concat DFs if we have any
    frames = [df]
    df_columns = df.columns.values

    # Append on the other failed samples for tracking
    # For the nanopore pipeline fail tracking
    # Rename to match what comes out of qc.py
    rename_columns = {
        'run': 'run_identifier',
        'ct': 'qpcr_ct',
        'date': 'collection_date'
    }
    if args.count_failure_tsv:
        count_failure_df = pd.read_csv(args.count_failure_tsv, sep='\t', dtype=object)
        count_failure_df.rename(columns={key: val for key, val in rename_columns.items() if val in df_columns}, inplace=True)
        frames.append(count_failure_df)
    if args.filter_failure_tsv:
        filter_failure_df = pd.read_csv(args.filter_failure_tsv, sep='\t', dtype=object)
        filter_failure_df.rename(columns={key: val for key, val in rename_columns.items() if val in df_columns}, inplace=True)
        frames.append(filter_failure_df)
    ### End appending failing samples ###

    ### Create final concated df and then output ###
    final_df = pd.concat(frames)
    # Fill blank numeric columns with 0 and then other columns with NA
    final_df[numeric_columns] = final_df[numeric_columns].fillna(value=0)
    final_df.fillna('NA', inplace=True)
    # Rearrange final output columns to output similar order each time
    cols = list(final_df.columns)
    key_cols = ['sample', 'run_identifier', 'barcode', 'project_id', 'num_aligned_reads', 'num_consensus_n', 'lineage', 'scorpio_call', 'variants', 'protein_variants']
    # Some key_cols are not mandatory so remove those while keeping the order of the csv file
    key_cols = [x for x in key_cols if x in cols]
    extra_cols = [x for x in cols if (x not in key_cols) and (x not in negative_columns)]
    final_df = final_df[key_cols+extra_cols+negative_columns]

    # Drop any additional columns
    final_df.drop(['nextflow_qc_pass'], axis=1, inplace=True)

    # Final Output
    final_df.to_csv('{}.qc.csv'.format(args.output_prefix), index=False)

if __name__ == "__main__":
    main()
