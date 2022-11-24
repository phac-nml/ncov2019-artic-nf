#!/usr/bin/env python3

import argparse
import os
import pathlib
import re
import subprocess
import pandas as pd

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sample_dir',
        required=False,
        default=False,
        help='Path to barcoded fastq sample directory to output SampleList.csv for the directory'
    )
    parser.add_argument(
        '--sample_info',
        required=True,
        help='File path to tab separated sample sheet with the following columns ["sample", "run", "barcode", "project_id", "ct", "date"]'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--fastq', action='store_true', help="Create SampleList for nanopore fastq files")
    group.add_argument('--fasta', action='store_true', help="Create SampleList for consensus fasta files")
    return parser

def create_sample_file_df(sample_tsv, sample_dir, file_type='', file_list=[]):
    '''
    Take input sample sheet, directory, and file type to create a df of values needed for IRIDA Uploads
        Rename the files to include the run name for better tracking
    '''
    # Read in input TSV file
    df = pd.read_csv(sample_tsv, sep='\t')
    df['sample'] = df['sample'].astype("string")

    # Get all of the fasta files in the input folder to a list to check against our sample_tsv df
    if file_type == 'fastq':
        FILE_WANTED_REGEX = re.compile(r'^.+\.fastq$')
    else:
        FILE_WANTED_REGEX = re.compile(r'^.+\.fasta$')
    for file_found in os.listdir(sample_dir):
        if FILE_WANTED_REGEX.match(file_found):
            file_list.append(file_found)

    # Generate Out DF from values in input table that have matching files found
    row_list = []
    for sample, proj_id, run in zip(df['sample'], df['project_id'], df['run']):
        # Check if sample string a is a substring of any files and then add it to df_out if so
        # Check so that only wanted samples are uploaded
        for file_name in file_list:
            if sample in file_name:
                # Rename file keeping the same extension
                ext = ''.join(pathlib.Path(file_name).suffixes)
                new_file_name = '{}_{}{}'.format(sample, run, ext)
                subprocess.run('mv {}/{} {}/{}'.format(sample_dir, file_name, sample_dir, new_file_name), shell=True)
                row_dict = {
                    'Sample_Name': sample,
                    'Project_ID': proj_id, 
                    'File_Forward': new_file_name, 
                    'File_Reverse': ''
                }
                row_list.append(row_dict)
    df_out = pd.DataFrame(row_list, columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])
    return df_out


def main():
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Create df_out using directory files and the sample TSV file based on the input file type wanted
    if args.fastq:
        df_out = create_sample_file_df(args.sample_info, args.sample_dir, file_type='fastq')
    else:
        df_out = create_sample_file_df(args.sample_info, args.sample_dir, file_type='fasta')

    # Output
    with open('{}/SampleList.csv'.format(args.sample_dir), 'w') as handle:
        handle.write('[Data]\n')
    
    df_out.to_csv('{}/SampleList.csv'.format(args.sample_dir), mode='a', header=True, index=False)

if __name__ == "__main__":
    main()
