#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sample_dir',
        required=True,
        help='Path to directory with consensus fasta files'
    )
    parser.add_argument(
        '--sample_info',
        required=True,
        help='File path to tab separated sampleinfo file with format [Name, Run, Barcode, Project_id, ct]'
    )
    
    return parser

def generate_samplelist(sample_tsv, directory):

    df_out = pd.DataFrame(columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])
    input_format = ['sample', 'run', 'barcode', 'project_id', 'ct']

    with open(sample_tsv) as input_handle:

        for index, line in enumerate(input_handle):

            current_line_list = line.strip('\n').split('\t') # Order is [sample, run, barcode, project_id, ct]


            # Error checking for formatting
            if index == 0:
                if len(current_line_list) != 5:
                    if len(current_line_list) <= 4:
                        print('ERROR: Header formated incorrectly. Please address this by matching the format: {}'.format(input_format))
                        quit()

                    del current_line_list[5:len(current_line_list)]

                if current_line_list != input_format:
                    print('ERROR: Header formated incorrectly. Please address this by matching the format: {}'.format(input_format))
                    quit()

                else:
                    continue

            # For all rows, check if they are 5 
            if len(current_line_list) == 5:
                pass

            else:
                if len(current_line_list) <= 4:
                    print('ERROR: Line {} is formatted incorrectly! Please address this by matching the format: {}'.format(index + 1, input_format))
                    quit()

                del current_line_list[5:len(current_line_list)]


            # File Checking
            file_path = '{}/{}.consensus.fasta'.format(directory, current_line_list[0])

            if os.path.exists(file_path):

                # Set DataFrame for easy output csv
                df_out.at[index, 'Sample_Name'] = current_line_list[0]
                df_out.at[index, 'Project_ID'] = current_line_list[3]
                df_out.at[index, 'File_Forward'] = '{}.consensus.fasta'.format(current_line_list[0])

            else:
                print('WARNING: File {}.consensus.fasta not found in {} directory.'.format(current_line_list[0], directory))

    return df_out


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    sample_tsv = args.sample_info
    directory = args.sample_dir


    # Parse samplesheet to make dataframe
    df_out = generate_samplelist(sample_tsv, directory)

    # Create output SampleList
    with open('{}/SampleList.csv'.format(directory), 'w') as handle:
        handle.write('[Data]\n')
        
    df_out.to_csv('{}/SampleList.csv'.format(directory), mode='a', header=True, index=False)

if __name__ == "__main__":
    main()
