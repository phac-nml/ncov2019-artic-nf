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

    with open(sample_tsv) as input_handle:

        for index, line in enumerate(input_handle):

            current_line_list = line.strip('\n').split('\t') # Order is [sample, run, barcode, project_id, ct]

            # Error Checking Formatting
            if index == 0:
                if str(current_line_list[0]) != 'sample':
                    print('ERROR: First column of the header line must be called "sample". Exiting')
                    quit() 
                else:
                    continue
            
            if len(current_line_list) != 5:
                print('ERROR: Line {} of file {} is formatted incorrectly! Please address this by matching the format: [sample, run, barcode, project_id, ct]'.format(index + 1, sample_tsv))
                quit()


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
