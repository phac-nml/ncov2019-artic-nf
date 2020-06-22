#!/usr/bin/env python3

import argparse
import os
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
        '--fastq',
        required=False,
        default=False,
        help='Path to single fastq input file'
    )
    parser.add_argument(
        '--sample_info',
        required=True,
        help='Sampleinfo.csv file path'
    )
    parser.add_argument(
        '--prefix',
        required=True,
        help='Pass prefix string'
    )
    parser.add_argument(
        '--illumina',
        default=False,
        action='store_true',
        required=False,
        help='Specify if input is illumina data'
    )
    
    return parser


def parse_sample_csv(sample_csv, prefix, sample_dir, fastq):
    '''
    Take input sampleinfo.csv file and turn it into a usable format
    '''

    # Set up dataframe
    df_out = pd.DataFrame(columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])

    # Open file and populate df
    with open(sample_csv) as input_handle:
        for index, line in enumerate(input_handle):

            if index == 0:
                continue
            
            current_line_list = line.strip().split(',') # Order is [Name, run, barcode, project_id]

            if len(current_line_list) != 4:
                print('ERROR: Line {} of file {} is formatted incorrectly! Please address this by matching the format: [Name, Run, Barcode, Project_id]'.format(index + 1, sample_csv))
                quit()


            # Set correct file name and check on barcode formatting
            if int(current_line_list[2]) not in range(1,25):
                print('ERROR: Line {} of file {} does not contain an allowed barcode in range 1-24'.format(index + 1, sample_csv))
                quit()

            if len(current_line_list[2]) != 2: # Checking that barcode is 2 digits
                barcode = '0{}'.format(current_line_list[2])

                if len(barcode) != 2: # If its somehow not still...
                    print('ERROR: Line {} of file {} does not contain an allowed barcode in range 1-24'.format(index + 1, sample_csv))
                    quit()
            
            else:
                barcode = current_line_list[2]

            file_name = '{}_barcode{}.fastq'.format(prefix, barcode)
            new_file_name = '{}.fastq'.format(current_line_list[0])
            

            # Directory input with SampleList.csv output
            if sample_dir:

                # Check that File is found and rename it
                file_path = '{}/{}'.format(sample_dir, file_name)
                new_file_path = '{}/{}'.format(sample_dir, new_file_name)
                if os.path.exists(file_path):
                    subprocess.run('mv {} {}'.format(file_path, new_file_path), shell=True)

                else:
                    print('WARN: File {} not found in {}'.format(file_name, sample_dir))
                    continue


                # Set DataFrame for easy output csv
                df_out.at[index, 'Sample_Name'] = current_line_list[0] # Name from input sample info file
                df_out.at[index, 'Project_ID'] = current_line_list[3] # Project number from sample info file
                df_out.at[index, 'File_Forward'] = new_file_name

            # Fastq rename
            else:
                if str(fastq) == file_name:
                    if os.path.exists(file_name):
                        subprocess.run('mv {} {}'.format(fastq, new_file_name), shell=True)
                        return None # Exit program as it should match only once

                    else:    
                        print('ERROR: No file found matching {} in current directory. Exiting'.format(fastq))
                        quit()
    return df_out


def main():
    '''
    - Create Irida SampleList.csv with input directory
    - Rename singular fastq file with fastq input
    '''
    
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    sample_csv = args.sample_info
    prefix = args.prefix

    # Inputs
    if args.sample_dir and args.fastq:
        print('ERROR: Please specify only one input of either --sample_dir or --fastq')
        quit()

    elif args.sample_dir:
        sample_dir = args.sample_dir
        fastq = False
    
    elif args.fastq:
        fastq = args.fastq
        sample_dir = False

    else:
        print('ERROR: No input specified. Please specify either --sample_dir or --fastq')
        quit()


    # Start processing inputs
    if args.illumina: # Illumina can be paired or single, needs different support than nanopore
        print('Illumina processing is not supported at the moment')
        quit()


    else: # Nanopore data is single end

        df_out = parse_sample_csv(sample_csv, prefix, sample_dir, fastq)
        

        # Output
        if args.sample_dir:
            with open('{}/SampleList.csv'.format(sample_dir), 'w') as handle:
                handle.write('[Data]\n')
            
            df_out.to_csv('{}/SampleList.csv'.format(sample_dir), mode='a', header=True, index=False)



if __name__ == "__main__":
    main()
