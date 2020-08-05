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
        help='File path to tab separated sampleinfo file with format [Name, Run, Barcode, Project_id, ct]'
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


def parse_sample_tsv(sample_tsv, prefix, sample_dir, fastq):
    '''
    Take input sampleinfo file and turn it into a usable format
    '''

    # Set up dataframe
    df_out = pd.DataFrame(columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])
    input_format = ['sample', 'run', 'barcode', 'project_id', 'ct']

    # Open file and populate df
    with open(sample_tsv) as input_handle:

        file_info_list = []
        file_names_set = set()

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


            # Set correct file name and check on barcode formatting
            if int(current_line_list[2]) not in range(1,25):
                print('ERROR: Line {} of file {} does not contain an allowed barcode in range 1-24'.format(index + 1, sample_tsv))
                quit()

            if len(current_line_list[2]) != 2: # Checking that barcode is 2 digits
                barcode = '0{}'.format(current_line_list[2])

                if len(barcode) != 2: # If its somehow not still...
                    print('ERROR: Line {} of file {} does not contain an allowed barcode in range 1-24'.format(index + 1, sample_tsv))
                    quit()
            
            else:
                barcode = current_line_list[2]

            file_name = '{}_barcode{}.fastq'.format(prefix, barcode)
            current_line_list[2] = barcode

            file_info_list.append(current_line_list)
            file_names_set.add(file_name)


    for index, sample_list in enumerate(file_info_list): # Order is [sample, run, barcode, project_id, ct]
        
        file_name = '{}_barcode{}.fastq'.format(prefix, sample_list[2])
        new_file_name = '{}.fastq'.format(sample_list[0])

        # For directory input we create SampleList.csv as the output and must make it
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
            df_out.at[index, 'Sample_Name'] = sample_list[0] # Name from input sample info file
            df_out.at[index, 'Project_ID'] = sample_list[3] # Project number from sample info file
            df_out.at[index, 'File_Forward'] = new_file_name


        # Else we just rename the singular input fastq file
        else:
            if str(fastq) not in file_names_set:
                if os.path.exists(fastq):
                    subprocess.run('mv {} extra_{}'.format(fastq, fastq), shell=True)
                    return None # Exit program as it should match only once
                
                else:    
                    print('ERROR: No file found matching {} in current directory. Exiting'.format(fastq))
                    quit()

            elif str(fastq) == file_name:
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

    sample_tsv = args.sample_info
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

        df_out = parse_sample_tsv(sample_tsv, prefix, sample_dir, fastq)
        

        # Output
        if args.sample_dir:
            with open('{}/SampleList.csv'.format(sample_dir), 'w') as handle:
                handle.write('[Data]\n')
            
            df_out.to_csv('{}/SampleList.csv'.format(sample_dir), mode='a', header=True, index=False)



if __name__ == "__main__":
    main()
