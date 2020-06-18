#!/usr/bin/env python3

import argparse
import subprocess
import os

import pandas as pd

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fastq',
        required=True,
        help='Path to fastq file'
    )
    parser.add_argument(
        '--sample_csv',
        required=True,
        help='Path to csv sample list'
    )
    parser.add_argument(
        '--prefix',
        required=True,
        help='Prefix'
    )
    
    return parser


def parse_sample_csv(sample_csv, prefix, fastq):
    '''
    Take input sampleinfo.csv file and turn it into a usable format
    '''

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

            if str(fastq) == file_name:
                if os.path.exists(file_name):
                    subprocess.run('mv {} {}'.format(fastq, new_file_name), shell=True)
                    return None # Exit program as it should match only once

                else:    
                    print('ERROR: No file found matching {} in current directory. Exiting'.format(fastq))
                    quit()

        # Nothing found, sample file may be off or GuppyPlex has done something weird
        print('WARNING: No files found that match to the contents of {}'.format(sample_csv))


def main():
    '''
    Main script
    '''
    
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    sample_csv = args.sample_csv
    prefix = args.prefix
    fastq = args.fastq


    # Process
    parse_sample_csv(sample_csv, prefix, fastq)

if __name__ == "__main__":
    main()
