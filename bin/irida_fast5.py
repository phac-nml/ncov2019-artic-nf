#!/usr/bin/env python3

import argparse
import os
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
        required=True,
        help='Path to directory with input fast5 files. Format should be dir/barcode##/samples.fast5'
    )
    parser.add_argument(
        '--sample_info',
        required=True,
        help='File path to tab separated sample sheet with the following columns ["sample", "run", "barcode", "project_id", "ct", "date"]'
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='File path to wanted output directory'
    )
    return parser

def create_sample_archive_df(sample_tsv, directory, output_dir, row_list=[]):
    '''
    Take input sample sheet, input directory, and output directory to create archives of fast5 files and the df of values needed for IRIDA Uploads
    '''
    df = pd.read_csv(sample_tsv, sep='\t')
    for sample, proj_id, barcode, run in zip(df['sample'], df['project_id'], df['barcode'], df['run']):
        # Check that barcode is at least 2 digits long
        if len(str(barcode)) < 2:
            barcode_str = '0{}'.format(barcode)
        else:
            barcode_str = str(barcode)

        # With barcode, check if that fast5 directory exists and if so zip all files together and move them to the output directory
        fast5_location = '{}/barcode{}/'.format(directory, barcode_str)
        if os.path.exists(fast5_location):
            file_name_list = os.listdir(fast5_location)
            mapping_found = False
            sed_cmd_list = []

            for file_name in file_name_list:
                
                # Need filename_mapping.txt or its not dehosted
                if re.search(r'filename_mapping', file_name):
                    subprocess.run('cp {}filename_mapping.txt filename_mapping.txt'.format(fast5_location), shell=True)
                    mapping_found = True
                    continue

                # Only want to rename fast5 files
                if re.search(r'\.fast5', file_name):
                    # Set new file name with the correct number
                    # All files end with _## or _# and we need to keep that
                    file_number = re.findall(r'_\d+', os.path.splitext(file_name)[0])[-1]
                    new_file_name = '{}{}.fast5'.format(sample, file_number)
                    subprocess.run('ln -s {}{} {}'.format(fast5_location, file_name, new_file_name), shell=True)
                    sed_cmd_list.append("sed -i -e 's/{}/{}/' filename_mapping.txt".format(file_name, new_file_name))

            # Actually run it
            # This part should never trigger in pipeline, should remove/revisit
            if not mapping_found:
                print('ERROR: Missing filename_mapping.txt in barcode{} which may indicate fast5 not dehosted for human reads.'.format(barcode))
                exit(1)
            else:
                for cmd in sed_cmd_list:
                    subprocess.run(cmd, shell=True)

                subprocess.run('mv filename_mapping.txt {}_filename_mapping.txt'.format(sample), shell=True)
                subprocess.run('tar -cvh --use-compress-program=pigz -f {}/{}_{}.tar.gz *.fast5 *_filename_mapping.txt && rm *.fast5 *.txt'.format(output_dir, sample, run), shell=True)

            # To output
            row_dict = {
                'Sample_Name': sample,
                'Project_ID': proj_id, 
                'File_Forward': '{}_{}.tar.gz'.format(sample, run), 
                'File_Reverse': ''
            }
            row_list.append(row_dict)
        else:
            print('WARNING: No files found for barcode{}'.format(barcode))
            pass

    df_out = pd.DataFrame(row_list, columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])
    return df_out


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    sample_tsv = args.sample_info
    directory = args.sample_dir
    output_dir = args.output_dir

    # Make output directory
    subprocess.run('mkdir -p {}'.format(output_dir), shell=True)
    # Parse samplesheet to make dataframe along with zipping up fast5 files
    df_out = create_sample_archive_df(sample_tsv, directory, output_dir)

    # Create output SampleList
    with open('{}/SampleList.csv'.format(output_dir), 'w') as handle:
        handle.write('[Data]\n')
        
    df_out.to_csv('{}/SampleList.csv'.format(output_dir), mode='a', header=True, index=False)

if __name__ == "__main__":
    main()
