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
        help='File path to tab separated sampleinfo file with format [Name, Run, Barcode, Project_id, ct]'
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='File path to wanted output directory'
    )
    
    return parser

def generate_samplelist(sample_tsv, directory, output_dir):

    df_out = pd.DataFrame(columns=['Sample_Name', 'Project_ID', 'File_Forward', 'File_Reverse'])
    input_format = ['sample', 'run', 'barcode', 'project_id', 'ct']

    with open(sample_tsv) as input_handle:

        for index, line in enumerate(input_handle):

            current_line_list = line.strip('\n').split('\t') # Order is [sample, run, barcode, project_id, ct]


            # Error checking for formatting
            if index == 0:
                if len(current_line_list) != 5:
                    if len(current_line_list) <= 4:
                        print('ERROR: Header formated incorrectly. Please address this by matching the format {} for the first columns'.format(input_format))
                        quit()

                    # Trim other rows to allow pipeline to continue if it passes other checks
                    # Allows end users to have additional information if they match the first 5 columns correctly
                    del current_line_list[5:len(current_line_list)]

                if current_line_list != input_format:
                    print('ERROR: Header formated incorrectly. Please address this by matching the format {} for the first columns'.format(input_format))
                    quit()

                else:
                    continue

            # For all rows, check if they are 5 
            if len(current_line_list) == 5:
                pass

            else:
                if len(current_line_list) <= 4:
                    print('ERROR: Line {} is formatted incorrectly! Please address this by matching the format {} for the first columns'.format(index + 1, input_format))
                    quit()
                
                # Same as above, cut off extra data if there is any. No need to fail if format is correct in the header.
                del current_line_list[5:len(current_line_list)]


            # Set correct file name and check on barcode formatting
            if len(current_line_list[2]) < 2: # Checking that barcode is 2 digits
                barcode = '0{}'.format(current_line_list[2])

            else:
                barcode = current_line_list[2]

            # With barcode, check if that fast5 directory exists and if so zip all files together and move them
            # to the output directory
            subprocess.run('mkdir -p {}'.format(output_dir), shell=True)

            fast5_location = '{}/barcode{}/'.format(directory, barcode)
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
                        new_file_name = '{}{}.fast5'.format(current_line_list[0], file_number)
                        subprocess.run('ln -s {}{} {}'.format(fast5_location, file_name, new_file_name), shell=True)
                        sed_cmd_list.append("sed -i -e 's/{}/{}/' filename_mapping.txt".format(file_name, new_file_name))

                # Actually run it
                if not mapping_found:
                    print('Skipping fast5 uploading. Missing filename_mapping.txt in barcode{} which may indicate fast5 not dehosted for human reads.'.format(barcode))
                    return df_out

                else:
                    for cmd in sed_cmd_list:
                        subprocess.run(cmd, shell=True)

                    subprocess.run('echo "read_id\tfilename" > {}_sequencing_summary.txt && cat  filename_mapping.txt >> {}_sequencing_summary.txt'.format(current_line_list[0], current_line_list[0]), shell=True)
                    subprocess.run('tar -cvh --use-compress-program=pigz -f {}/{}.tar.gz *.fast5 *_sequencing_summary.txt && rm *.fast5 *.txt'.format(output_dir, current_line_list[0]), shell=True)


                df_out.at[index, 'Sample_Name'] = current_line_list[0] # Name from input sample info file
                df_out.at[index, 'Project_ID'] = current_line_list[3] # Project number from sample info file
                df_out.at[index, 'File_Forward'] = '{}.tar.gz'.format(current_line_list[0]) # File name created by subprocess

            else:
                print('WARNING: No files found for barcode{}'.format(barcode))
                pass

    return df_out


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    sample_tsv = args.sample_info
    directory = args.sample_dir
    output_dir = args.output_dir


    # Parse samplesheet to make dataframe
    df_out = generate_samplelist(sample_tsv, directory, output_dir)

    # Create output SampleList
    with open('{}/SampleList.csv'.format(output_dir), 'w') as handle:
        handle.write('[Data]\n')
        
    df_out.to_csv('{}/SampleList.csv'.format(output_dir), mode='a', header=True, index=False)

if __name__ == "__main__":
    main()
