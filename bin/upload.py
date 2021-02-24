#!/usr/bin/env python3

import argparse
import configparser
import csv
import re

from iridauploader import api, model


def init_parser():
    '''
    Specify command line arguments for automatic upload to irida

    Returns command line parser with inputs
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--config',
        required=True,
        help='Irida Uploader config file to parse'
    )
    parser.add_argument(
        '-m',
        '--metadata_csv',
        required=True,
        help='Metadata qc csv for upload containing matching ids to SampleList.csv in the FIRST column and a project id in the SECOND'
    )
    return parser


def generate_api_instance(config_in):
    '''
    PURPOSE:
        Generate api instance for IRIDA platform based on input keys found in csv file

    INPUTS:
        - KEY_CONFIG --> IRIDA Uploader config file

    RETURNS:
        - IRIDA_API --> API class instance based on the input keys given in the config
    '''

    config = configparser.ConfigParser()
    config.read(config_in)

    settings = config['Settings']

    irida_api = api.ApiCalls(settings['client_id'], settings['client_secret'], settings['base_url'], settings['username'], settings['password'])

    return irida_api
                

def send_metadata(api_instance, metadata_csv):
    '''
    PURPOSE:
        Send metadata from qc.csv to IRIDA for each sample

    INPUTS:
        - API_INSTANCE --> Irida API instance from generate_api_instance

        - METADATA_CSV --> CSV file that contains all of the metadata along with the sample name and project id
            in the first and second columns respectively
    '''

    with open(metadata_csv, 'r') as input_handle:
        reader = csv.reader(input_handle)

        for index, row in enumerate(reader):
            
            if index == 0:
                header = row
                continue

            else:
                # Create dictionary of each row for creation of dictionary to upload to irida
                metadata = {}
                passing = True

            for i in range(len(row)):

                if i == 0 and re.search('sample', header[0]):
                    # Get sample name from row 1
                    sample_name = row[i]
                
                # Get the project ID from the correct header no matter where it is
                elif re.search('project_id', header[i]):

                    if row[i] == 'Unknown' or row[i] == 'NA':
                        passing = False
                        break

                    # Get project id from row 2
                    project_id = row[i]

                else:
                    # Put metadata into metadata dictionary for upload
                    metadata[header[i]] = row[i]
            if passing:
                # Check that sample exists and make it if not
                if api_instance.sample_exists(sample_name=sample_name, project_id=project_id):
                    pass
                else:
                    irida_sample = model.Sample(sample_name=sample_name)
                    api_instance.send_sample(sample=irida_sample, project_id=project_id)

                upload_metadata = model.Metadata(metadata=metadata, project_id=project_id, sample_name=sample_name)
                status = api_instance.send_metadata(upload_metadata, upload_metadata.project_id, upload_metadata.sample_name )

                print(status, '\n')

            else:
                print('Unknown sample data for {}, moving to next sample'.format(sample_name))


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    key_config = args.config
    metadata_csv = args.metadata_csv

    irida_api = generate_api_instance(key_config)
    send_metadata(irida_api, metadata_csv)


if __name__ == "__main__":
    main()
