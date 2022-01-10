#!/usr/bin/env python3

import argparse
import configparser
import re
import pandas as pd

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
    parser.add_argument(
        '--no_sample_creation',
        required=False,
        default=False,
        action='store_true',
        help='Turns off sample creation for sample names not found in IRIDA'
    )
    return parser

def generate_api_instance(config_in):
    '''
    PURPOSE:
        Generate api instance for IRIDA platform based on input keys found in csv file
    INPUTS:
        - KEY_CONFIG: IRIDA Uploader config file to be used to connect to IRIDA
    RETURNS:
        - IRIDA_API: API class instance based on the input keys given in the config
    '''
    config = configparser.ConfigParser()
    config.read(config_in)
    settings = config['Settings']
    irida_api = api.ApiCalls(settings['client_id'], settings['client_secret'], settings['base_url'], settings['username'], settings['password'])
    return irida_api

def _create_track_dict(sample, project_id, uploaded, status):
    """
    Create and Return dictionary for tracking status of uploaded metadata samples
    """
    d = {
        'sample': sample,
        'project_id': project_id,
        'uploaded': uploaded,
        'status': status
    }
    return d

def send_metadata(api_instance, metadata_csv, no_sample_creation):
    '''
    PURPOSE:
        Send metadata for each sample in the input CSV file given a valid project ID
    INPUTS:
        - API_INSTANCE: Irida API instance from generate_api_instance
        - METADATA_CSV: CSV file that contains all of the metadata along with the sample name and project id
            in the first and second columns respectively
        - NO_SAMPLE_CREATION: Boolean value, if true, samples not already in IRIDA are skipped instead of created
    RETURNS:
        TRACKING_DICT_LIST: List of dictionarys containing the status of each samples upload
    '''
    df_in_dict = pd.read_csv(metadata_csv).fillna('NA').to_dict(orient='records')
    tracking_dict_list = [] # List to append dicts that track the status of samples for data uploads
    for metadata_dict in df_in_dict:
        sample_name = metadata_dict.pop('sample')
    
        # Check if we have a project ID that is an integer for IRIDA
        try:
            project_id = int(metadata_dict.pop('project_id'))
        except ValueError:
            print('Unknown Project ID {} for sample {}, moving to next sample'.format(project_id, sample_name))
            tracking_dict_list.append(_create_track_dict(sample_name, project_id, False, 'Unknown Project ID'.format()))
            continue

        # Check that sample exists and if the new values are better than previous
        sample_id = api_instance.get_sample_id(sample_name=sample_name, project_id=project_id)
        if sample_id:
            irida_metadata = api_instance.get_metadata(sample_id)
            # If there is metadata, check that the new data is better than the old data first with N count, and then with aligned read count
            if irida_metadata != {}:
                irida_n_count = int(irida_metadata['num_consensus_n']['value'])
                if irida_n_count < metadata_dict['num_consensus_n']:
                    print('Skipped sample {} metadata upload as IRIDA N Count {} < {} New Sample N Count'.format(sample_name, irida_n_count, metadata_dict['num_consensus_n']))
                    tracking_dict_list.append(_create_track_dict(sample_name, project_id, False, 'IRIDA N Count {} < {} New Sample N Count'.format(irida_n_count, metadata_dict['num_consensus_n'])))
                    continue
                irida_aligned_reads = int(irida_metadata['num_aligned_reads']['value'])
                if irida_aligned_reads >= metadata_dict['num_aligned_reads']:
                    print('Skipped sample {} metadata upload as IRIDA aligned read count {} >= {} New Sample aligned read count'.format(sample_name, irida_aligned_reads, metadata_dict['num_aligned_reads']))
                    tracking_dict_list.append(_create_track_dict(sample_name, project_id, False, 'IRIDA aligned read count {} >= {} New Sample aligned read count'.format(irida_aligned_reads, metadata_dict['num_aligned_reads'])))
                    continue
        # If given, skip samples that are not in IRIDA already
        elif no_sample_creation:
            print('Skipped sample {} metadata as it is not in IRIDA and --no_sample_creation arg passed'.format(sample_name))
            tracking_dict_list.append(_create_track_dict(sample_name, project_id, False, 'Not in IRIDA and --no_sample_creation arg passed'))
            continue
        # If sample does not already exist, create it and grab its sample_id
        else:
            irida_sample = model.Sample(sample_name=sample_name)
            response = api_instance.send_sample(sample=irida_sample, project_id=project_id)
            sample_id = response['resource']['identifier']

        # Actually Upload data
        upload_metadata = model.Metadata(metadata=metadata_dict, project_id=project_id, sample_name=sample_name)
        status = api_instance.send_metadata(upload_metadata, sample_id)
        print('Uploaded {} metadata to {}'.format(sample_name, project_id))
        tracking_dict_list.append(_create_track_dict(sample_name, project_id, True, ''))

    return tracking_dict_list

def main():
    
    # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Connect and upload metadata
    irida_api = generate_api_instance(args.config)
    tracking_dict_list = send_metadata(irida_api, args.metadata_csv, args.no_sample_creation)

    # Output tracking df
    df = pd.DataFrame.from_dict(tracking_dict_list)
    df.to_csv('metadata_upload_status.csv', index=False)

if __name__ == "__main__":
    main()
