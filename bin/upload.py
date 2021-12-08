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

def send_metadata(api_instance, metadata_csv):
    '''
    PURPOSE:
        Send metadata from qc.csv to IRIDA for each sample

    INPUTS:
        - API_INSTANCE --> Irida API instance from generate_api_instance

        - METADATA_CSV --> CSV file that contains all of the metadata along with the sample name and project id
            in the first and second columns respectively

    RETURNS:
        TRACKING_DICT_LIST --> List of dictionarys containing the status of each samples upload
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
        if api_instance.sample_exists(sample_name=sample_name, project_id=project_id):
            irida_metadata = api_instance.get_metadata(sample_name, project_id)
            # If no metadata yet in IRIDA pass and just upload
            if irida_metadata == {}:
                pass
            # If there is metadata, check that the new data is better than the old data
            else:
                irida_n_count = int(irida_metadata['num_consensus_n']['value'])
                if irida_n_count <= metadata_dict['num_consensus_n']:
                    print('Skipped sample {} metadata upload as IRIDA N Count {} <= {} New Sample N Count'.format(sample_name, irida_n_count, metadata_dict['num_consensus_n']))
                    tracking_dict_list.append(_create_track_dict(sample_name, project_id, False, 'IRIDA N Count {} <= {} New Sample N Count'.format(irida_n_count, metadata_dict['num_consensus_n'])))
                    continue
        # If sample does not exist, make it exist
        else:
            irida_sample = model.Sample(sample_name=sample_name)
            api_instance.send_sample(sample=irida_sample, project_id=project_id)

        upload_metadata = model.Metadata(metadata=metadata_dict, project_id=project_id, sample_name=sample_name)
        status = api_instance.send_metadata(upload_metadata, upload_metadata.project_id, upload_metadata.sample_name)
        print('Uploaded {} metadata to {}'.format(sample_name, project_id))
        tracking_dict_list.append(_create_track_dict(sample_name, project_id, True, ''))
            
    return tracking_dict_list

def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    irida_api = generate_api_instance(args.config)
    tracking_dict_list = send_metadata(irida_api, args.metadata_csv)

    df = pd.DataFrame.from_dict(tracking_dict_list)
    df.to_csv('metadata_upload_status.csv', index=False)

if __name__ == "__main__":
    main()
