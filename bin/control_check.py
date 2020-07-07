#!/usr/bin/env python3

import argparse
import pandas as pd

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--qc_csv',
        required=True,
        help='Nextflow *.qc.csv file'
    )
    parser.add_argument(
        '--control_names',
        required=True,
        help='Names of the negative contol samples to check against'
    )
    
    return parser


def check_controls(qc_csv, negative_names):

    negative_name_list = negative_names.replace(' ', '').split(',')
    qc_df =  pd.read_csv(qc_csv)
    spot = qc_df.shape[1] -1 # Want second last position to add control check

    for index, sample_name in enumerate(qc_df['sample'].tolist()):
        if sample_name in negative_name_list and qc_df['qc_pass'][index] == True:
        
            qc_df.insert(spot, 'control', 'FAIL')
            return qc_df

    qc_df.insert(spot, 'control', 'Pass')
    return qc_df
    

def main():

    parser = init_parser()
    args = parser.parse_args()

    qc_csv = args.qc_csv
    negative_names = args.control_names

    qc_df = check_controls(qc_csv, negative_names)

    qc_df.to_csv('out.csv', header=True, index=False)


if __name__ == "__main__":
    main()