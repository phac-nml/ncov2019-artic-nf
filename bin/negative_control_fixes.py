#!/usr/bin/env python3

import argparse
import re
import pandas as pd
import numpy as np

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--qc_csv',
        required=True,
        help='Path to qc_csv to check for negative control status'
    )

    return parser


def main():
    
     # Init Parser and set arguments
    parser = init_parser()
    args = parser.parse_args()

    # Read in qc csv file
    df = pd.read_csv(args.qc_csv, sep=',')

    # Grab the negative control headers
    # They have negative in them as it was set earlier
    columns = df.columns
    negative_columns = [x for x in columns if re.search('negative', x)]

    for column in negative_columns:
        # Get integer locations of the non-null data in the negative columns
        spot = np.where(df[column].notnull())[0].tolist()
        print(spot)



if __name__ == "__main__":
    main()
