'''
Sets up bootstrap values for dN-dS-dStop. Expects data in Usher tree format.
Mandatory inputs:
 --muts, tsv generated from matUtils summary --translate)
 --clade, optional only if you want to include clade info
 --output, path to save muts for bootstrap tsv
'''

import argparse
import pandas as pd
import numpy as np
import os

def recode_nodeID(df):
    '''
    Recodes node id as an integer.
    '''
    ids = df['id'].unique()

    encoder = dict(zip(ids,range(len(ids))))

    df['id'] = df['id'].map(encoder)
    return df[['id','gene','mut_type','clade']]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--muts', type=str, required=True, help = 'path to input muts tsv')
    parser.add_argument('--clade', type=str, required=True, help = 'clade')
    parser.add_argument('--output', type = str, required=True, help = 'Path to save muts for clade bootstrap')
    args = parser.parse_args()

    with open(args.muts, 'r') as f:
        df = pd.read_csv(f, sep='\t')
    filt = df[df.clade==args.clade].copy()
    encoded = recode_nodeID(filt)

    with open(args.output, 'w') as f:
        encoded.to_csv(f,index=False,sep='\t')
