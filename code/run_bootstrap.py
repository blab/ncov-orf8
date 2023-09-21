'''
Sets up bootstrap values for dN-dS-dStop. Expects data in Usher tree format.
Mandatory inputs:
 --muts, tsv with node_id, gene, and mut_type
 --samples, json with samples for each iteration
 --counts, path to save counts tsv
'''
import argparse
import os
import pandas as pd
from random import choices

def draw_samples(length,n):
    '''
    Draws random samples of ids with replacement.
    '''
    sample_ids = {iterat:choices(range(length),k=length) for iterat in range(n)}
    return sample_ids

def draw_counts(samples,df):
    '''
    Returns df of counts per gene by mutation type for each iteration.
    '''
    for i, k in enumerate(samples.keys()):
        filt = df.loc[samples[k],:]
        count = filt.groupby(by=['gene','mut_type']).size().reset_index()
        count.rename(columns={0:'count'},inplace=True)
        count['iteration'] = k
        if i == 0:
            bootstrap = count
        else:
            bootstrap = pd.concat([bootstrap,count])
    return bootstrap

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--muts', type=str, required=True, help = 'path to tsv with node_id, gene, and mut_type')
    parser.add_argument('--size',type=int,default=20, help = 'number of bootstrap iterations to run per chunk')
    parser.add_argument('--counts', type = str, required=True, help = 'path to save counts tsv')
    args = parser.parse_args()

    # Check if directory exists
    pathList = args.counts.split('/')
    directory = '/'.join(pathList[0:2])

    check_directory = os.path.isdir(directory)

    # If folder doesn't exist, then create it.
    if not check_directory:
        os.makedirs(directory)

    ## Load df
    with open(args.muts, 'r') as f:
        df = pd.read_csv(f, sep='\t').set_index('id')
    df[['gene','mut_type']] = df[['gene','mut_type']].astype("category")

    #with open(args.size, 'r') as f:
    #    size = pd.read_table(f, header=None,names='chunkSize')

    length = len(df.index.unique())

    chunksize = args.size

    sampled = draw_samples(length,chunksize)

    ## Get counts
    counted = draw_counts(sampled, df)

    ## Save counts
    with open(args.counts, 'w') as f:
        counted.to_csv(f,sep='\t',index=False)
