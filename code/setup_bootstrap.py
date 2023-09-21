'''
Sets up bootstrap values for dN-dS-dStop. Expects data in Usher tree format.
Mandatory inputs:
 --muts, tsv generated from matUtils summary --translate)
 --output, path to save muts for bootstrap tsv
 --chunkSizes, folder to save .txt with chunk size for each chunk
'''

import argparse
#import os
import pandas as pd
import numpy as np

def make_mutations_df(path):
    '''
    Constructs dataframe with all mutations in Usher tree for each gene.
    '''
    with open(path,'r') as f:
        muts = pd.read_csv(f,sep='\t',usecols=['node_id','aa_mutations'])

    ## Split up nodes with multiple mutations into single mutations
    muts['aa_mutations'] = muts.aa_mutations.str.split(';')
    all_muts = muts.explode(['aa_mutations'],ignore_index=True)

    ## Label with gene
    all_muts[['gene','aa_mutation']] = all_muts['aa_mutations'].str.split(':',expand=True)
    all_muts.drop(columns=['aa_mutations'],inplace=True)

    ## Get residue
    all_muts['residue'] = pd.to_numeric(all_muts.aa_mutation.str[1:-1])

    ## Classify mutation type
    all_muts['mut_type'] = np.where(all_muts.aa_mutation.str[-1]==all_muts.aa_mutation.str[0],'synonymous','missense')
    all_muts['mut_type'] = np.where(all_muts.aa_mutation.str[-1]=='*','nonsense',all_muts['mut_type'])
    all_muts['mut_type'] = np.where(all_muts.aa_mutation.str[0]=='*','undoStop',all_muts['mut_type'])

    ## Add S1
    S1 = all_muts[(all_muts.gene=='S')&(all_muts.residue >= 13)&(all_muts.residue<=685)].reset_index(drop=True)
    S1['gene'] = 'S1'
    final_muts = pd.concat([all_muts,S1])
    final_muts[['gene','mut_type']] = final_muts[['gene','mut_type']].astype("category")

    return final_muts[['node_id','gene','mut_type']]

def recode_nodeID(df):
    '''
    Recodes node id as an integer.
    '''
    ids = df['node_id'].unique()

    encoder = dict(zip(ids,range(len(ids))))

    df['id'] = df['node_id'].map(encoder)

    return df[['id', 'gene','mut_type']]

#def draw_samples(ids,n,chunksize):
#    '''
#    Draws random samples with replacement of node ids for bootstrap.
#    '''
#    sample_ids = {iterat:choices(ids,k=chunksize) for iterat in range(n)}
#    return sample_ids#

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--muts', type=str, required=True, help = 'path to input tsv with translations. Output of matUtils summary --translate')
    parser.add_argument('--output', type = str, required=True, help = 'Path to save muts for bootstrap tsv')
    #parser.add_argument('--chunkSizes', type = str, required=True, help = 'folder to save chunk sizes in')
    args = parser.parse_args()

    #size = 10000

    #CHECK_FOLDER = os.path.isdir(args.chunkSizes))

    # If folder doesn't exist, then create it.
    #if not CHECK_FOLDER:
    #    os.makedirs(args.chunkSizes)

    df = make_mutations_df(args.muts)

    encoded = recode_nodeID(df)

    #del df

    with open(args.output, 'w') as f:
        encoded.to_csv(f,index=False,sep='\t')

    #ids = encoded['id'].unique()
#    length = len(ids)#

#    del encoded
#    del ids#

#    n_chunks = int(np.floor(length/size))
#    remain = length % size#

#    chunkSizes = [size] * n_chunks + [remain]
#    for i, chunk in enumerate(chunkSizes):
#        with open(args.chunkSizes + 'chunk' + str(i)+'.txt', 'w') as fp:
#                # write each item on a new line
#                fp.write("%s\n" % chunk)
#    print('Done writing chunkSizes')
