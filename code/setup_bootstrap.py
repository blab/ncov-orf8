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
from Bio import Phylo

def load_dates(path):
    df = pd.read_csv(path,sep='\t',usecols=['strain','date'],compression='gzip')
    df = df[~df.date.isna()]
    df = df[~df.date.str.contains('?',regex=False)]
    df['date'] = pd.to_datetime(df['date'])
    df.set_index('strain',inplace=True)
    return df

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if not clade.name:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

def get_leaves(node, named):
    leaves = named[node].get_terminals()
    names = [leaf.name for leaf in leaves]
    return names

def get_time(node,named,date_df):
    leaves = get_leaves(node,named)
    dates = date_df.loc[leaves]['date']
    first = min(dates)
    time = max(dates) - first
    return time.days, first

def get_times(df,named,date_df):
    time_vect = np.vectorize(get_time,excluded=[1,2])
    df['days_circulated'],df['date_observed'] = time_vect(df['node_id'],named,date_df)
    return df

def make_mutations_df(path,named,dates):
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

    ## Add time
    time_muts = get_times(all_muts,named,dates)

    ## Add S1
    S1 = time_muts[(time_muts.gene=='S')&(time_muts.residue >= 13)&(time_muts.residue<=685)].reset_index(drop=True)
    S1['gene'] = 'S1'
    final_muts = pd.concat([time_muts,S1])
    final_muts[['gene','mut_type']] = final_muts[['gene','mut_type']].astype("category")

    return final_muts[['node_id','gene','mut_type']]

def recode_nodeID(df):
    '''
    Recodes node id as an integer.
    '''
    ids = df['node_id'].unique()

    encoder = dict(zip(ids,range(len(ids))))

    df['id'] = df['node_id'].map(encoder)

    return df[['id', 'gene','mut_type','date']]

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
    parser.add_argument('--tree', type=str, required=True, help = 'path to input nwk tree')
    parser.add_argument('--dates', type=str, required=True, help = 'path to input tsv.gz containing strain & dates columns for samples in the usher tree')
    parser.add_argument('--output', type = str, required=True, help = 'Path to save muts for bootstrap tsv')
    args = parser.parse_args()

    ## Open tree
    with open(args.tree, 'r') as f:
        tree = Phylo.read(f,'newick')

    ## Load dates
    dates = load_dates(args.dates)

    ## Name tree
    named = tabulate_names(tree)

    df = make_mutations_df(args.muts,named,dates)

    encoded = recode_nodeID(df)



    with open(args.output, 'w') as f:
        encoded.to_csv(f,index=False,sep='\t')
