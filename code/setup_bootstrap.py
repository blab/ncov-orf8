'''
Sets up bootstrap values for dN-dS-dStop. Expects data in Usher tree format.
Mandatory inputs:
 --muts, tsv generated from matUtils summary --translate)
 --clades, optional only if you want to include clade info
 --output, path to save muts for bootstrap tsv
'''

import argparse
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

def recode_nodeID(df, clades=False):
    '''
    Recodes node id as an integer.
    '''
    ids = df['node_id'].unique()

    encoder = dict(zip(ids,range(len(ids))))

    df['id'] = df['node_id'].map(encoder)
    if clades == False:
        return df[['id', 'gene','mut_type']]
    elif clades == True:
        return df[['id','gene','mut_type','clade']]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--muts', type=str, required=True, help = 'path to input tsv with translations. Output of matUtils summary --translate')
    parser.add_argument('--clades', type=str, help = 'path to input tsv with clade designation for all nodes')
    parser.add_argument('--output', type = str, required=True, help = 'Path to save muts for bootstrap tsv')
    args = parser.parse_args()

    df = make_mutations_df(args.muts)
    if args.clades:
        with open(args.clades,'r') as f:
            clades = pd.read_csv(f,sep='\t')
        df_nodes = df.merge(clades, left_on='node_id',right_on='node',how='left')
        encoded = recode_nodeID(df_nodes,clades=True)
    else:
        encoded = recode_nodeID(df)

    with open(args.output, 'w') as f:
        encoded.to_csv(f,index=False,sep='\t')
