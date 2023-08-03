'''
Creates tsv(s) for permutation test which tests if clades with each type of
mutation (synonymous, missense, nonsense) are the same size as their sibling
clades on the same genetic backbone.

Inputs:
 --cladesDir, directory containing clade files generated from
 getMutationClusters.py, which has mutation clusters
 --stats,  usher nodes .tsv (generated from matUtils summary --node-stats)
 --tree, usher nwk tree
 --dates, df containing all the dates
 --outputdir, AKA where this stuff going to go at the end
'''

import argparse
import os
import pandas as pd
from Bio import Phylo
import numpy as np
import getMutationClusters as gmc

def load_clades(path):
    with open(path, 'r') as f:
        df = pd.read_csv(f, sep='\t',usecols=['node_id','cluster','parent','mut_type'])
    df['mut_type'] = df['mut_type'].astype("category")
    df["cluster"] = df["cluster"].apply(pd.to_numeric,downcast="unsigned")
    return df

def get_children(df,named):
    ids = []
    mutType = []
    clusters = []
    cladeTypes = []
    for mut in ['synonymous','missense','nonsense','undoStop']:
        filt = df[df.mut_type==mut]
        observed = filt['node_id']
        parented = filt['parent']
        clustered = filt['cluster']
        for node,parent,cluster in zip(observed,parented,clustered):
            parent_tree = named[parent]
            kids = [child.name for child in parent_tree.clades]
            loc = kids.index(node)
            n = len(kids)
            cladeType = ['nomut']*n
            cladeType[loc] = 'mut'
            ids.extend(kids)
            mutType.extend([mut]*n)
            clusters.extend([cluster]*n)
            cladeTypes.extend(cladeType)
    children = pd.DataFrame({'node_id':ids,'mutType':mutType,'cluster':clusters,'cladeType':cladeTypes})
    return children

def remove_noSiblings(df):
    remove = set()
    for mut in df.mutType.unique():
        filt = df[df.mutType==mut]
        for cluster in filt.cluster.unique():
            clustered = filt[filt.cluster==cluster]
            if len(clustered) == 1:
                remove.update(clustered.index)
    dropped = df.drop(index=list(remove))
    return dropped

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--cladesDir', type=str, required=True, help = 'path to directory containing output of getMutationClusters.py')
    parser.add_argument('--tree', type=str, required=True, help = 'path to input nwk tree')
    parser.add_argument('--stats', type=str, required=True, help = 'path to input tsv with node stats. Output of matUtils summary --node-stats')
    parser.add_argument('--muts', type=str, required=True, help = 'path to input tsv with translations. Output of matUtils summary --translate')
    parser.add_argument('--dates', type=str, required=True, help = 'path to input tsv containing strain & dates columns for samples in the usher tree')
    parser.add_argument('--outputdir', type=str, required=True, help = 'path to output directory')
    args = parser.parse_args()

    ## Load inputs
    with open(args.tree, 'r') as f:
        tree = Phylo.read(f,'newick')

    with open(args.stats, 'r') as f:
        stats = pd.read_csv(f,sep='\t')

    with open(args.muts, 'r') as f:
        muts = pd.read_csv(f,sep='\t')

    dates = gmc.load_dates(args.dates)

    # Map node ids to tree object
    named = gmc.tabulate_names(tree)

    #for gene in ['ORF1a','ORF1b','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF9b']: ## ADD code for S1
    for gene in ['ORF8']:
        print('Starting ' + gene)

        # Load clade df
        clade = load_clades(args.cladesDir + gene + '_clades.tsv')

        # Get children
        all_children = get_children(clade, named)

        # Add on node stats
        children_size = all_children.merge(stats,how='left', left_on = 'node_id', right_on = 'node')

        # remove branches of length == 0
        children_noPolytomy = children_size[children_size.mut_count>0]

        # Remove nodes without any sibs
        children_shrunk = remove_noSiblings(children_noPolytomy)

        # Add times
        children_complete = gmc.get_times(children_shrunk,named,dates)

        # Add mutation info
        children_muts = children_complete.merge(muts,how='inner',on='node_id')
        children_muts['branch_muts_count'] = children_muts['nt_mutations'].str.count(';') + 1

        # Save
        children_muts.to_csv(args.outputdir + gene + '_permutations.tsv', sep='\t',index=False)
