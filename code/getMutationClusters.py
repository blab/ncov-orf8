'''
Takes usher .nwk, usher nodes .tsv (generated from matUtils summary --node-stats)
& usher translate .tsv (generated from matUtils summary --translate)
(Above options generated from matUtils v0.6.2)

If --nested, will call all mutation clusters for that gene. Otherwise, finds
non-nested mutation clusters for each gene, i.e. clusters with a mutation in
that gene & only one mutation.

Returns TSV that at min contains starting node of that cluster, number of
samples in cluster, days circulated, and the parent.
'''

import argparse
from Bio import Phylo
import pandas as pd
import numpy as np
import pathlib

def load_dates(path):
    df = pd.read_csv(path,sep='\t',usecols=['strain','date'],compression='gzip')
    df = df[~df.date.isna()]
    df = df[~df.date.str.contains('?',regex=False)]
    df['date'] = pd.to_datetime(df['date'])
    df.set_index('strain',inplace=True)
    return df

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child.name] = clade.name
    return parents

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if not clade.name:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

def find_nodes(muts, gene):
    '''
    Finds nodes containing mutation in gene
    '''
    aa = muts['aa_mutations']
    nt = muts['nt_mutations']
    codons = muts['codon_changes']
    idxs = [idx for idx in aa.index if gene+':' in aa[idx]]
    counts = [aa[idx].count(gene) for idx in idxs]
    gene_idxs = [idx for idx, c in zip(idxs,counts) if c == 1]
    split_aa = [aa[idx].split(';') if ';' in aa[idx] else [aa[idx]]for idx in gene_idxs]
    locations = [i for lst in split_aa for i,j in enumerate(lst) if gene in j]
    AA = [v[i] for i,v in zip(locations,split_aa)]
    gene_aa = [v.split(':')[1] for v in AA]
    ogs = [g_aa[0] for g_aa in gene_aa]
    news = [g_aa[-1] for g_aa in gene_aa]
    mut_type = []
    for og, new in zip(ogs, news):
        if og == new:
            mut_type.append('synonymous')
        elif new == '*':
            mut_type.append('nonsense')
        elif og != '*':
            mut_type.append('missense')
        else:
            mut_type.append('undoStop')
    branch_muts = nt[gene_idxs].str.count(';') + 1
    nodes = pd.DataFrame({'node_id': muts['node_id'][gene_idxs], 'aa_mutations':aa[gene_idxs], 'nt_mutations':nt[gene_idxs], 'codon_change':codons[gene_idxs],'mut_type':mut_type,'branch_muts_counts': branch_muts, 'gene':gene}).reset_index(drop=True)
    return nodes

def get_ancestry(node, parentDict):
    parentals = set()
    parent = node
    while parent in parentDict.keys():
        parent = parentDict[parent]
        parentals.add(parent)
    return parentals

def find_nested(df, parentDict):
    nodes = set(df['node_id'])
    remove = set()

    indexed = pd.Series(df['node_id'].index, index=df['node_id'])

    for idx in df.index:
        anc = get_ancestry(df['node_id'][idx],parentDict)
        nested = anc.intersection(nodes)
        if len(nested):
            nested_idx = indexed[list(nested)]
            remove.update(nested_idx)
    return remove

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

def add_parents(parents, df):
    the_parents = [parents[node] for node in df['node_id']]
    df['parent'] = the_parents
    return df

def add_cluster(df):
    length = len(df)
    df['cluster'] = np.arange(1,length+1)
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--tree', type=str, required=True, help = 'path to input nwk tree')
    parser.add_argument('--stats', type=str, required=True, help = 'path to input tsv with node stats. Output of matUtils summary --node-stats')
    parser.add_argument('--muts', type=str, required=True, help = 'path to input tsv with translations. Output of matUtils summary --translate')
    parser.add_argument('--dates', type=str, required=True, help = 'path to input tsv containing strain & dates columns for samples in the usher tree')
    parser.add_argument('--nested', action="store_true", help = 'Identifies nested clades')
    parser.add_argument('--genes', nargs='+', type=str, default='ORF1a ORF1b S ORF3a E M ORF6 ORF7a ORF7b ORF8 N ORF9b', help = 'Genes for which to generate mutation clusters')
    parser.add_argument('--outputdir', type=str, required=True, help = 'path to output directory')
    args = parser.parse_args()

    ## Load inputs
    with open(args.tree, 'r') as f:
        tree = Phylo.read(f,'newick')

    with open(args.muts, 'r') as f:
        muts = pd.read_csv(f,sep='\t')

    with open(args.stats, 'r') as f:
        stats = pd.read_csv(f,sep='\t')

    dates = load_dates(args.dates)

    # Map node ids to tree object
    named = tabulate_names(tree)

    # Create dictionary of parent for each node
    parents = all_parents(tree)


    path = pathlib.Path(args.outputdir)
    path.mkdir(parents=True, exist_ok=True)

    for gene in args.genes: ## ADD code for S1
        print('Starting ' + gene)

        # Find nodes with a mutation in that gene & classify mutation type
        nodes = find_nodes(muts,gene)

        # Add in size of that node + other stats
        merged = nodes.merge(stats,how='left',left_on='node_id',right_on='node').drop(columns=['node'])

        # Identify nested nodes & remove them
        if args.nested:
            df = merged
        else:
            to_remove = find_nested(merged,parents)
            no_nested = merged.drop(index=list(to_remove))
            df = no_nested

        # Add in circulation times
        plus_times = get_times(df,named,dates)

        # Add parents
        plus_parents = add_parents(parents,plus_times)

        # Add cluster
        clustered = add_cluster(plus_parents)

        # Save output
        with open(args.outputdir + '/' + gene + '_clades.tsv', 'w') as f:
            clustered.to_csv(f,sep='\t',index=False)
