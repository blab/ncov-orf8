'''
Labels nodes in tree with clade.
Inputs:
    --tree, .nwk tree file
    --clades, .tsv containing clade labels for every sample in tree
Output:
    --output, .tsv containing clade labels for every node in tree
'''
import argparse
from Bio import Phylo
import pandas as pd
import numpy as np

def all_parents(tree):
    parents= {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child.name] = clade
    return parents

def get_mapping(clades):
    return {clade:i for i,clade in enumerate(clades)}

def get_oldest(lineages,mapping):
    lineages = list(lineages)
    ordered = [mapping[lin]for lin in lineages]
    oldest = lineages[np.argmin(ordered)]
    return set([oldest])

def label_up(child,parents,labels,mapped):
    if child.name in parents.keys():
        parent = parents[child.name]
        if parent.name not in labels.keys():
            labels[parent.name] = labels[child.name]
            label_up(parent,parents,labels,mapped)
        elif not labels[child.name].issubset(labels[parent.name]):
            summed = labels[parent.name].union(labels[child.name])
            labels[parent.name] = get_oldest(summed,mapped)
            label_up(parent,parents,labels,mapped)

def tabulate_labels(tree,parents,clades,mapped):
    labels = {}
    leaves = tree.get_terminals()
    for leaf in leaves:
        labels[leaf.name] = set([clades[leaf.name]])
        label_up(leaf,parents,labels,mapped)
    return labels

#def resolve_labels(tree, labels,parents):
#    for clade in tree.find_clades():
#        if len(labels[clade.name])>1:
#            if clade.name in parents.keys():
#                parent = parents[clade.name]
#                labels[clade.name] = labels[parent.name]
#    return labels

def to_df(labelled):
    nodes = []
    labels = []
    for node in labelled.keys():
        if len(labelled[node])==1:
            for e in labelled[node]:
                nodes.append(node)
                labels.append(e)
    new_df = pd.DataFrame({'node':nodes,'clade':labels})
    return new_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--tree', type=str, required=True, help = 'path to .nwk tree')
    parser.add_argument('--clades', type=str, required=True, help = 'path to tsv with clade labels for each sample')
    parser.add_argument('--output', type=str, required=True, help = 'path to output tsv')
    args = parser.parse_args()

    # Nextstrain clades
    nextstrain_clades = [
        '19A',
        '19B',
        '20A',
        '20B',
        '20C',
        '20D',
        '20E (EU1)',
        '20F',
        '20G',
        '20H (Beta,V2)',
        '20I (Alpha,V1)',
        '20J (Gamma,V3)',
        '21A (Delta)',
        '21B (Kappa)',
        '21C (Epsilon)',
        '21D (Eta)',
        '21E (Theta)',
        '21F (Iota)',
        '21G (Lambda)',
        '21H (Mu)',
        '21I (Delta)',
        '21J (Delta)',
        '21K (Omicron)',
        '21L (Omicron)',
        '21M (Omicron)',
        '22A (Omicron)',
        '22B (Omicron)',
        '22C (Omicron)',
        '22D (Omicron)',
        '22E (Omicron)',
        '22F (Omicron)'
    ]

    # Load files
    with open(args.tree, 'r') as f:
        tree = Phylo.read(f,'newick')
    with open(args.clades,'r') as f:
        df = pd.read_csv(f, sep='\t')

    # Get parents & clades
    parents = all_parents(tree)
    clades = {df.at[row,'sample']:df.at[row,'annotation_1'] for row in df.index}
    mapped = get_mapping(nextstrain_clades)

    # Backwards traversal
    labelled = tabulate_labels(tree,parents,clades,mapped)

    # Forwards traversal
    #labelled = resolve_labels(tree,labelled_up,parents)

    # Write to df
    label_df = to_df(labelled)
    with open(args.output,'w') as f:
        label_df.to_csv(f,sep='\t',index=False)
