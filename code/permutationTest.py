'''
Runs permutation test on tips in an Usher tree for cluster circulation time & size.
Input a tsv file with the following columns:
 - leaf_count, number of tips in clade
 - days_circulated, days in between the first sample from that clade & the last sample.
 - cladeType, nomut or mut. Mut is what was observed, nomut is sibling clades.
 - mutType, stops, nonsynonymous, nonsynonymous
 - cluster
 - gene

 Outputs a dictionary with test results.
'''

import argparse
import pandas as pd
import numpy as np
from scipy import stats as st
import json

def load_df(path):
    with open(path, 'r') as f:
        df = pd.read_csv(f, sep='\t',usecols=['cluster','cladeType','leaf_count','days_circulated','mutType','node_id','mut_count'])
    df['mutType'] = df['mutType'].astype("category")
    df['cladeType'] = df['cladeType'].astype("category")
    df[["cluster", "leaf_count","mut_count"]] = df[["cluster", "leaf_count","mut_count"]].apply(pd.to_numeric,downcast="unsigned")
    df["days_circulated"] = pd.to_numeric(df["days_circulated"],downcast="float")
    return df

def get_permutations(df, size):
    indices = [df[df.cluster==cluster].index for cluster in df.cluster.unique()]
    choices = np.array([np.random.default_rng().choice(idx,size=size,replace=True) for idx in indices],dtype=int)
    choicesList = [tuple(choices[:,col]) for col in range(size)]
    deduped = list(set(choicesList))
    arr = np.array([list(values) for values in deduped]).astype(int)
    return arr

def permutation_test_absolute(df,size,var,direct='greater'):
    trueRows = df.loc[df.cladeType=='mut'].index
    trueVar = df.loc[trueRows,var]
    options = get_permutations(df,int(size*1.5))
    values = np.vectorize(df[var].get)(options[:size,:])
    if var == 'leaf_count':
        trueMean = st.gmean(trueVar)
        means = st.gmean(values,axis=1) ### FIGURE OUT WAYS TO SAVE ON size
    else:
        trueMean = np.mean(trueVar)
        means = np.mean(values, axis=1) ## FIGURE OUT WAYS TO SAVE
    if direct == 'greater':
        total = len(means[means > trueMean])
    else:
        total = len(means[means<trueMean])
    pvalue = total/size
    return pvalue, means, trueMean #,prop

def permut_gene(df, size):
    d = {}
    new_df = df[df.mut_count>0]
    for mut in ['nonsense','missense','synonymous']:
        filt = new_df[new_df.mutType==mut]
        d[mut] = {}
        for var in ['leaf_count','days_circulated']:
            d[mut][var] = {}
            p,means,trueMean = permutation_test_absolute(filt,size,var)
            d[mut][var]['pvalue'] = float(p)
            d[mut][var]['means'] = means.tolist()
            d[mut][var]['mean'] = float(trueMean)
    return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', type=str, required=True, help = 'path to input tsv')
    parser.add_argument('--permuts', type=int, default=1000, help= 'Number of permutations to do.')
    parser.add_argument('--output', type=str, required=True, help = 'path to output jsson')
    args = parser.parse_args()

    ## load TSV
    df = load_df(args.input)

    result = permut_gene(df,args.permuts)

    with open(args.output, 'w') as f:
        json.dump(result,f)
