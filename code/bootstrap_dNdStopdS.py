'''
Runs bootstrap of dN-dS-dStop. Expects data in Usher tree format.
Mandatory inputs:
 --ref, reference Genbank file of gene locations
 --muts, tsv generated from matUtils summary --translate)
 --counts, path to save bootstrap counts tsv
 --diffs, path to save dN-dS-dStop tsv
 Optional inputs:
  --number,  number of bootstrap iterations
  --ci, confidene interval to use
'''

import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
from random import choices
from itertools import compress

def get_locations(ref):
    '''
    Returns dictionary of all gene locations.
    '''
    location_by_gene = {}
    with open(ref,'r') as reference_handle:
        for record in SeqIO.parse(reference_handle, "genbank"):
            wuhan_seq = record.seq
            for feature in record.features:
                if feature.type == 'CDS':
                    gene = feature.qualifiers['gene'][0]
                    location = feature.location
                    location_by_gene[gene] = location
    return location_by_gene,wuhan_seq

def get_codons(ldict,sequence):
    '''
    Returns codons by gene.
    Extracts sequence of each gene & finds constituent codon sequences
    '''
    codons_by_gene = {}
    for g, l in ldict.items():
        gene_sequence = l.extract(sequence)
        #split into codons
        #indexing is 0-based, but codons count from 1.
        # So use codons_by_gene[gene][codon-1]
        gene_codons = [gene_sequence[i:i+3] for i in range(0, len(gene_sequence), 3)]
        codons_by_gene[g] = gene_codons
    return codons_by_gene

def calculate_expected_sites(cdict,submatrix,nts):
    '''
    Walks through codons in each gene and returns every possible mutation
    * probability of mutatation.
    '''
    syn_denominators ={}
    nonsyn_denominators ={}
    stop_denominators = {}

    for gene, codons in cdict.items():
        nonsyn_sites_per_gene = 0
        syn_sites_per_gene = 0
        stop_sites_per_gene = 0

        #walk through codons
        for codon in codons:
            original_aa = codon.translate()

            #for each position in the codon
            for i in range(len(codon)):
                #wuhan nt
                nt = codon[i]
                #introduce each other mutation
                for mut in [x for x in nts if x!=nt]:
                    mut_codon = codon[:i]+mut+codon[i+1:]
                    mut_aa = mut_codon.translate()
                    #find whether nonsynonymous
                    if mut_aa != original_aa:
                        #find whether this is a stop codon
                        if mut_aa =='*':
                            #get the probability of this mutation
                            #add to the total number of stop sites
                            stop_sites_per_gene+=submatrix[nt][mut]
                        else:
                            #or the total number of nonsyn sites
                            nonsyn_sites_per_gene+=submatrix[nt][mut]
                    elif mut_aa == original_aa:
                        syn_sites_per_gene+=submatrix[nt][mut]

        syn_denominators[gene] = syn_sites_per_gene
        nonsyn_denominators[gene] = nonsyn_sites_per_gene
        stop_denominators[gene] = stop_sites_per_gene
        denominators = {'synonymous':syn_denominators,'missense':nonsyn_denominators,'nonsense':stop_denominators}
    return denominators

def make_mutations_df(path):
    '''
    Constructs dataframe with all mutations in Usher tree for each gene.
    '''
    with open(path,'r') as f:
        muts = pd.read_csv(f,sep='\t')

    ## Split up nodes with multiple mutations into single mutations
    muts['aa_mutations'] = muts.aa_mutations.str.split(';')
    muts['nt_mutations'] = muts.nt_mutations.str.split(';')
    muts['codon_changes'] = muts.codon_changes.str.split(';')
    all_muts = muts.explode(['aa_mutations','nt_mutations','codon_changes'],ignore_index=True)

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

    return final_muts

def draw_samples(df,n):
    '''
    Draws random samples with replacement of node ids for bootstrap.
    Returns indexes of random samples.
    '''
    ids = df['node_id'].unique()
    idxs = {ID:np.flatnonzero(df['node_id']==ID) for ID in ids}
    length  = len(ids)
    sample_ids = {iterat:choices(ids,k=length) for iterat in range(n)}
    sample_idxs = {iterat:np.concatenate([idxs[ID] for ID in sample_ids[iterat]]) for iterat in range(n)}
    return sample_idxs

def draw_genes(samples,genes,df):
    '''
    Returns indexes of random samples for each gene.
    '''
    bootstrap_genes = {}
    for k in samples.keys():
        bootstrap_genes[k] = {}
        for gene in genes:
            filt = df.iloc[samples[k],4].values == gene
            subset = compress(samples[k],filt)
            bootstrap_genes[k][gene] = list(subset)
    return bootstrap_genes

def draw_counts(gene_samples,genes,df):
    '''
    Returns df of counts per gene by mutation type for each iteration.
    '''
    iterations = []
    the_genes = []
    mutTypes = []
    counts = []
    for k in gene_samples.keys():
        iterations.extend([k]*len(genes)*3)
        for gene in genes:
            the_genes.extend([gene]*3)
            for mut_type in ['synonymous','missense','nonsense']:
                filt = df.iloc[gene_samples[k][gene],7].values == mut_type
                subset = compress(gene_samples[k][gene],filt)
                count = sum(1 for x in subset)
                mutTypes.append(mut_type)
                counts.append(count)
    bootstrap = pd.DataFrame({'bootstrap':iterations,'gene':the_genes,'mutType':mutTypes,'count':counts})
    return bootstrap

def get_counts(df, genes):
    '''
    Get observed counts
    '''
    all_syn = {gene:len(df[(df.gene==gene)&(df.mut_type=='synonymous')]) for gene in genes}
    all_missense = {gene:len(df[(df.gene==gene)&(df.mut_type=='missense')]) for gene in genes}
    all_nonsense = {gene:len(df[(df.gene==gene)&(df.mut_type=='nonsense')]) for gene in genes}
    counts = {'synonymous':all_syn,'missense':all_missense,'nonsense':all_nonsense}
    return counts

def get_differences(counts,denominators,genes):
    '''
    Given counts & denominators, get dS, dN, dStop values.
    '''
    diffs = {}
    for mutType in counts.keys():
        diffs[mutType] = {gene:counts[mutType][gene]/denominators[mutType][gene] for gene in genes}
    return diffs

def get_ci_counts(df,genes,ci):
    '''
    Given df of bootstrap counts, returns dictionary of min & max ci counts.
    '''
    min_counts = {}
    max_counts = {}
    for mutType in df['mutType'].unique():
        min_counts[mutType] = {}
        max_counts[mutType] = {}
        for gene in genes:
            values = df[(df.gene==gene) & (df.mutType==mutType)]['count']
            minimum = 100-ci/2
            maximum = 100-minimum
            min_counts[mutType][gene] = np.percentile(values,minimum)
            max_counts[mutType][gene] = np.percentile(values,maximum)
    return min_counts, max_counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--ref', type=str, required=True, help = 'path to reference genbank file')
    parser.add_argument('--muts', type=str, required=True, help = 'path to input tsv with translations. Output of matUtils summary --translate')
    parser.add_argument('--number',type=int, default=1000, help = 'Number of bootstrap iterations to perform')
    parser.add_argument('--ci', type=int, default=95, help = 'confidence interval for bootstrap')
    parser.add_argument('--counts', type = str, required=True, help = 'Path to save counts tsv')
    parser.add_argument('--diffs', type = str, required=True, help = 'Path to save dS-dN-dStop tsv')
    args = parser.parse_args()

    #substitution matrix from https://www.sciencedirect.com/science/article/pii/S1567134821004287
    #format is {from_A:{to_T:x, to_C:x, to_G:x}, from_T:{to_A:x, ...}, ...}

    sub_matrix = {'A': {'T': 0.0383, 'C': 0.0219, 'G': 0.0747},
          'T': {'A': 0.0356, 'C': 0.2085, 'G': 0.0234},
          'C': {'A': 0.0356, 'T': 0.3648, 'G': 0.0234},
          'G': {'A': 0.1138, 'T': 0.0383, 'C': 0.0219}}

    all_nts = ['A', 'T', 'C', 'G']

    # Load genbank file and get gene locations
    locs,refseq = get_locations(args.ref)

    all_genes = list(locs.keys())

    # Create dictionary with codons for each gene
    codons = get_codons(locs,refseq)

    # Calculate expected number of sites
    denominators = calculate_expected_sites(codons,sub_matrix, all_nts)

    # Load and clean mutations df
    the_muts = make_mutations_df(args.muts)

    # Draw sample iterations for bootstrap
    sampled = draw_samples(the_muts,args.number)

    # Split out sample idxs by gene
    sampled_by_gene = draw_genes(sampled,all_genes,the_muts)

    # To save memory, remove no longer needed dictionary
    del sampled

    # Get counts for each bootstrap iteration
    counted = draw_counts(sampled_by_gene,all_genes,the_muts)

    # To save momery, remove no longer needed dictionary
    del sampled_by_gene

    # Save bootstrap to TSV
    with open(args.counts,'w') as f:
        counted.to_csv(f,sep='\t',index=False)

    # Get ci counts from bootstrap
    minimum_counts, maximum_counts = get_ci_counts(counted,all_genes,args.ci)

    # To save memory, remove no longer needed df
    del counted

    # Get counts actually observed
    observed_counts = get_counts(the_muts,all_genes)

    # Get the observed dN-dS-dStop values
    observed_diffs = get_differences(observed_counts,denominators,all_genes)

    # Get CI dN-dS-dStop values
    minimum_diffs = get_differences(minimum_counts,denominators,all_genes)
    maximum_diffs = get_differences(maximum_counts,denominators,all_genes)

    to_save = []
    for gene in all_genes:
        to_save.append({
        'gene':gene,
        'dS': observed_diffs['missense'][gene],
        'dN': observed_diffs['missense'][gene],
        'dStop': observed_diffs['nonsense'][gene],
        'dN_dS':observed_diffs['missense'][gene]/observed_diffs['synonymous'][gene],
        'dStop_dS':observed_diffs['nonsense'][gene]/observed_diffs['missense'][gene],
        'dS_min': minimum_diffs['missense'][gene],
        'dN_min': minimum_diffs['missense'][gene],
        'dStop_min': minimum_diffs['nonsense'][gene],
        'dN_dS_min':minimum_diffs['missense'][gene]/minimum_diffs['synonymous'][gene],
        'dStop_dS_min':minimum_diffs['nonsense'][gene]/minimum_diffs['missense'][gene],
        'dS_max': maximum_diffs['missense'][gene],
        'dN_max': maximum_diffs['missense'][gene],
        'dStop_max': maximum_diffs['nonsense'][gene],
        'dN_dS_max':maximum_diffs['missense'][gene]/maximum_diffs['synonymous'][gene],
        'dStop_dS_max':maximum_diffs['nonsense'][gene]/maximum_diffs['missense'][gene],
        'CI':args.ci
        })

    df_to_save = pd.DataFrame(to_save)

    df_to_save.to_csv(args.diffs,sep='\t',index=False)
