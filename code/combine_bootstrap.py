'''
Runs bootstrap of dN-dS-dStop. Expects data in Usher tree format.
Mandatory inputs:
 --ref, reference Genbank file of gene locations
 --muts, tsv generated from matUtils summary --translate)
 --counts, path to save bootstrap counts tsv
 --diffs, path to save dN-dS-dStop tsv
 Optional inputs:
  --number,  number of bootstrap iterations
  --ci, confidence interval to use
'''

import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
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

def get_counts(df, genes):
    '''
    Get observed counts
    '''
    geneList = []
    mutTypes = []
    counts = []
    for gene in genes:
        for mut_type in df.mut_type.unique():
            mutTypes.append(mut_type)
            count = len(df[(df.gene==gene)&(df.mut_type==mut_type)])
            counts.append(count)
            geneList.append(gene)
    result = pd.DataFrame({'gene':geneList,'mut_type':mutTypes,'count':counts})
    return result

def get_differences(counts,denominators, genes):
    '''
    Given counts & denominators, get dS, dN, dStop & dN_dS and dStop_dS values.
    '''
    df = counts.copy()
    cols = list(df.columns)
    cols.remove('mut_type')
    cols.remove('count')
    df = df.pivot_table(index=cols, columns='mut_type',values='count').reset_index()
    for gene in genes:
        df.loc[(df.gene==gene),'dS'] = df.loc[(df.gene==gene),'synonymous']/denominators['synonymous'][gene]
        df.loc[(df.gene==gene),'dN'] = df.loc[(df.gene==gene),'missense']/denominators['missense'][gene]
        if 'nonsense' in df.columns:
            df.loc[(df.gene==gene),'dStop'] = df.loc[(df.gene==gene),'nonsense']/denominators['nonsense'][gene]
        else:
            df.loc[(df.gene==gene),'dStop'] = 0/denominators['nonsense'][gene]
    df['dN_dS'] = df['dN']/df['dS']
    df['dStop_dS'] = df['dStop']/df['dS']
    return df

def get_percentile(df,p,suffix):
    '''
    Returns percentile with suffix for all genes.
    '''
    transformed = df.groupby('gene').quantile(p).add_suffix(suffix).reset_index()
    return transformed

def get_ci(df,ci):
    '''
    Returns confidence intervals from bootstrap diffs.
    '''
    if 'undoStop' in df.columns:
        if 'nonsense' in df.columns:
            new_df = df.drop(columns=['iteration','missense','synonymous','nonsense','undoStop'])
    elif 'nonsense' in df.columns:
        new_df = df.drop(columns=['iteration','missense','synonymous','nonsense'])
    else:
        new_df = df.drop(columns=['iteration','missense','synonymous'])
    minimum = (100-ci)/200
    maximum = 1-minimum
    minimums = get_percentile(new_df,minimum,'_min')
    maximums = get_percentile(new_df,maximum,'_max')
    intervals = minimums.merge(maximums, on='gene')
    return intervals


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--ref', type=str, required=True, help = 'path to reference genbank file')
    parser.add_argument('--muts', type=str, required=True, help = 'path to input muts by mutType')
    parser.add_argument('--ci', type=int, default=95, help = 'confidence interval for bootstrap')
    parser.add_argument('--counts', nargs='+', type = str, required=True, help = 'Path to save counts tsv')
    parser.add_argument('--degenerate', default=False, action='store_true', help = 'flag to use substitution matrix from 4-fold degenerate sites')
    parser.add_argument('--diffs', type = str, required=True, help = 'Path to save dS-dN-dStop tsv')
    args = parser.parse_args()

    if args.degenerate == True:
        # substitution matrix from 4-fold degenerate expected counts from Bloom & Neher 2023:
        # https://www.biorxiv.org/content/10.1101/2023.01.30.526314v2.full
        # substitution matrix just normalized from expected counts
        sub_matrix = {'A': {'T': 0.01427,'C': 0.00628,'G': 0.05130},
            'T': {'A': 0.00978,'C': 0.04658,'G': 0.00609},
            'C': {'A': 0.02716,'T': 0.39887,'G': 0.00839},
            'G': {'A': 0.11016,'T': 0.30833,'C': 0.01279}}
    else:
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

    # Aggregate bootstrap counts
    for i, count in enumerate(args.counts):
        with open(count, 'r') as f:
            countedF = pd.read_csv(f,sep='\t')
        if i == 0:
            counted = countedF
        else:
            lastIT = max(counted['iteration'])
            countedF['iteration'] = countedF['iteration'] + lastIT + 1
            counted = pd.concat([counted,countedF])

    bootstrap_diffs = get_differences(counted,denominators,all_genes)

    bootstrap_cis = get_ci(bootstrap_diffs, args.ci)

    ## Load muts
    with open(args.muts, 'r') as f:
        the_muts = pd.read_csv(f, sep='\t')

    # Get counts actually observed
    observed_counts = get_counts(the_muts,all_genes)

    # Get the observed dN-dS-dStop values
    observed_diffs = get_differences(observed_counts,denominators,all_genes)

    # Combine observed with CIs
    to_save = observed_diffs.merge(bootstrap_cis,on='gene')
    to_save['CI'] = args.ci

    to_save.to_csv(args.diffs,sep='\t',index=False)
