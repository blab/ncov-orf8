'''
This script uses a .nwk tree & tsv containing gene KO's called by find_ko.py to
call conservative KO clusters from trees.
'''

import argparse
import pandas as pd
import numpy as np
from Bio import Phylo

def replaceNaNs(df,genes):
    for gene in genes:
        df[gene+'_koType'] = np.where(df[gene+'_koType'].isnull(),'None',df[gene+'_koType'])
        df[gene+'_misStops'] = np.where(df[gene+'_misStops'].isnull(),'None',df[gene+'_misStops'])
        df[gene+'_deletions'] = df[gene+'_deletions'].str.replace(",\s+(?=[^\(]*\))",':', regex=True)
        df[gene+'_deletions'] = df[gene+'_deletions'].str.replace('[\[\]()]', '', regex=True)
    return df

def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if not clade.name:
            clade.name = str(idx)
        names[clade.name] = clade
    return names

def compareDels(deleted):
    linear = sorted(list(deleted))
    rolled = np.roll(linear,1)
    diffed = np.where((linear-rolled)>1)[0]
    if len(diffed):
        starts = []
        ends = []
        sizes = []
        nstart = linear[0]
        for value in diffed:
            nend = linear[value-1]
            sizes.append(nend-nstart+1)
            starts.append(nstart)
            ends.append(nend)
            nstart = linear[value]
        j = np.argmax(sizes)
        start = starts[j]
        end = ends[j]
    elif len(deleted):
        start = linear[0]
        end = linear[-1]
    else:
        start = 0
        end = 0
    length = int(end-start + 1)
    return length,start,end

def extractDels(meta,clade,gene):
    locs = meta.loc[meta.strain==clade.name,gene+'_deletions']
    dels = [i for i in locs if i]
    missing = set()
    if len(dels):
        dels = dels[0].split(', ')
        for loc in dels:
            split = loc.split(':')
            asInt = [int(v) for v in split]
            if asInt[1] - asInt[0] >= 30:
                missing.update(range(asInt[0],asInt[1]+1))
    else:
        missing = {-1}
    return missing

def getNodeDels(subtree,meta,gene,states,leafs):
    for clade in subtree.find_clades():
        if clade in leafs:
            missing = extractDels(meta,clade,gene)
            states[clade.name] = missing
            return states[clade.name]
        else:
            flag = 0
            cluster = 0
            last = []
            current = set(range(0,29903+1))
            for child in clade:
                childState = getNodeDels(child,meta,gene,states,leafs)
                if len(childState) <= 1:
                    flag = 1
                else:
                    deleted = current.intersection(childState)
                    length,start,end = compareDels(deleted)
                    if length >= 30:
                        current = deleted
                        last.append(current)
                        cluster = 1
                    else:
                        flag = 1
                        if len(last):
                            current = last[-1]
                        else:
                            current = set(range(0,29903+1))
            if flag == 0:
                states[clade.name] = current
                return current
            elif cluster == 1:
                current.add(-1)
                states[clade.name] = current
            else:
                states[clade.name] = {-1}
            return {-1}

def getClustersDels(states,subtree):
    clusters = {}
    count = 1
    for clade in subtree.find_clades():
        if clade not in leafs:
            state = states[clade.name]
            if len(state) > 1:
                compared = state - {-1}
                clusters[count] = {}
                clusters[count]['strains'] = []
                for child in clade.find_clades():
                    if child in clade:
                        childState = states[child.name]
                        if -1 not in childState:
                            deleted = compared.intersection(childState)
                            length,start,end = compareDels(deleted)
                            if length >= 30:
                                compared = deleted
                                terms = child.get_terminals()
                                named = [leaf.name for leaf in terms]
                                clusters[count]['strains'].extend(named)
                                clusters[count]['deletion'] = str(start)+':'+str(end)
                count += 1
    return clusters

def getNodeStops(subtree,meta,gene,states,leafs):
    for clade in subtree.find_clades():
        if clade in leafs:
            if meta.loc[meta.strain==clade.name][gene+'_koType'].values[0] == 'earlyStop':
                states[clade.name] = [meta.loc[meta.strain==clade.name,gene+'_misStops'].values[0]]
            else:
                states[clade.name] = ['None']
            return states[clade.name]
        else:
            childStates = set()
            for child in clade:
                childStates.update(set(getNodeStops(child,meta,gene,states,leafs)))
            if len(childStates)==1:
                states[clade.name] = list(childStates)
                return list(childStates)
            else:
                states[clade.name] = list(childStates)
                return ['None']

def getClustersStops(states,subtree):
    clusters = {}
    count = 1
    for clade in subtree.find_clades():
        if clade not in leafs:
            state = states[clade.name]
            if len(state) > 1:
                for i in state:
                    if i != 'None':
                        clusters[count] = {}
                        clusters[count]['strains'] = []
                        for child in clade.find_clades():
                            if child in clade:
                                childState = states[child.name]
                                if len(childState) == 1:
                                    if i == childState[0]:
                                        terms = child.get_terminals()
                                        named = [leaf.name for leaf in terms]
                                        clusters[count]['strains'].extend(named)
                                        clusters[count]['mut']=i
                        count += 1
    return clusters

def getNodeFinch(subtree,meta,gene,states,leafs):
    for clade in subtree.find_clades():
        if clade in leafs:
            if meta.loc[meta.strain==clade.name][gene+'_koType'].values[0] == 'earlyStop':
                states[clade.name] = [meta.loc[meta.strain==clade.name,gene+'_misStops'].values[0]]
            else:
                states[clade.name] = ['None']
            return states[clade.name]
        else:
            childStates = []
            for child in clade:
                childStops  = set(getNodeFinch(child,meta,gene,states,leafs))
                childStates.append(childStops)

            intersect = set.intersection(*childStates)
            if not len(intersect):
                union = set.union(*childStates)
                nodeState = union
            else:
                nodeState = intersect
            states[clade.name] = list(nodeState)
            return list(nodeState)

def getClusterLeafs(value,states,subtree,named):
    for child in subtree.find_clades():
        if child in subtree:
            childState = states[child.name]
            if value in childState:
                if child in leafs:
                    named.append(child.name)
                else:
                    getClusterLeafs(value,states,child,named)


def getClustersFinch(states,subtree):
    clusters = {}
    count = 1
    for clade in subtree.find_clades():
        if clade not in leafs:
            state = states[clade.name]
            for i in state:
                if i != 'None':
                    clusters[count] = {}
                    clusters[count]['strains'] = []
                    clusters[count]['mut'] = i
                    getClusterLeafs(i,states,clade,clusters[count]['strains'])
                    count += 1
    return clusters

def trimClusters(clusters):
    clustered = set()
    rid = []
    for k,v in clusters.items():
        values = set(v['strains'])
        if values.issubset(clustered):
            rid.append(k)
        else:
            clustered.update(values)
    for k in rid:
        clusters.pop(k)
    return clusters

def saveClusters(stopC, delC,path,gene):
    clusts = []
    strains = []
    muts = []
    dels = []
    kotype = []
    for i,k in enumerate(stopC.keys()):
        samples = stopC[k]['strains']
        mut = stopC[k]['mut']
        length = len(samples)
        muts.extend([mut]*length)
        dels.extend(['None']*length)
        strains.extend(samples)
        clusts.extend([i+1]*length)
        kotype.extend(['earlyStop']*length)
    maxClust = i+1
    for j,l in enumerate(delC.keys()):
        samples = delC[l]['strains']
        delet = delC[l]['deletion']
        length = len(samples)
        dels.extend([delet]*length)
        muts.extend(['None']*length)
        strains.extend(samples)
        clusts.extend([maxClust+j+1]*length)
        kotype.extend(['bigDeletion']*length)
    df = pd.DataFrame({'strain':strains,'cluster':clusts,gene+'_misStops':muts,gene+'_deletions':dels,gene+'_koType':kotype})
    with open(path, 'w') as f:
        df.to_csv(f,sep='\t',index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--tree', type=str, required=True, help = 'path to nwk tree')
    parser.add_argument('--ko', type=str, required=True, help = 'path to tsv containing KOs')
    parser.add_argument('--outdir', type=str, required=True, help = 'path to output directory')
    args = parser.parse_args()

    genes = ['ORF1a','ORF1b','S','ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF9b']

    with open(args.tree,'r') as f:
        tree = Phylo.read(f, "newick")

    with open(args.ko,'r') as f:
        meta = pd.read_csv(f,sep='\t')

    leafs = tree.get_terminals()

    meta = replaceNaNs(meta,genes)
    names = tabulate_names(tree)

    for gene in genes:

        # Call all node deletions across internal nodes
        dels = {}
        getNodeDels(tree,meta,gene,dels,leafs)

        # Call deletion clusters
        delClusters = getClustersDels(dels,tree)

        # call all early stops across internal nodes
        stops = {}
        #getNodeStops(tree,meta,gene,stops,leafs)
        getNodeFinch(tree,meta,gene,stops,leafs)

        # call stop clusters
        #stopClusters = getClustersStops(stops,tree)
        stopClusters = getClustersFinch(stops,tree)
        trimClusters(stopClusters)

        # combine and write out
        outpath = args.outdir + 'clusters_' + gene + '.tsv'
        saveClusters(stopClusters,delClusters,outpath,gene)
