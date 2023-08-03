'''
This script takes a .fasta alignment and finds all truncations & deletions in
all genes. A gene is called a KO if it is missing 10 amino acids or has a
deletion 30+ bp long.
'''

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import re
import numpy as np

def getNsGaps(align):
    with open(align,'r') as f:
        records = list(SeqIO.parse(f, 'fasta'))
    ids = []
    ns = []
    gaps = []

    for record in records:
        seq = record.seq
        ids.append(record.id)
        n = [(match.start(), match.end()) for match in re.finditer('N+',str(seq))]
        ns.append(n)
        gap = [(match.start(), match.end()) for match in re.finditer('-+',str(seq))]
        gaps.append(gap)

    df = pd.DataFrame({'strain': ids, 'Ns': ns, 'gaps': gaps})

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--align', type=str, required=True, help = 'path to alignment file as fasta')
    parser.add_argument('--ref', type=str, required=True, help = 'genbank file')
    parser.add_argument('--amplicon', default=False, action='store_true', help = 'flag to call known amplicon dropout sites as KO')
    parser.add_argument('--output', type=str, required=True, help = 'path to output TSV')
    args = parser.parse_args()

    with open(args.align,'r') as f:
        records = list(SeqIO.parse(f, 'fasta'))

    ref = SeqIO.read(args.ref,'gb')

    df = getNsGaps(args.align)

    for feature in ref.features:
        if feature.type == "CDS":
            locs = []
            glocs = []
            missing = []
            bigDels = []
            bigGaps = []
            plengths = []
            aa_stops = []
            aa_stoploss = []
            fshifts = []
            missenses = []
            kos = []
            types = []

            start = feature.location.start
            end = feature.location.end
            name = feature.qualifiers['gene'][0]
            refgene = ref.seq[start:end]
            refprotein = refgene.translate(to_stop=True)
            refplength = len(refprotein)

            for nlocs, gaplocs, record in zip(df['Ns'],df['gaps'],records):
                loc = []
                gloc = []
                fshift = []
                gone = 0
                for i in range(len(nlocs)):
                    nstart,nend = nlocs[i]
                    if nstart != 0 and nend != 29903:
                        if nstart >= start and nstart <= end:
                            loc.append((nstart,nend))
                            if nend > end:
                                gone += (end - nstart + 1)
                            else:
                                gone += (nend - nstart + 1)
                        elif nstart < start and nend>=start:
                            loc.append((nstart,nend))
                            gone += (nend - start + 1)
                for i in range(len(gaplocs)):
                    nstart,nend = gaplocs[i]
                    gapLength = nend-nstart
                    if nstart != 0 and nend != 29903:
                        if nstart >= start and nstart <= end:
                            loc.append((nstart,nend))
                            gloc.append((nstart,nend))
                            if nend > end:
                                gone += (end - nstart + 1)
                            else:
                                gone += (nend - nstart + 1)
                        #    if (gapLength)%3 != 0:
                        #        fshift.append(str(nstart)+'-'+str(nend)+':'+refgene[nstart:nend]+'>'+''.join(['-'*gapLength]))
                        elif nstart < start and nend>=start:
                            loc.append((nstart,nend))
                            gloc.append((nstart,nend))
                            gone += (nend - nstart + 1)
                        if (gapLength)%3 != 0:
                            fshift.append(str(nstart)+'-'+str(nend)+':'+refgene[nstart:nend]+'>'+''.join(['-'*gapLength]))

                if len(loc):
                    if args.amplicon == False:
                        if record.seq[27807-1]=='T':
                            drop = []
                            for i,(cstart,cend) in enumerate(loc):
                                if cstart >= 27808 and cstart < 27854:
                                    drop.append(i)
                            if len(drop)>0:
                                for idx in reversed(drop):
                                    del loc[idx]

                    sizes = [cend - cstart for cstart,cend in loc]
                    if len(loc):
                        biggest = max(sizes)
                    else:
                        biggest=0
                else:
                    biggest = 0

                if len(gloc):
                    gsizes = [cend - cstart for cstart,cend in gloc]
                    biggestg = max(gsizes)
                else:
                    biggestg = 0

                locs.append(sorted(loc))
                glocs.append(gloc)
                missing.append(gone)
                bigDels.append(biggest)
                bigGaps.append(biggestg)
                fshifts.append(fshift)

                missense = ''
                seq = record.seq
                if name == 'ORF1a':
                    if len(gaplocs):
                        gstart,gend = gaplocs[0]
                        if gstart == 0:
                            seq = Seq('N'*(gend+1))+record.seq[gend+1:]

                if name == 'N':
                    if len(gaplocs):
                        gstart,gend = gaplocs[-1]
                        if gend == 29903:
                            seq = record.seq[:gstart-1]+Seq('N'*(gend-gstart+1))
                genePotent = seq[start:]
                no_gaps = Seq(''.join([bp for bp in genePotent if bp != '-']))
                protein = no_gaps.translate(to_stop=True)
                length = len(protein)
                plengths.append(length)
                if length < refplength:
                    bpStart = (length*3)
                    bpEnd = length*3+3
                    codon = no_gaps[bpStart:bpEnd]
                    codonref = refgene[bpStart:bpEnd]
                    if codon != codonref:
                        missense =str(start+bpStart)+'-'+str(start+bpEnd-1)+':'+codonref+'>'+codon
                missenses.append(missense)

            df[name+'_deletions'] = locs
            df[name+'_gaps'] = glocs
            df[name+'_nDeleted'] = missing
            df[name+'_maxDeletion'] = bigDels
            df[name+'_maxGap'] = bigGaps
            df[name+'_frameShifts'] = fshifts
            df[name+'_proteinLength'] = plengths
            df[name+'_misStops'] = missenses
            df[name+'_ko'] = np.where((df[name+'_maxDeletion']>=30)|((refplength - df[name+'_proteinLength'])>= 10),"Yes","No")
            df[name+'_koType'] = np.where(df[name+'_maxDeletion']>=30,'BigDeletion','')
            df[name+'_koType'] = np.where((refplength - df[name+'_proteinLength'])>= 10 ,'earlyStop',df[name+'_koType'])
            df[name+'_koType'] = np.where(df[name+'_maxDeletion']>=30,'BigDeletion',df[name+'_koType'])
    df.to_csv(args.output, sep = '\t', index=False)
