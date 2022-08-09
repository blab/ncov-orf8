'''
This script takes a .fasta alignment and finds genes with either deletions or truncations.
It defaults to look for knockouts in SARS-CoV-2 ORF8 according to alignment to 
[SARS2 ref sequence](https://github.com/nextstrain/ncov/blob/master/defaults/reference_seq.gb).

Must provide alignment file, --align, and output file, --output. If want to look
for KO in other genes, provide the start, --start, and end, --end of the gene.
'''

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="check sequences for anomalies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--align', type=str, required=True, help = 'path to alignment file as fasta')
    parser.add_argument('--start', type=int, default = 27894, help = 'starting base of gene')
    parser.add_argument('--end', type=int, default = 28259, help = 'end base of gene')
    parser.add_argument('--output', type=str, required=True, help = 'path to output TSV')
    args = parser.parse_args()

    ids = []
    ns = []
    gaps = []
    protein_lengths = []

    records = list(SeqIO.parse(args.align, 'fasta'))

    for record in records:
        gene = record.seq[args.start-1:args.end]
        ids.append(record.id)
        ns.append(gene.count("N"))
        gaps.append(gene.count("-"))
        no_gaps = Seq(''.join([bp for bp in gene if bp != '-']))
        protein = no_gaps.translate(to_stop=True)
        protein_lengths.append(len(protein))

    df = pd.DataFrame({'strain': ids, 'Ns': ns, 'gap': gaps, 'protein_length': protein_lengths})

    df.to_csv(args.output, sep = '\t', index=False)
