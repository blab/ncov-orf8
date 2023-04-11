#!/bin/bash
#SBATCH --job-name=refine

  augur refine \
    --tree WA_20K_tree_raw.nwk \
    --alignment WA_20K.fasta\
    --metadata WA_20K_metadata.tsv \
    --root Wuhan/Hu-1/2019 \
    --coalescent opt \
    --divergence-unit mutations \
    --keep-polytomies \
    --output-tree WA_20K_tree.nwk \
    --output-node-data WA_20K_branch_lengths.json \
    --timetree \
    --clock-rate 0.0008 \
    --clock-std-dev 0.0004 \
    --date-inference marginal \
    --date-confidence \
    --no-covariance \
