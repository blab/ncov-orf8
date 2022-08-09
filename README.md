## Configuring environment
This code requires having biopython & pandas installed. If these are not already available, you can configure the appropriate conda environment from `ko.yaml`.
```
conda env create -f ko.yaml
conda activate knockout
```

## Running the code
`find_orf8_ko.py` takes a .fasta alignment and finds genes with either deletions or truncations.
The script defaults to look for knockouts in SARS-CoV-2 ORF8 according to alignment to 
[SARS2 ref sequence](https://github.com/nextstrain/ncov/blob/master/defaults/reference_seq.gb).
To run the script, please provide alignment file (--align) and output file (--output). If want to look
for KOs in other genes, provide the start (--start) and end (--end) of the gene.
```
python code/find_orf8_ko.py --align [FASTA alignment] --output [Path to output file]
```
