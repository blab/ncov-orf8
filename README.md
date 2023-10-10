# Positive selection underlies repeated knockout of ORF8 in SARS-CoV-2 evolution

**Cassia Wagner**<sup>1,2,*</sup>, **Kathryn E. Kistler** <sup>2,3</sup>, **Garrett A. Perchetti** <sup>4</sup>, **Noah Baker** <sup>4</sup>, **Lauren A. Frisbie** <sup>5</sup>, **Laura Marcela Torres** <sup>5</sup>, **Frank Aragona** <sup>5</sup>, **Cory Yun** <sup>5</sup>, **Marlin Figgins** <sup>2,6</sup>, **Alexander L. Greninger** <sup>2,4</sup>,   **Alex Cox** <sup>5</sup>, **Hanna N. Oltean** <sup>5</sup>, **Pavitra Roychoudhury** <sup>2,4</sup>, **Trevor Bedford** <sup>1,2,3</sup> <br />

<sup>1</sup> *Department of Genome Sciences, University of Washington, Seattle, WA, USA;* <br />
<sup>2</sup> *Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA;* <br />
<sup>3</sup> *Howard Hughes Medical Institute, Seattle, WA, USA;* <br />
<sup>4</sup> *Department of Laboratory Medicine and Pathology, University of Washington, Seattle, Washington, USA;* <br />
<sup>5</sup> *Washington State Department of Health, Shoreline, Washington, USA;* <br />
<sup>6</sup> *Department of Applied Mathematics, University of Washington, Seattle, Washington, USA.* <br />
<sup>*</sup> *Corresponding author: cassiasw@uw.edu* <br />

Knockout of the ORF8 protein has repeatedly spread through the global viral population during SARS-CoV-2 evolution. Here we use both regional and global pathogen sequencing to explore the selection pressures underlying its loss. In Washington State, we identified transmission clusters with ORF8 knockout throughout SARS-CoV-2 evolution, not just on novel, high fitness viral backbones. Indeed, ORF8 is truncated more frequently and knockouts circulate for longer than for any other gene. Using a global phylogeny, we find evidence of positive selection to explain this phenomenon: nonsense mutations resulting in shortened protein products occur more frequently and are associated with faster clade growth rates than synonymous mutations in ORF8. Loss of ORF8 is also associated with reduced clinical severity, highlighting the diverse clinical impacts of SARS-CoV-2 evolution.

---

#### Structure of this repository
This repository includes the code for the analyses and figures for the above manuscript.
Clinical data from Washington State Disease Reporting System is not included as this data is derived from confidential medical records.
GISAID metadata and sequenced used in the analysis may be accessed at [gisaid.org/EPI_SET_230921by](https://www.epicov.org/epi3/epi_set/EPI_SET_230921by?main=true).  
The SARS-CoV-2 UShER phylogeny is available from [UShER](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2023/05/01).
- `code` contains the scripts for all analyses.
- `data` contains simulated clinical data for use in severity analysis & mutation annotations from [Obermeyer et al](https://www.science.org/doi/full/10.1126/science.abm1208). Please access the GISAID sequences and SARS-CoV-2 UShER phylogeny using the above links.
- `nextstrain_build` contains the identified clusters and the configurations for building the nextstrain trees to identify transmission clusters of gene knockouts in Washington State.
- `envs` contains the conda config files for python code & notebooks and for matUtils.
- `figs` contains the figures for the manuscript.
- `notebooks` contains jupyter notebooks for plotting results and initial analyses.
- `params` includes the SARS-CoV-2 reference genomes used in analyses & the config file for snakemake pipeline.
- `usher` contains results from analyses using the usher phylogeny.

#### Setup & installation
Use [mamba](https://anaconda.org/conda-forge/mamba) to quickly (~5 min) install matUtils & python notebooks environments.
The environment for python scripts & notebooks can be set up & activated using:
```
# Install
mamba env create -f envs/orf8ko.yaml

# Activate
mamba activate orf8ko
```
The environment for matUtils can be set up & activated using:
```
# Install
mamba env create -f envs/usher-env.yaml

# Activate
mamba activate usher-env
```

Rscripts were run in RStudio using R version 4.1.2. The R environment dependencies are listed in `envs/renv.lock`. To use this environment:
```
# Install renv
install.packages("renv")

# Create and activate renvironment
renv::restore(lockfile = 'envs/renv.lock')
```
This process should take a few minutes.

#### Running the analyses
- Run `code/find_ko.py` on .fasta alignment of WA sequences to call potential gene knockouts. See above to access sequences and metadata from GISAID.
- Build and call transmission clusters using `nextstrain_build`
- Calculate dN/dS using the snakemake workflow: `code/dNdS_snakefile`. See above to download the UShER tree for this analysis.
- Call mutation clusters from UShER tree using `code/getMutationClusters.py`
- Model cluster growth rates using: `code/clusterSize_regression.R`
- `code/combineClinicalData.R` is used to generate the dataframe for clinical analysis.
- Use `code/Fig4.R` to run the clinical severity analysis. Although we cannot share clinical data to protect patient privacy, we have provided `data/clinical_example.tsv` as a demo dataset.
