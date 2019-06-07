# INTERVAL RNA seq eQTL analysis
This is a repository of scripts used in the analysis of the RNA seq data from the INTERVAL cohort. Currently this is a work in progress, but this README will be updated as the analysis progresses.

## Data
### Phenotype
Batches 1-4 of seq data have been downloaded from the Sanger server with globus, these are currently stored in `/home/jm2294/projects/RNAseq/globus`. CRAM files have not been downloaded. 
### Covariates
Covariate data is stored in a QC-ready file `INTERVAL_RNA_batch1_5_covariates_release31MAY2019.csv`. This file includes technical covariates from the RNA sequencing run as well as age and sysmex data from the same timepoint as the blood sample used for RNA seq. Derivation of this file is detailed [here](https://github.com/JonMarten/RNAseq/blob/master/covariates/README.md#phenotypes-in-the-interval-study).

### Genotype
The 'master' bgen files are stored in `/scratch/curated_genetic_data/interval/imputed/impute_[CHR]_interval.bgen`, symlinked in `/home/jm2294/GENETIC_DATA/INTERVAL/master`
These are unsuitable for limix analysis as the pipeline requires bgen v1.2, which contains sample identifiers embedded in the file. 
To test the limix pipeline, a new version of the bgen files has been created with qctool v2 (see [this script](make_test_bgen.sh) for details).
In the first instance, the file was pruned to just chr22:23500000-24500000, and limited to the 188 individuals in batch 1. 
It became apparent that duplicate SNP ids in the file were causing problems with the pipeline, so the current version prunes out any SNP with nonunique identifiers. These include SNPS without rsids coded as "." as well as triallelic SNPs and indels ([This script](get_duplicate_rsids.R) documents how these were identified). 
## Pipeline
The Limix python pipeline as written by Marc Jan Bonder is currently being tested. Limix is installed on Cardio and is currently in the process of being set up on CSD3. It currently works on phased test data provided with the code, but not on our own data.
### Installation
The 'HipSci Pipeline' is currently only available to download from Marc's Google Drive, please email me for the link.

The file limix_install.txt is mostly adapted from a file contained within the Google Drive. I have detailed the commands I used to get the pipeline running below:
1. Download `hipsci_pipeline.zip` from Google Drive and unzip to cardio. This is currently stored in `/home/jm2294/projects/RNAseq/hipsci_pipeline`.
2. Set up a new **conda** environment to ensure the correct versions of each package are installed with `conda create -n limix_qtl python=3 anaconda`
3. Activate the environment with `source activate limix_qtl`
4. Install dependencies: `conda install -c anaconda pytest pytables`
5. Install limix: `pip install limix==2.0.3`. **Note that version 2.0.3 of limix is required for the pipeline to work, as limix 3.0 outputs different file structures.** This will also install the required versions of bgen-reader and other dependencies.
6. If using snakemake, install with `pip install snakemake==4.5.0` (the current version of snakemake cannot be installed due to problems with the `datrie` depedent package. 

Limix can now be called directly from the command line, but Marc's pipeline is implemented in a series of python scripts. QTLs are called with `hipsci_pipeline/limix_QTL_pipeline\run_QTL_analysis.py`. His implementation uses a **snakemake** file to manage the workflow and chunk the file into manageable pieces. This is currently being investigated.
