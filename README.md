# INTERVAL RNA seq eQTL analysis
This is a repository of scripts used in the analysis of the RNA seq data from the INTERVAL cohort. Currently this is a work in progress, but this README will be updated as the analysis progresses.
## Current state of analysis
### Data
#### Phenotype
Batches 1-3 of seq data have been downloaded from the Sanger server with globus, these are currently stored in `/home/jm2294/projects/RNAseq/globus`. CRAM files have not been downloaded. 
#### Covariate
Covariate data has been requested for the whole cohort from the data management team. Identifier mapping is currently available for batches 1 & 2 (see [the appropriate R script](make_id_mapper_file.R) for creation of the mapping file).	
#### Genotype
The 'master' bgen files are stored in `/scratch/curated_genetic_data/interval/imputed/impute_[CHR]_interval.bgen`, symlinked in `/home/jm2294/GENETIC_DATA/INTERVAL/master`
These are unsuitable for limix analysis as the pipeline requires bgen v1.2, which contains sample identifiers embedded in the file. 
To test the limix pipeline, a new version of the bgen files has been created with qctool v2 (see [this script](make_test_bgen.sh) for details).
In the first instance, the file was pruned to just chr22:23500000-24500000, and limited to the 188 individuals in batch 1. 
It became apparent that duplicate SNP ids in the file were causing problems with the pipeline, so the current version prunes out any SNP with nonunique identifiers. These include SNPS without rsids coded as "." as well as triallelic SNPs and indels ([This script](get_duplicate_rsids.R) documents how these were identified). 
### Pipeline
The Limix python pipeline as written by Marc Jan Bonder is currently being tested. Limix is installed on Cardio and is currently in the process of being set up on CSD3. It currently works on phased test data provided with the code, but not on our own data.
#### Installation
The 'HipSci Pipeline' is currently only available to download from Marc's Google Drive, please email me for the link.

The file limix_install.txt is mostly adapted from a file contained within the Google Drive. I have detailed the commands I used to get the pipeline running below:
1. Download `hipsci_pipeline.zip` from Google Drive and unzip to cardio. This is currently stored in `/home/jm2294/projects/RNAseq/hipsci_pipeline`.
2. Set up a new **conda** environment to ensure the correct versions of each package are installed with `conda create -n limix_qtl python=3 anaconda`
3. Activate the environment with `source activate limix_qtl`
4. Install dependencies: `conda install -c anaconda pytest pytables`
5. Install limix: `pip install limix==2.0.3`. **Note that version 2.0.3 of limix is required for the pipeline to work, as limix 3.0 outputs different file structures.** This will also install the required versions of bgen-reader and other dependencies.

Limix can now be called directly from the command line, but Marc's pipeline is implemented by calling  `hipsci_pipeline/limix_QTL_pipeline\run_QTL_analysis.py`
