# INTERVAL RNA-seq eQTL analysis
## Scripts
This is a repository of scripts used in the analysis of the RNA seq data from the INTERVAL cohort. Current generation scripts are stored in the root folder, with the [initial Limix pipeline](01_limix_pipeline) and the [tensorQTL phase I analysis](02_tensorqtl_phase_1) scripts stored in subfolders for reference purposes.

## Data
### Phenotype
#### Raw data
RNA-seq data is downloaded from the Sanger HGI service Globus server. This is accessed from https://app.globus.org/file-manager and requires an endpoint to be set up on CSD3 to transfer files. The commands required to configure this are stored in [globus_config_for_csd3.txt](globus_config_for_csd3.txt).

Initial Phase I data is stored in `/rds/project/jmmh2/rds-jmmh2-pre_qc_data/interval/rna_seq/raw_data/globus` in subfolders by batch. Within each subfolder, `.cram` files are stored tar files in the `data` folders and processed gene counts are stored in the `results-study5591-*` folders.

Phase I was recalled together with batches 9-12 of Phase II to bring INTERVAL into closer alignment with the BioAid and GAINS studies. This latest data release is stored on the globus server, with some files downloaded to `/rds/rds-jmmh2-projects/interval_rna_seq/globus_phase2_recalled`.

### Annotation
Genomic positions of genes are obtained from Ensembl. 58,394 features are mapped in the original RNA seq counts. Of these, 92 are listed with "ERCC" identifiers, which I believe are spike-ins. 382 features are retired in the current release of Ensembl and are not included in the annotation file. In the final annotation file, chromosome and position information is incldued for 57,861 Ensembl genes.

### Covariates
Covariate data is stored in a QC-ready file `INTERVAL_RNA_batch1_5_covariates_release31MAY2019.csv`. This file includes technical covariates from the RNA sequencing run as well as age and sysmex data from the same timepoint as the blood sample used for RNA seq. Derivation of this file is detailed [here](https://github.com/JonMarten/RNAseq/blob/master/covariates/README.md#phenotypes-in-the-interval-study).

### Genotype
The 'master' bgen files are stored in `/scratch/curated_genetic_data/interval/imputed/impute_[CHR]_interval.bgen`, symlinked in `/home/jm2294/GENETIC_DATA/INTERVAL/master`
These are unsuitable for limix analysis as the pipeline requires bgen v1.2, which contains sample identifiers embedded in the file. 
To test the limix pipeline, a new version of the bgen files has been created with qctool v2 (see [this script](make_test_bgen.sh) for details).
In the first instance, the file was pruned to just chr22:23500000-24500000, and limited to the 188 individuals in batch 1. 
It became apparent that duplicate SNP ids in the file were causing problems with the pipeline, so the current version prunes out any SNP with nonunique identifiers. These include SNPS without rsids coded as "." as well as triallelic SNPs and indels ([This script](get_duplicate_rsids.R) documents how these were identified). 

Filter of MAF > 0.2% (roughly corresponding to MAC > 10 in the full sample of ~2,500) was applied before running the analysis. No imputation quality filter was applied pre-analysis as this is not implemented in the pipeline and makes little difference on the number of SNPs after filtering on MAF. A post-GWAS filter will be applied.

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
