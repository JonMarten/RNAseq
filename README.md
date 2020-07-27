# INTERVAL RNA-seq eQTL analysis

## Scripts
This is a repository of scripts used in the analysis of the RNA seq data from the INTERVAL cohort. Current generation scripts are stored in the root folder, with the [initial Limix pipeline](01_limix_pipeline) and the [tensorQTL phase I analysis](02_tensorqtl_phase_1) scripts stored in subfolders for reference purposes.

Where possible, I've tried to name files in a way that makes it obvious what order they have to be run in.

Unless otherwise specified, all file paths are given relative to `/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/`.

## Data

### Phenotype

#### Raw data
RNA-seq data is downloaded from the Sanger HGI service Globus server. This is accessed from https://app.globus.org/file-manager and requires an endpoint to be set up on CSD3 to transfer files. The commands required to configure this are stored in [globus_config_for_csd3.txt](globus_config_for_csd3.txt).

Initial Phase I data is stored in `/rds/project/jmmh2/rds-jmmh2-pre_qc_data/interval/rna_seq/raw_data/globus` in subfolders by batch. Within each subfolder, `.cram` files are stored tar files in the `data` folders and processed gene counts are stored in the `results-study5591-*` folders.

Phase I was recalled together with batches 9-12 of Phase II to bring INTERVAL into closer alignment with the BioAid and GAINS studies. This latest data release is stored on the globus server, with some files downloaded to `globus_phase2_recalled`. Not all files have been copied,  as many of these files are symlinks on the Sanger side and are ignored by the globus sync functions. This is something that needs to be worked out with Guillaume Noell.

#### Reformatted data
Gene counts from globus have been filtered and TMM-normalised by Artika Nath - a version of the script used for this is stored [here](artika_TMM_normalisation.R) but she will have the most up-to-date version. The file used as input for phenotype generation is the same one used for peer and is consequently stored in `analysis/04_phase2_full_analysis/peer_factors/peer_InputFiles/GeneExpr_PEER_TmmInvRankNormalised_swapsSwapped_mismatchRemoved.csv`. Phenotype input files for TensorQTL are stored in `analysis/04_phase2_full_analysis/phenotypes` and are generated with [this script](3_1_make_tensorQTL_input_phase2.R).

#### Annotation
Genomic positions of genes are obtained from Ensembl. Annotation files are stored in `analysis/04_phase2_full_analysis/annotation_file`, with the raw BioMart output stored in `19_9_5_mart_export.txt` and the reformatted annotation file in `Feature_Annotation_Ensembl_gene_ids_autosomesPlusChrX_b38.txt`. This is a throwback to the Limix pipeline where annotation was supplied as a separate file, but this information is now stored within the TensorQTL .bed input. However, the phenotype generation scriptes still require the annotation file to be in the Limix format. The script to reformat this is stored [here](3_0_make_annotation_file_autosomes_plus_x.R).

### Covariates
Scripts for covariate handing are stored in [their own subfolder](covariates). The master file containing all covariates is `analysis/04_phase2_full_analysis/covariates/processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_07_01.csv`. The analysis ready covariate file is stored in `analysis/04_phase2_full_analysis/covariates`.

### PEER factors
PEER factors were also generated by Artika. The scripts for doing this are duplicated [here](covariates/PEER). 20 factors were generated as more than this contribute little to explaining variance in the RNA-seq data. Biological covariates explicitly adjusted for are detailed in `analysis/04_phase2_full_analysis/peer_factors/peer_InputFiles/Covariates_for_PEER_mismatchRemoved.csv`. Additional technical covariates are included in the regression model in TensorQTL. The rationale behind this was to explicitly correct for those covariates where there might be an interest in partitioning variance, e.g. the cell count traits. 

The Sysmex traits included as covariates were all the traits used in the [Astle et al Cell GWAS](https://pubmed.ncbi.nlm.nih.gov/27863252/) minus the cell counts that were also represented as a percentage (eg Lymphocyte percentage was retained while Lymphocyte count was not. This was because the former was used to derive the latter).

### Genotype
TensorQTL takes plink genotype format as input. This does mean imputed genotypes are hard-called, but I consider the trade off in speed to be worth it. 

The scripts in the [genotypes](genotypes) subfolder describe the filters applied. Genetic PCs are calculated using [this script](genotypes/1_7_plink_get_PCs.sh) which is based on Ben's code from the SomaLogic analysis.

ChrX files were created as part of the [COVID-19 subproject](covid-19). The files `INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2_b38_rsids_deduplicated_MAF0.005.*` were copied to `analysis/04_phase2_full_analysis/genotypes` and renamed `INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr23.*`


## Pipeline
Initial cis-eQTL mapping was performed using the [Limix Pipeline](01_limix_pipeline). Experimental trans-eQTL mapping trialed in Limix, but this was switched to TensorQTL, which has been used for all subsequent analyses.

Current scripts are as follows:
* [3_0_make_annotation_file_autosomes_plus_x.R](3_0_make_annotation_file_autosomes_plus_x.R): Converts BioMart annotation file into the right format for the next step.
* [3_1_make_tensorQTL_input_phase2.R](3_1_make_tensorQTL_input_phase2.R): Output phenotype `.bed` files for use in TensorQTL.
* [3_2_index_bed.sh](3_2_index_bed.sh): Compress and index `.bed` files from previous step.
* [3_3a_map_cis_eQTLs_submissions_script.sh](3_3a_map_cis_eQTLs_submissions_script.sh) and [3_3b_map_cis_eQTLs.py](3_3b_map_cis_eQTLs.py): The python script that runs the cis-eQTL mapping in TensorQTL, and the shell script for submission to CSD3.

## Files
All files on CSD3 are currently located in the project folder, `/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/` (i.e. not in the GWASqc folder, since they're not yet finalised).

