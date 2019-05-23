# INTERVAL RNA seq eQTL analysis
This is a repository of scripts used in the analysis of the RNA seq data from the INTERVAL cohort. Currently this is a work in progress, but this README will be updated as the analysis progresses.

## Data
### Phenotype
Batches 1-4 of seq data have been downloaded from the Sanger server with globus, these are currently stored in `/home/jm2294/projects/RNAseq/globus`. CRAM files have not been downloaded. 
### Covariate
Covariate data has been provided for the whole cohort by the data management team; the current data relase is stored in `V:\INTERVAL\INTERVAL_MasterMerge\p1103\20190513\Data`. Phenotypes are stored in `INTERVALdata_13MAY2019.csv` and `INTERVALdata_P3_13MAY2019.csv`. `omicsMap.csv` and `omicsMap_P3.csv` map the randomly-generated 7 digit individual identifiers used by the data management team to the identifiers used on the various omic datasets, including the RNA-seq data. 

The INTERVAL study is split into phases which do not directly affect the RNA seq data but do affect how phenotype data is stored. **Each individual has only one timepoint at which RNA was measured.** This is preferentially two years after recruitment but may be later if no blood was available at this timepoint.

[This script](make_id_mapper_file.R) processes the omicsmap files to create a single file, `rna_id_mapper.csv`. This maps the phenotype '**identifier**' to the INT_RNA identifier used for the RNA seq files. It also lists '**batch**' (corresponding to RNA seq data release batch from the Sanger) and '**phase**', corresponding to the phase of the sample used for sequencing. This is obtained from the omicsmap file - the column containing the INT_RNA identifier specifies the phase. Where an id is included in multiple columns it is assigned one on the following order of preference: 24m > 48m > p3. Phase information is used to calculate the correct age for the RNA seq blood sample.

#### Calculating Ages
Ages are calculated using date of birth (assuming the 15th of the month) and the appointment date corresponding to blood sample, as obtained above. The plots below show how age varies with phase and batch, both absolute age and time between initial appointment and the blood sample used for RNA. Note that all batch 1 samples (the pilot study) are in the 24 month phase. Currently, phentoype data is only mapped to INT_RNA identifiers for batches 1 and 2.

![age plot](https://github.com/JonMarten/RNAseq/blob/master/bothplot.png?raw=true)


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
