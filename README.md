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
It became apparent that duplicate SNP ids in the file were causing problems with the pipeline, so the current version prunes out any SNP with nonunique identifiers. These include SNPS without rsids coded as "." as well as triallelic SNPs and indels (See [this script](
### Pipeline
The Limix python pipeline as written by Marc Jan Bonder is currently being tested. Limix is installed on Cardio and is currently in the process of being set up on CSD3.  
