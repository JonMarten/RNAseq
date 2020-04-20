#!/bin/bash

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load plink-1.9-gcc-5.4.0-sm3ojoi

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes

plink\
 --bfile merged_cleaned_chrx\
 --keep /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/processing/rna_seq_phase1-2_affy_ids.txt\
 --make-bed\
 --out merged_cleaned_chrx_RNAseq_phase1-2