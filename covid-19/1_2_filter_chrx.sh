#!/bin/bash

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load plink-1.9-gcc-5.4.0-sm3ojoi

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes

cp /home/jm2294/rds/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data/INT_X_filt_merged.vcf.gz .
gunzip INT_X_filt_merged.vcf.gz

plink\
 --vcf INT_X_filt_merged.vcf\
 --chr 23\
 --keep /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/processing/rna_seq_phase1-2_affy_ids.txt\
 --make-bed\
 --out INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2