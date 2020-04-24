#!/bin/bash
#SBATCH -A PAUL-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 50G
#SBATCH --job-name=filter_chrx
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes/plink_filter_log_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load plink-1.9-gcc-5.4.0-sm3ojoi

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes

plink\
 --vcf /home/jm2294/rds/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data/INT_X_filt_merged.vcf.gz\
 --chr 23\
 --keep /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/processing/rna_seq_phase1-2_affy_ids.txt\
 --make-bed\
 --out INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2