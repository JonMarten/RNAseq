#!/bin/bash
#SBATCH -A PAUL-SL2-CPU
#SBATCH -p skylake
#SBATCH --mem 150G
#SBATCH --job-name=plink_merge_1_22
#SBATCH --time=36:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/plink_merge_chromosomes_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load plink-1.9-gcc-5.4.0-sm3ojoi

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes

plink\
 --merge-list /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/genotypes/plink_autosomes_file_list.txt\
 --make-bed\
 --memory 191000\
 --out INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_AllAutosomes
