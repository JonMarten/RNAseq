#!/bin/bash
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 20G
#SBATCH --job-name=plink_convert_b38
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/genotype_processing/plink_convert_b38_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load ceuadmin/plink/2.0_09_09_18

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes

plink2 \
 --bgen processing/b38_bgen/filtered/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered_no0_rnaSeqPhase1.bgen\
 --make-bed\
 --out INTERVAL_chr${SLURM_ARRAY_TASK_ID}_imputed_b38_filtered_rnaSeq