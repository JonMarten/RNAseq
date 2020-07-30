#!/bin/bash
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 20G
#SBATCH --job-name=plink_convert_b38
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/plink/logs/plink_convert_b38_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load ceuadmin/plink/2.0_09_09_18

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered

plink2 \
 --bgen impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered_no0_rnaSeqPhase1.bgen\
 --make-bed\
 --out plink/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered_no0_rnaSeqPhase1