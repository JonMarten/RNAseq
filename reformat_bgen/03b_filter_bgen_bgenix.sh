#!/bin/bash
#SBATCH --job-name=filter_bgen_bgenix
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 12G
#SBATCH --time=12:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/update_bgen_bgenix_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load bgen

BGEN=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38.bgen
INCL=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/c${SLURM_ARRAY_TASK_ID}_b38_filter_snps.txt
OUT=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.bgen

# extract variants in lists
bgenix -g $BGEN -incl-rsids $INCL > $OUT