#!/bin/bash
#SBATCH --job-name=snpmap
#SBATCH --time=3:0:0
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/make_snp_map_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH --mem 10G

module load R
Rscript /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/reformat_bgen/02_generate_snp_ids.R $SLURM_ARRAY_TASK_ID