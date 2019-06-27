#!/bin/bash
#SBATCH --job-name=snpmap
#SBATCH --time=3:0:0
#SBATCH --partition=short 
#SBATCH --output=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/make_snp_map_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

module load R
Rscript /home/jm2294/projects/RNAseq/scripts/RNAseq/reformat_bgen/02_generate_snp_ids.R $SLURM_ARRAY_TASK_ID