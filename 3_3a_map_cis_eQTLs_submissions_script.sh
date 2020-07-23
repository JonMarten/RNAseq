#!/bin/bash
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/tensorqtl_cis_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH -J cis_eQTL
#SBATCH -A PAUL-SL3-GPU
#SBATCH --time=12:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --gres=gpu:4
#SBATCH -p pascal

#! ############################################################
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-gpu
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz 
module load cuda/9.2
module load r-3.6.0-gcc-5.4.0-bzuuksv

source activate tensorQTL

CHR=$SLURM_ARRAY_TASK_ID 

python /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/3_3b_map_cis_eQTLs.py $CHR
