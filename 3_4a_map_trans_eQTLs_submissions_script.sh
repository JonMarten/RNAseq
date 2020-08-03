#!/bin/bash
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/tensorqtl_trans_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH -J trans_eQTL
#SBATCH -A INOUYE-SL2-GPU
#SBATCH --time=0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
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

python /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/3_4b_map_trans_eQTLs.py $CHR
