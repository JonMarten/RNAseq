#!/bin/bash
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/logs/TensorQTL_pythonmodule_test_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH -J TensorQTL_pythonmodule_test
#SBATCH -A INOUYE-SL2-GPU
#SBATCH --time=36:00:00
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

python /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/2_2_tensorqtl_subset.py
