#!/bin/bash
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/logs/tensorqtl_cond_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH -J TensorQTL_test
#SBATCH -A PAUL-SL3-GPU
#SBATCH --time=12:00:00

#! How many whole nodes should be allocated?
#SBATCH --nodes=2
#! How many (MPI) tasks will there be in total? (Note probably this should not exceed the total number of GPUs in use.)
#SBATCH --ntasks=4
#! Specify the number of GPUs per node (between 1 and 4; must be 4 if nodes>1).
#! Note that the job submission script will enforce no more than 3 cpus per GPU.
#SBATCH --gres=gpu:4
#SBATCH -p pascal

#! ############################################################
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-gpu
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz 
module load cuda/9.2

source activate tensorQTL

DIR=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl

GPATH=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.005
PHEPATH=${DIR}/phenotypes/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz
OPATH=${DIR}/results/tensorqtl_cis_MAF0.005
COVPATH=${DIR}/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_PC10_PEER20.txt

python\
 -m tensorqtl ${GPATH} ${PHEPATH} ${OPATH}_cis_independent\
 --covariates ${COVPATH}\
 --cis_output /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results/tensorqtl_cis_MAF0.005_cisnominal.tensorQTL.cis_nominal_manualPython.csv\
 --mode cis_independent
