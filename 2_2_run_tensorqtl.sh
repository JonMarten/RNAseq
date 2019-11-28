#!/bin/bash
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/logs/tensorqtl_test_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH -J TensorQTL_test
#SBATCH -A PAUL-SL3-GPU
#SBATCH --time=12:00:00

#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (Note probably this should not exceed the total number of GPUs in use.)
#SBATCH --ntasks=4
#! Specify the number of GPUs per node (between 1 and 4; must be 4 if nodes>1).
#! Note that the job submission script will enforce no more than 3 cpus per GPU.
#SBATCH --gres=gpu:4
#SBATCH -p pascal

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-gpu
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz 
module load cuda/9.2

source activate tensorQTL

DIR=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl

GPATH=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/plink/impute_22_interval_b38_filtered_no0_rnaSeqPhase1
PHEPATH=${DIR}/phenotypes/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed
OPATH=${DIR}/results/tqtl_c22_test
COVPATH=${DIR}/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_PC10_PEER20.txt

python3 \
 -m tensorqtl ${GPATH} ${PHEPATH} ${OPATH}\
 --covariates ${COVPATH}\
 --mode cis\
 --maf_threshold 0.005
