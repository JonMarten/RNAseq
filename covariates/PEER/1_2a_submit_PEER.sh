#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --mem 30G
#SBATCH --time=36:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/runPEER_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH --job-name=runPEER
#SBATCH -A PAUL-SL2-CPU

#!/bin/bash
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz 
source activate PEER

Rscript /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/covariates/PEER/1_2b_run_PEER.R