#!/bin/bash
#SBATCH -A PAUL-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --job-name=eqtl_merge_postprocessed_output
#SBATCH --time=12:0:0
#SBATCH --mem 100G
#SBATCH --output=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/logs/eqtl_merge_postprocessed_output_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load R

Rscript /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/scripts/RNAseq/04_merge_postprocessed_output.R $1