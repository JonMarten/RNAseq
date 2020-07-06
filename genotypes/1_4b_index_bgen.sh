#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --mem 10G
#SBATCH --time=12:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/index_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH --job-name=index_bgen
#SBATCH -A PAUL-SL3-CPU

#!/bin/bash
module load ceuadmin/bgenix/1.0.2

# Create bgi index file for bgen files
cd /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/processing/b38_bgen/filtered
CHR=$SLURM_ARRAY_TASK_ID
bgenix -g impute_${CHR}_interval_b38_filtered_no0_rnaSeqPhase1_2.bgen -index 