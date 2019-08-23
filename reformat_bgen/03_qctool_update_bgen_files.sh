#!/bin/bash
#SBATCH --job-name=update_bgen
#SBATCH -A PAUL-SL2-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 12G
#SBATCH --time=36:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/update_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

# Remake BGEN files with new positions and duplicates removed

# get start time
start=$(date +%s.%N)

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load qctool

# Update position and add new unique identifiers
qctool\
 -g /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
 -s /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples\
 -map-id-data /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr${SLURM_ARRAY_TASK_ID}_b37_to_b38_map.txt\
 -og /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_no0_tinytest.bgen

end=$(date +%s.%N)

runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"