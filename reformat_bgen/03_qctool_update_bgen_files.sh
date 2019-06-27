#!/bin/sh
#SBATCH -p long
#SBATCH --job-name=update_bgen
#SBATCH -t 72:0:0
#SBATCH --output=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/update_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH -a 1-22
# submit with --dependency=afterok:1640947 to run after step 02b

# Remake BGEN files with new positions and duplicates removed

# get start time
start=$(date +%s.%N)

module load gcc/5.2.0 
module load qctool2/rc4-6.8

# Note that sample file can just be the chr22 file for any chromosome as this does not change.
# Update position and make new unique identifiers
qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -map-id-data /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/INTERVAL_chr${SLURM_ARRAY_TASK_ID}_b37_to_b38_map.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38.bgen

# use unique identifiers to retain only non-duplicated SNPs mapped to b38
qctool
 -g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38.bgen\
 -incl-snpids /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/c22_b38_filter_snps.txt\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.bgen
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"