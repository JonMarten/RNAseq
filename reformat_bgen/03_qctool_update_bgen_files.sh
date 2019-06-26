#!/bin/sh
#SBATCH -p long
#SBATCH --job-name=update_bgen
#SBATCH -t 36:0:0

# Remake BGEN files with new positions and duplicates removed

# get start time
start=$(date +%s.%N)

module load gcc/5.2.0 
module load qctool2/rc4-6.8

# Note that sample file can just be the chr22 file for any chromosome as this does not change.
qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -excl-positions /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/plink/impute_${SLURM_ARRAY_TASK_ID}_duplicate_positions.txt\
 -map-id-data /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/INTERVAL_chr${SLURM_ARRAY_TASK_ID}_b37_to_b38_map.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38.bgen

end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"