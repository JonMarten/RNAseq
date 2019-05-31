#!/bin/sh
# Get SNP info from master bgen files
#SBATCH -p long
#SBATCH --job-name=update_bgen
#SBATCH -t 36:0:0

# get start time
start=$(date +%s.%N)

module load gcc/5.2.0 
module load qctool2/rc4-6.8

qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -map-id-data\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38.bgen
 
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"