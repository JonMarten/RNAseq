#!/bin/sh
#SBATCH -p long
#SBATCH --job-name=vcf_make
#SBATCH -t 36:0:0

# get start time
start=$(date +%s.%N)

module load gcc/5.2.0 
module load qctool2/rc4-6.8

# Note that sample file can just be the chr22 file for any chromosome as this does not change.
qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.bgen\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.vcf
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"
