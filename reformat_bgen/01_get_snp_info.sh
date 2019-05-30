#!/bin/sh
# Get SNP info from master bgen files
#SBATCH -p medium
#SBATCH --job-name=snpstats
#SBATCH -t 24:0:0
module load gcc/5.2.0 
module load qctool2/rc4-6.8

qctool\
   -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
   -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
   -snp-stats\
   -osnp /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/snp_stats/impute_${SLURM_ARRAY_TASK_ID}_interval_snp_stats_unfiltered.txt