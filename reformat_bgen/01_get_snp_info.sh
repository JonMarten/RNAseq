#!/bin/sh
# Get SNP info from master bgen files
module load gcc/5.2.0 
module load qctool2/rc4-6.8

for CHR in 1 .. 22
do
  qctool\
   -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${CHR}_interval.bgen\
   -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
   -snp-stats\
   -osnp /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/snp_stats/impute_${CHR}_interval_snp_stats_unfiltered.txt
done