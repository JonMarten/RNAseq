#!/bin/bash
module load plink/2.0_09_09_18

cd /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover

plink2 \
 --sample /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 --bgen /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
 --out /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/plink/impute_${SLURM_ARRAY_TASK_ID}_b37_unfiltered