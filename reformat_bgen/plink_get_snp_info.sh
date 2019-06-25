#!/bin/bash
module load plink/2.0_09_09_18

cd /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover

for i in {1..22}
do
  plink2 --sample /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample --bgen /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_${i}_interval.bgen --out /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/plink/impute_${i}_interval_b38_filtered_on_rsids --write-snplist
done