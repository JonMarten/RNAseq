module load gcc/5.2.0 
module load qctool2/rc4-6.8

qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_22_interval_b38_filtered.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -excl-rsids /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/c22_b38_filter_snps_by_rsids.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_22_interval_b38_filtered_on_rsids.bgen