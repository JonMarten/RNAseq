
module load gcc/5.2.0 
module load qctool2/rc4-6.8
qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -incl-samples /home/jm2294/projects/RNAseq/test_run/batch1_ids.txt\
 -incl-range 22:23500000-24500000\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile.bgen
 
 qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -incl-samples /home/jm2294/projects/RNAseq/test_run/batch1_ids.txt\
 -incl-range 22:23500000-24500000\
 -excl-snpids /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/filter_snps.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids.bgen
 
 # Filter out SNPs with no unique identifier ## this loses sample info, have to remake from scratch
 #qctool\
 #-g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile.bgen\
 #-excl-snpids filter_snps.txt\
 #-og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids.bgen
 
 
# check snp lists
qctool -g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile.bgen -snp-stats -osnp test_bgen_snp_stats.txt
qctool -g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids.bgen -snp-stats -osnp test_bgen_snp_stats_unique.txt

# Make new test bgen for whole c22
# Get SNP Stats
module load gcc/5.2.0 
module load qctool2/rc4-6.8
qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -incl-samples /home/jm2294/projects/RNAseq/test_run/batch1_ids.txt\
 -snp-stats\
 -osnp c22_test_bgen_snp_stats.txt
qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -incl-samples /home/jm2294/projects/RNAseq/test_run/batch1_ids.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile.bgen

qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -incl-samples /home/jm2294/projects/RNAseq/test_run/batch1_ids.txt\
 -incl-range 22:23500000-24500000\
 -excl-snpids /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/c22_filter_snps.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/impute_22_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids.bgen
