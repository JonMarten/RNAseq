#!/bin/bash 
#SBATCH --ntasks=1 
#SBATCH --job-name=limixtest
#SBATCH --time=1:0:0 #or –t, max wall time, in min
#SBATCH --cpus-per-task=8 #or –c, max CPUs for the job
#SBATCH --partition=long #or –p, partition to run job on 
#SBATCH --output=limixtest.err # or –o, file name to log output to 
#SBATCH --error=limixtest.out # or –e, file name to log error to

source activate limix_qtl
GENPATH=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq
PHEPATH=/home/jm2294/projects/RNAseq/test_run
OUTPATH=/home/jm2294/projects/RNAseq/test_run/output

python /home/jm2294/projects/RNAseq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen ${GENPATH}/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids\
 -af ${PHEPATH}/Feature_Annotation_Ensembl_gene_ids_autosomes.txt\
 -pf ${PHEPATH}/phenotype_5281-fc-genecounts.txt\
 -od ${OUTPATH}\
 --sample_mapping_file ${PHEPATH}/sample_mapping_file_gt_to_phe.txt\
 -c\
 -np 10000\
 -maf 0.000000001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w 250000\
 --block_size 50