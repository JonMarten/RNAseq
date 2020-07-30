module purge                  
module load rhel7/default-peta4
module load ceuadmin/plink/2.0_09_09_18

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes

plink2\
 --bgen processing/b38_bgen/filtered/impute_6_interval_b38_filtered_no0_rnaSeqPhase1_2.bgen\
 --make-bed\
 --exclude processing/snp_inclusion_filters/plink_multiallelic_filter_SNPs.txt\
 --out INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_chr6_nofilters