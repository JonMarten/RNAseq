for i in {2..21}
do
  plink --bfile /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/plink/impute_${i}_interval_b38_filtered_no0_rnaSeqPhase1\
   --exclude /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/multiallelic_filter_SNPs.txt\
   --make-bed\
   --out /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_c${i}
done