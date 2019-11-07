# Submit limix jobs for chr22
sbatch -a 1134-1167 02_submit_limix.sh age_sex_rin_batch_PC10
sbatch -a 1134-1167 02_submit_limix.sh age_sex_rin_batch_PC10_PEER20
sbatch -a 1134-1167 02_submit_limix.sh age_sex_rin_batch_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT

# Postprocess results
sbatch 03_postprocess_limix.sh /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10
sbatch 03_postprocess_limix.sh /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20
sbatch 03_postprocess_limix.sh /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT

# Merge results into single file
Rscript 04_merge_postprocessed_output.R

# Get significant eSNPs and comparison to eQTLgen
Rscript 05_eQTLgen_comparison.R



sbatch -a 1-1133 02_submit_limix.sh age_sex_rin_batch_PC10_PEER20