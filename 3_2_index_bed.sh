#!/bin/bash
###########################################
## NOTE: The bed file must be compressed and indexed with the commands below:
############################################

module load ceuadmin/tabix/0.2.6
cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/phenotypes

bgzip INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed && tabix -p bed INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz

# Compress and index per-chromosome files for conditional analysis
for i in {1..22}
do
  bgzip INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr${i}.bed && tabix -p bed INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr${i}.bed.gz
done