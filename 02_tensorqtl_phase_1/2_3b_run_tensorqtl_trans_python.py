import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

dir = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/"

phenotype_bed_file =  dir + "phenotypes/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz"

covariates_minimal = dir + "covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10.txt"
covariates_peer = dir + "covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_PEER20.txt"
covariates_cell = dir + "covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.txt"

outdir = dir + "results/python_module_method/"

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_min_df = pd.read_csv(covariates_minimal, sep='\t', index_col=0).T  # samples x covariates
covariates_peer_df = pd.read_csv(covariates_peer, sep='\t', index_col=0).T  # samples x covariates
covariates_cell_df = pd.read_csv(covariates_cell, sep='\t', index_col=0).T  # samples x covariates

interaction_s = pd.read_csv("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/covariates/INTERVAL_RNAseq_phase1_GxE_neutPCT.txt", sep = "\t", index_col=0, squeeze = True).T

for i in range(1,23):
  print(i)
  plink_prefix_path = "/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_chr" + str(i) + "_MAF0.005"
  print(plink_prefix_path)
  pr = genotypeio.PlinkReader(plink_prefix_path)
  genotype_df = pd.DataFrame(pr.get_all_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
  variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
  
  trans_min_df = trans.map_trans(genotype_df, phenotype_df, covariates_min_df, return_sparse=True, maf_threshold = 0.005)
  trans_min_df.to_csv(outdir + "tensorqtl_trans_MAF0.005_chr" + str(i) + "_age_sex_rin_batch_readDepth_PC10.csv")
  
  trans_peer_df = trans.map_trans(genotype_df, phenotype_df, covariates_peer_df, return_sparse=True, maf_threshold = 0.005)
  trans_peer_df.to_csv(outdir + "tensorqtl_trans_MAF0.005_chr" + str(i) + "_age_sex_rin_batch_readDepth_PC10_PEER20.csv")
  
  trans_cell_df = trans.map_trans(genotype_df, phenotype_df, covariates_cell_df, return_sparse=True, maf_threshold = 0.005)
  trans_cell_df.to_csv(outdir + "tensorqtl_trans_MAF0.005_chr" + str(i) + "_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT.csv")

