###pip install git+https://github.com/broadinstitute/tensorqtl.git#egg=tensorqtl

#stty kill ^U
#stty erase ^H
#. /etc/profile.d/modules.sh
#module purge
#module load rhel7/default-gpu
#module load miniconda3-4.5.4-gcc-5.4.0-hivczbz 
#module load cuda/9.2
#module load r-3.6.0-gcc-5.4.0-bzuuksv
#source activate tensorQTL
#python 

import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

# Function to re-add RSids as these are no longer in the vcf file
def add_rsid(data_df, rsid_df):
  data_df['pheno_id'] = data_df.index
  out = pd.merge(data_df, rsid_df, on = 'variant_id', how = 'left')
  return out

rsid = pd.read_csv("~/rds/rds-jmmh2-projects/covid/ace2/interval_genetic_data/interval_imputed_data/INTERVAL_imp_rsIDs.txt", sep=" ")
rsid.columns = ['variant_id', "rsid"]

covdir = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19"

phenotype_file =  covdir + "/phenotypes/INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed.gz"
covariates_file = covdir + "/covariates/INTERVAL_RNAseq_COVID19_covariates.txt"

phenotype_file_ace2 =  covdir + "/phenotypes/INTERVAL_RNAseq_phase1-2_UNfilteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_ACE2_no_zeros.bed.gz"
covariates_file_ace2 = covdir + "/covariates/INTERVAL_RNAseq_COVID19_covariates_ACE2_no_zeros.txt"

outdir = covdir + "/results"

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_file)
ace2_phenotype_df, ace2_phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_file_ace2)

# drop all phenotypes but ACE2
phenotype_df = phenotype_df.drop(phenotype_df.index[[0,1,2]])
phenotype_pos_df = phenotype_pos_df.drop(phenotype_pos_df.index[[0,1,2]])

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates
ace2_covariates_df = pd.read_csv(covariates_file_ace2, sep='\t', index_col=0).T  # samples x covariates

rootpath = '/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/'

for i in range(1,23):
  print(i)
  plink_prefix_path = rootpath + 'genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr' + str(i)
  print(plink_prefix_path)
  gw_pr = genotypeio.PlinkReader(plink_prefix_path)
  gw_genotype_df = pd.DataFrame(gw_pr.load_genotypes(), index=gw_pr.bim['snp'], columns=gw_pr.fam['iid'])
  gw_variant_df = gw_pr.bim.set_index('snp')[['chrom', 'pos']]
  MAF_filter = 0.005
  gw_trans_df = trans.map_trans(gw_genotype_df, phenotype_df, covariates_df, return_sparse = True, return_r2=True, maf_threshold = MAF_filter, batch_size = gw_variant_df.shape[0])
  gw_trans_df.to_csv(outdir + "/tensorqtl_trans_MAF" + str(MAF_filter) + "_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_CHR" + str(i) + ".csv")
  
# chrX trans
plink_prefix_path_x = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes/INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2_deduplicated_MAF0.005"
x_pr = genotypeio.PlinkReader(plink_prefix_path_x)
x_genotype_df = pd.DataFrame(x_pr.load_genotypes(), index=x_pr.bim['snp'], columns=x_pr.fam['iid'])
x_variant_df = x_pr.bim.set_index('snp')[['chrom', 'pos']]
MAF_filter = 0.005
x_trans_df = trans.map_trans(x_genotype_df, phenotype_df, covariates_df, return_sparse = True, return_r2=True, maf_threshold = MAF_filter, batch_size = x_variant_df.shape[0])
x_trans_df.to_csv(outdir + "/tensorqtl_trans_MAF" + str(MAF_filter) + "_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_CHRX.csv")
  
  
  
  
  
gw_perm = trans.map_permutations(gw_genotype_df, covariates_df, permutations=None, chr_s=None, nperms=10000, maf_threshold=0,batch_size=20000, logger=None, seed=None, verbose=True)
  

gw_genotype_df.to_csv("genotype_df.csv")
gw_variant_df.to_csv("variant_df.csv")
  

# Permutations?
gw_perm = trans.map_permutations(gw_genotype_df, covariates_df, permutations=None, chr_s=None, nperms=10000, maf_threshold=0,batch_size=20000, logger=None, seed=None, verbose=True)
app_perm = trans.apply_permutations(gw_perm, gw_trans_df)

# Find dodgy SNPs
import numpy as np
gw_genotype_df.index.get_loc("rs202244392")
bad_snps = gw_genotype_df[1608:1609]
bad_snps = bad_snps.rename(index = {"rs202244392" : "badsnp"})
bad_cov = covariates_df[["age_RNA","sex"]]
bad_cov['c1'] = np.random.randint(0,5, size=len(bad_cov))
bad_cov['c2'] = np.random.randint(0,5, size=len(bad_cov))
bad_cov = bad_cov[["c1","c2"]]

bad_phe = phenotype_df
bad_phe = bad_phe.rename(index = {"ENSG00000130234" : "badphe"})
#bad_phe = phenotype_df.T
#bad_phe['p1'] = np.random.uniform(-11.939618,-5.938530, size=len(bad_phe))
#bad_phe = bad_phe["p1"]
#bad_phe = pd.DataFrame(bad_phe).T

bad_trans_df = trans.map_trans(bad_snps,  bad_phe, bad_cov, return_sparse=False)
smol_bad_snps = bad_snps.T[1:10]
gw_trans_df = trans.map_trans(smol_bad_snps, phenotype_df, covariates_df, return_sparse=True, maf_threshold = MAF_filter)

