#. /etc/profile.d/modules.sh
#module purge
#module load rhel7/default-gpu
#module load miniconda3-4.5.4-gcc-5.4.0-hivczbz 
#module load cuda/9.2
#module load r-3.6.0-gcc-5.4.0-bzuuksv
#
#source activate tensorQTL
#
#python 

import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

dir = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/"
covdir = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/"

phenotype_bed_file =  covdir + "INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_COVID19.bed.gz"
covariates_peer = covdir + "INTERVAL_RNAseq_COVID19_covariates.txt"
outdir = covdir

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_peer_df = pd.read_csv(covariates_peer, sep='\t', index_col=0).T  # samples x covariates

interaction_s = pd.read_csv("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/INTERVAL_RNAseq_COVID19_neutPCT_GxE.txt", sep = "\t", index_col=0, squeeze = True).T
interaction_s = interaction_s.squeeze()

plink_prefix_path = "/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.005"
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pd.DataFrame(pr.get_all_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# cis
# Cis gene-level mapping
pheno_df_noACE2 = phenotype_df.drop("ENSG00000130234")
phenopos_df_noACE2 = phenotype_pos_df.drop("ENSG00000130234")

pheno_df_noACE2 = pheno_df_noACE2.drop("ENSG00000184012")
phenopos_df_noACE2 = phenopos_df_noACE2.drop("ENSG00000184012")

cis_df = cis.map_cis(genotype_df, variant_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0)
cis_df.to_csv(outdir + "tensorqtl_cis_MAF0.005_cisPerGene_chr1.csv", index=True, index_label = "Phenotype")

# Cis nominal mapping
cisnom_df = cis.map_nominal(genotype_df, variant_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df, prefix= covdir + "tensorqtl_cis_MAF0.005_cisNominal_covid19")

# Conditional analysis
indep_df = cis.map_independent(genotype_df, variant_df, cis_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df)
indep_df.to_csv(covdir + "tensorqtl_cis_MAF0.005_cisIndependent_covid19.csv", index=True, index_label = "Phenotype")

# GxE
cisGxE_df = cis.map_nominal(genotype_df, variant_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df, prefix= covdir + "Test_gxe", interaction_s=interaction_s)
cis.map_nominal(genotype_df, variant_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df, prefix="tensorqtl_cis_MAF0.005_cisGxE_covid19",interaction_s=interaction_s, maf_threshold_interaction=0.005,group_s=None, run_eigenmt=True, output_dir=covdir)

for i in [8,9,21]:
  df = pd.read_parquet(covdir + "tensorqtl_cis_MAF0.005_cisGxE_covid19.cis_qtl_pairs." + str(i) + ".parquet")
  df.to_csv(covdir + "tensorqtl_cis_MAF0.005_cisGxE_covid19.cis_qtl_pairs." + str(i) + ".csv", index=False)

# trans
trans_peer_df = trans.map_trans(genotype_df, pheno_df_noACE2, covariates_peer_df, return_sparse=True, maf_threshold = 0.005)
trans_peer_df.to_csv(outdir + "tensorqtl_trans_MAF0.005_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19.csv")

#################################################################
# chrX
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

# Read genotypes
plink_prefix_path = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes/INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2_deduplicated_MAF0.005"
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pd.DataFrame(pr.get_all_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

ace2_plink_prefix_path = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/covid19/genotypes/INTERVAL_chrX_merged_cleaned_RNAseq_phase1-2_deduplicated_MAF0.005_ACE2nonzero"
ace2_pr = genotypeio.PlinkReader(ace2_plink_prefix_path)
ace2_genotype_df = pd.DataFrame(ace2_pr.get_all_genotypes(), index=ace2_pr.bim['snp'], columns=ace2_pr.fam['iid'])
ace2_variant_df = ace2_pr.bim.set_index('snp')[['chrom', 'pos']]

# Cis gene-level mapping
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
cis_df_out = add_rsid(cis_df, rsid)
cis_df_out.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cis.csv", index=True, index_label = "Phenotype")

ace2_cis_df = cis.map_cis(ace2_genotype_df, ace2_variant_df, ace2_phenotype_df, ace2_phenotype_pos_df, ace2_covariates_df)
ace2_cis_df_out = add_rsid(ace2_cis_df, rsid)
ace2_cis_df_out.to_csv(covdir + "/results/INTERVAL_ACE2_nonzeros_MAF0.005_cis.csv", index=True, index_label = "Phenotype")

# Cis nominal mapping
cisnom_df = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix = covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal")
cisnom_df = pd.read_parquet(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal.cis_qtl_pairs.23.parquet", )
cisnom_df = add_rsid(cisnom_df, rsid)
cisnom_df.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal.cis_qtl_pairs.23.csv", index=False)

ace2_cisnom_df = cis.map_nominal(ace2_genotype_df, ace2_variant_df, ace2_phenotype_df, ace2_phenotype_pos_df, ace2_covariates_df, prefix = covdir + "/results/INTERVAL_ACE2_nonzeros_MAF0.005_cisNominal")
ace2_cisnom_df = pd.read_parquet(covdir + "/results/INTERVAL_ACE2_nonzeros_MAF0.005_cisNominal.cis_qtl_pairs.23.parquet")
ace2_cisnom_df = add_rsid(ace2_cisnom_df, rsid)
ace2_cisnom_df.to_csv(covdir + "/results/INTERVAL_ACE2_nonzeros_MAF0.005_cisNominal.cis_qtl_pairs.23.csv", index=False)


# Conditional analysis (can't run on ace2 nonzeros as no significant P-values
# Add qval column manually to avoid breaking funtion (not intended for use on a single phenotype)
cis_df['qval'] = cis_df['pval_perm']
indep_df = cis.map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df)
indep_df = add_rsid(indep_df, rsid)
indep_df.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisConditional.csv", index=True, index_label = "Phenotype")

# trans
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, return_sparse=True, maf_threshold = 0.03)
trans_peer_df.to_csv(outdir + "tensorqtl_trans_MAF0.03_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_CHRX.csv")

# split by sex
male = covariates_df["sex"] == 1
maleids = pd.Series(covariates_df[male].index)
male_cis_df = cis.map_cis(genotype_df[maleids], variant_df, phenotype_df[maleids], phenotype_pos_df, covariates_df[male])
male_cis_df = add_rsid(male_cis_df, rsid)
male_cis_df.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cis_male.csv", index=True, index_label = "Phenotype")

male_cisnom_df = cis.map_nominal(genotype_df[maleids], variant_df, phenotype_df[maleids], phenotype_pos_df, covariates_df[male], prefix = covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal_male")
male_cisnom_df = pd.read_parquet(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal_male.cis_qtl_pairs.23.parquet", )
male_cisnom_df = add_rsid(male_cisnom_df, rsid)
male_cisnom_df.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal_male.cis_qtl_pairs.23.csv", index=False)


female = covariates_df["sex"] == 0
femaleids = pd.Series(covariates_df[female].index)
female_cis_df = cis.map_cis(genotype_df[femaleids], variant_df, phenotype_df[femaleids], phenotype_pos_df, covariates_df[female])
female_cis_df = add_rsid(female_cis_df, rsid)
female_cis_df.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cis_female.csv", index=True, index_label = "Phenotype")

female_cisnom_df = cis.map_nominal(genotype_df[femaleids], variant_df, phenotype_df[femaleids], phenotype_pos_df, covariates_df[female], prefix = covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal_female")
female_cisnom_df = pd.read_parquet(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal_female.cis_qtl_pairs.23.parquet", )
female_cisnom_df = add_rsid(female_cisnom_df, rsid)
female_cisnom_df.to_csv(covdir + "/results/INTERVAL_ACE2_MAF0.005_cisNominal_female.cis_qtl_pairs.23.csv", index=False)



# genome-wide trans


#gw_plink_prefix_path = "/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_AllAutosomes"
#gw_pr = genotypeio.PlinkReader(gw_plink_prefix_path)
#gw_genotype_df = pd.DataFrame(gw_pr.get_all_genotypes(), index=gw_pr.bim['snp'], columns=gw_pr.fam['iid'])
#gw_variant_df = gw_pr.bim.set_index('snp')[['chrom', 'pos']]
#
#gw_trans_df = trans.map_trans(gw_genotype_df, phenotype_df, covariates_df, return_sparse=True, maf_threshold = 0.03)


for i in range(1,23):
  print(i)
  plink_prefix_path = "/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr" + str(i)
  print(plink_prefix_path)
  gw_pr = genotypeio.PlinkReader(plink_prefix_path)
  gw_genotype_df = pd.DataFrame(gw_pr.get_all_genotypes(), index=gw_pr.bim['snp'], columns=gw_pr.fam['iid'])
  gw_variant_df = gw_pr.bim.set_index('snp')[['chrom', 'pos']]
  
  gw_trans_df = trans.map_trans(gw_genotype_df, phenotype_df, covariates_df, return_sparse=True, maf_threshold = 0.03)
  gw_trans_min_df.to_csv(outdir + "tensorqtl_trans_MAF0.03_all_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_CHR" + str(i) + ".csv")
  
  

gw_genotype_df.to_csv("genotype_df.csv")
gw_variant_df.to_csv("variant_df.csv")



# GxE
cisGxE_df = cis.map_nominal(genotype_df, variant_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df, prefix= covdir + "Test_gxe", interaction_s=interaction_s)
cis.map_nominal(genotype_df, variant_df, pheno_df_noACE2, phenopos_df_noACE2, covariates_peer_df, prefix="tensorqtl_cis_MAF0.005_cisGxE_covid19",interaction_s=interaction_s, maf_threshold_interaction=0.005,group_s=None, run_eigenmt=True, output_dir=covdir)
trans_peer_df = trans.map_trans(genotype_df, phenotype_df, covariates_peer_df, return_sparse=True, maf_threshold = 0.005, pval_threshold = 0.05)
trans_peer_df.to_csv(outdir + "results/tensorqtl_trans_MAF0.005_chrx_age_sex_rin_batch_readDepth_PC10_PEER20_COVID19_phase1.csv")
p = phenotype_df.drop(phenotype_df.index[[0,1,2]])
p_pos = phenotype_pos_df.drop(phenotype_pos_df.index[[0,1,2]])
cis_df = cis.map_cis(genotype_df, variant_df, p, p_pos, covariates_peer_df)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0)
cis_df.to_csv(outdir + "results/tensorqtl_cis_MAF0.005_cisPerGene_chr1.csv", index=True, index_label = "Phenotype")
cisnom_df = cis.map_nominal(genotype_df, variant_df, p, p_pos, covariates_peer_df, prefix= covdir + "results/tensorqtl_cis_MAF0.005_cisNominal_ACE2_phase1")
cisnom_df = pd.read_parquet(covdir + "results/tensorqtl_cis_MAF0.005_cisNominal_ACE2_phase1.cis_qtl_pairs.23.parquet")
cisnom_df.to_csv(covdir + "results/tensorqtl_cis_MAF0.005_cisNominal_ACE2_phase1.cis_qtl_pairs.23.csv", index=False)
