# Python script to run trans-eQTL mapping. Chromosome is obtained from command line arguments
import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

path = "/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/"
outpath = path + "side_projects/gmpr/gmpr_"
chr = "6"

phenotype_bed_file = path + "phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr" + chr + ".bed.gz"
gw_phenotype_bed_file = path + "phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.bed.gz"
covariates_file = path + "covariates/INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis.txt"
plink_prefix_path = path + "genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_chr6_nofilters"

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
gw_phenotype_df, gw_phenotype_pos_df = tensorqtl.read_phenotype_bed(gw_phenotype_bed_file)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates

# Read in genotypes 
pr = genotypeio.PlinkReader(plink_prefix_path)

# load genotypes and variants into data frames
genotype_df = pd.DataFrame(pr.load_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Subset to Dirk's GMPR SNPs of interest 
gmpr = pd.read_csv(outpath + "GMPR_variants_for_lookup_20200717_LS_sorted.txt", sep = '\t')
gmpr_genotype_df = genotype_df.loc[gmpr['rsID']]
gmpr_genotype_df = gmpr_genotype_df[gmpr_genotype_df['110000315494'].notnull()]
gmpr_variant_df = variant_df.loc[gmpr['rsID']]
gmpr_variant_df = gmpr_variant_df[gmpr_variant_df['chrom'].notnull()]

# Call cis-eQTLs
cis_df = cis.map_cis(gmpr_genotype_df, gmpr_variant_df, phenotype_df, phenotype_pos_df, covariates_df)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
cis_df.to_csv(outpath + "tensorqtl_cis_cisPerGene_chr" + chr + ".csv", index=True, index_label = "Phenotype")

# Cis nominal mapping
cisnom_df = cis.map_nominal(gmpr_genotype_df, gmpr_variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix=outpath + "tensorqtl_cis_cisNominal_chr" + chr)
cisnom_df2 = pd.read_parquet(outpath + "tensorqtl_cis_cisNominal_chr6.cis_qtl_pairs.6.parquet")
cisnom_df2.to_csv(outpath + "tensorqtl_cis_cisNominal_chr6.cis_qtl_pairs.6.csv", index = False)

# Call trans-eQTLs
trans_min_df = trans.map_trans(gmpr_genotype_df, gw_phenotype_df, covariates_df, return_sparse=True)
trans_min_df.to_csv(outpath + "tensorqtl_trans.csv", index = False)

# Conditional cis-analysis (may time out!)
#indep_df = cis.map_independent(gmpr_genotype_df, gmpr_variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df, nperm=10000)
#indep_df.to_csv(outpath + "tensorqtl_cis_cisIndependent_chr" + chr + ".csv", index=True, index_label = "Phenotype")