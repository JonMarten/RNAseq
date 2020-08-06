# Python script to run cis-eQTL mapping. Chromosome is obtained from command line arguments

import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

path = "/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/"
outpath = path + "results/cis_eQTLs/"
chr = "21"

phenotype_bed_file = path + "phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr" + chr + ".bed.gz"
covariates_file = path + "covariates/INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis.txt"
plink_prefix_path = path + "genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr" + chr

all_phenotype_df, all_phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)

# load genotypes and variants into data frames
genotype_df = pd.DataFrame(pr.load_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T  # samples x covariates
covariates_df = covariates_df[["age_RNA"]]
# Read in genotypes 
pr = genotypeio.PlinkReader(plink_prefix_path)


# Limit to 3 phenotypes to test conditional analysis
phenotype_df = all_phenotype_df[0:15]
phenotype_pos_df = all_phenotype_pos_df[0:15]

# Cis gene-level mapping
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0) # lambda of 0 is equivalent to BH correction
cis_df.to_csv(outpath + "tensorqtl_cis_MAF0.005_cisPerGene_3phenotest_chr" + chr + ".csv", index=True, index_label = "Phenotype")

# Cis nominal mapping
#cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, prefix=outpath + "tensorqtl_cis_MAF0.005_cisNominal_3phenotest_chr" + chr)
#cisnom_df = pd.read_parquet(outpath + "tensorqtl_cis_MAF0.005_cisNominal_3phenotest_chr21.cis_qtl_pairs.21.parquet")

# Conditional analysis
indep_df = cis.map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df, nperm=10000)
indep_df.to_csv(outpath + "tensorqtl_cis_MAF0.005_cisIndependent_3phenotest_chr" + chr + ".csv", index=True, index_label = "Phenotype")




# ENSG00000235609 on chr21 seems to break it
phenotype_df = all_phenotype_df.loc[['ENSG00000235609','ENSG00000180530','ENSG00000155313']]
phenotype_pos_df = all_phenotype_pos_df.loc[['ENSG00000235609','ENSG00000180530','ENSG00000155313']]
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
cis_df['qval'] = cis_df['pval_beta']
indep_df = cis.map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df, nperm=10000)

phenotype_df = all_phenotype_df.loc[['ENSG00000235609','ENSG00000180530']]
phenotype_pos_df = all_phenotype_pos_df.loc[['ENSG00000235609','ENSG00000180530']]
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
cis_df['qval'] = cis_df['pval_beta']
indep_df = cis.map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df, nperm=10000)


phenotype_pos_df = phenotype_pos_df.loc['ENSG00000235609']

