# Python script to run trans-eQTL mapping
import sys
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

# Specify root directory for analysis
path = "/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/"

# Specify output path (and prefix if desired). Can be set to any folder if you want it outside the analysis directory
outpath = path + "results/trans_eQTLs/"

# Get analysis chromosome from command line arguments, this in turn comes from SLURM array ID
chr = str(sys.argv[1]) 

# Set up file paths
phenotype_bed_file = path + "phenotypes/INTERVAL_RNAseq_phase1-2_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis_chr" + chr + ".bed.gz"
covariates_file = path + "covariates/INTERVAL_RNAseq_phase1-2_fullcovariates_foranalysis.txt"
plink_prefix_path = path + "genotypes/INTERVAL_RNAseq_Phase1-2_imputed_b38_biallelic_MAF0.005_chr" + chr

# Read in phenotypes
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)

# Read in covariates and make subset to only ids that are in the phenotype file
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
covariates_df = covariates_df[phenotype_df.columns].T

# Read in genotypes 
pr = genotypeio.PlinkReader(plink_prefix_path)

# load genotypes and variants into data frames
genotype_df = pd.DataFrame(pr.load_genotypes(), index=pr.bim['snp'], columns=pr.fam['iid'])
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Call trans-eQTLs
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, return_sparse=True, maf_threshold = 0.005)
trans_df.to_csv(outpath + "tensorqtl_trans_MAF0.005_chr" + chr + ".csv", index = False)
