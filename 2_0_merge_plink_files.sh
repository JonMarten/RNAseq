#!/bin/bash
#SBATCH -A PAUL-SL2-CPU
#SBATCH -p skylake
#SBATCH --mem 100G
#SBATCH --job-name=plink_merge_1_22
#SBATCH --time=36:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/logs/plink_merge_chromosomes_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load plink-1.9-gcc-5.4.0-sm3ojoi

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes

# SNP filter list was generated from a failed merge run 3 SNPs with multiple positions and 35k multiallelic SNPS. Maybe retain these at a later date?

plink\
 --merge-list /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/chr1_22_filelist.txt\
 --make-bed\
 --exclude /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/multiallelic_filter_SNPs.txt\
 --out /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes/INTERVAL_b38_autosomes_RNAseqPhase1_biallelic


