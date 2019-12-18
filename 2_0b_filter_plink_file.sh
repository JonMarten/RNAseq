#!/bin/bash
#SBATCH -A PAUL-SL3-CPU
#SBATCH -p skylake
#SBATCH --mem 150G
#SBATCH --job-name=plink_maf_filter
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/logs/plink_merge_chromosomes_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge                  
module load rhel7/default-peta4
module load plink-1.9-gcc-5.4.0-sm3ojoi

cd /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/genotypes

plink\
 --memory 191000\
 --bfile INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all\
 --make-bed\
 --maf 0.01\
 --out INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.01\
 --hwe 1e-6

