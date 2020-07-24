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

for i in {1..22}
do
  echo "Running chr ${i}"
  plink\
   --memory 191000\
   --bfile INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_all_MAF0.005\
   --chr $i\
   --make-bed\
   --out INTERVAL_b38_autosomes_RNAseqPhase1_biallelic_chr${i}_MAF0.005
done

