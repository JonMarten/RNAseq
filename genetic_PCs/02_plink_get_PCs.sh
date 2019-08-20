#!/bin/bash
#SBATCH --job-name=plink_calc_PCs
#SBATCH -p skylake-himem
#SBATCH -A PETERS-SL3-CPU
#SBATCH --time=12:0:0
#SBATCH --mem=15G
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/genetic_PCs/plink_calc_PCs_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load plink

plink2\
 --bfile /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/merged/merged_samplecleaned\
 --keep /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/genetic_PCs/rna_seq_2.7k_affy_ids.txt\
 --pca 20 'var-wts'\
 --make-rel\
 --out /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/genetic_PCs/INTERVAL_genotype_PCs_phase1