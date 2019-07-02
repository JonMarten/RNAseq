#!/bin/bash
#SBATCH --job-name=plink_calc_PCs
#SBATCH --time=24:0:0
#SBATCH --partition=medium 
#SBATCH --output=/home/jm2294/projects/RNAseq/genetic_PCs/plink_calc_PCs_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

module load plink/2.0_09_09_18

plink2\
 --bfile /scratch/curated_genetic_data/interval/genotyped/interval_qced_24.8.18\
 --keep /home/jm2294/projects/RNAseq/genetic_PCs/rna_seq_5k_affy_ids.txt\
 --pca 20 'var-wts'\
 --make-rel\
 --out /home/jm2294/projects/RNAseq/genetic_PCs/INTERVAL_genotype_PCs