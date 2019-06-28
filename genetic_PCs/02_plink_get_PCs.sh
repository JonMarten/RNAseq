#!/bin/bash
module load plink/2.0_09_09_18

plink2\
 --bfile /scratch/curated_genetic_data/interval/genotyped/interval_qced_24.8.18\
 --keep /home/jm2294/projects/RNAseq/genetic_PCs/rna_seq_5k_affy_ids.txt\
 --pca\
 --out /home/jm2294/projects/RNAseq/genetic_PCs/INTERVAL_genotype_PCs