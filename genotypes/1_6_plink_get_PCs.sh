#!/bin/bash
#SBATCH --job-name=plink_calc_PCs_benscript
#SBATCH -p skylake-himem
#SBATCH -A PETERS-SL3-CPU
#SBATCH --time=12:0:0
#SBATCH --mem=15G
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/genotype_processing/plink_calc_PCs_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk


#----------------------------------------------------------------------------------#
# Prepare files for PCA #
#----------------------------------------------------------------------------------# 

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load plink
module load plink-1.9-gcc-5.4.0-sm3ojoi

# plink files for 5000: autosomes only
INFILE=/home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/merged/merged_samplecleaned
OUTDIR=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes
  
# text file with positions of complex regions such as HLA
COMPLEX=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/PCs/processing/complex_regions.txt 

# Extract SNPs with MAF ≥1% and filter to retain only phase 1 RNA seq samples
plink2 \
 --bfile $INFILE\
 --keep ${OUTDIR}/processing/rna_seq_phase1-2_affy_ids.txt\
 --maf 0.01\
 --geno 0.03\
 --hwe 0.00001 midp\
 --make-bed --out ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5

# Exclude complex regions (these are regions chr and position so this file can be used with any b37 array data)

plink2 --bfile ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5 \
--exclude range $COMPLEX \
--make-bed --out ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl

#LD prune, this gernates two lists of SNPs, one that's kept (in) and one that's pruned (out)

plink --bfile ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl \
--indep 50 5 1.5 \
--out ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning

#Extract only the prune in SNPs

plink2 --bfile ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl \
--extract ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.prune.in \
--make-bed --out ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning

#Check summary stats on final bed files used for PCA
plink2 --bfile ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning \
--freq \
--missing \
--hardy midp \
--out ${OUTDIR}/PCs/processing/rnaseq_phase1.pca.file.summary

#----------------------------------------------------------------------------------#
# Perform MDS clustering for merged files #
#----------------------------------------------------------------------------------#

#Run genome

plink --bfile ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning \
--genome \
--out  ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome

#Run PCA with first 20PC

# cluster and reaf-genome may be unnecessary ## JCLM: I think it is. PCA alone is all you need.

plink2 --bfile ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning \
--pca 20 var-wts 'vcols=+pos' \
--make-rel 'square' \
--out ${OUTDIR}/PCs/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca
  
#Count SNPs in each file for the log
wc -l ${OUTDIR}/PCs/processing/rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning* >>  ${OUTDIR}/PCs/processing/file_lengths.txt
  