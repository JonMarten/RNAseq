#!/bin/bash

#----------------------------------------------------------------------------------#
# Prepare files for PCA #
#----------------------------------------------------------------------------------# 

module load plink-bgi

# plink files for 5000: autosomes only
INFILE=/scratch/bp406/data_sets/interval_subset_olink/genotype_files/plink_format/related_5000_original_base/interval_qced_2.2.16.olink_related_subset.autosomal_only

OUTDIR=/scratch/jp549/olink5000.plink/

# text file with positions of complex regions such as HLA
COMPLEX=/home/jp549/post-doc/genetics/complex_regions.txt 

# Extract SNPs with MAF ≥1%

# plink --bfile $INFILE --missing

plink --bfile $INFILE \
      --allow-no-sex \
      --maf 0.01 \
      --geno 0.03 \
      --hwe 0.00001 midp include-nonctrl \
      --make-bed --out ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5
  
# Exclude complex regions (these are regions chr and position so this file can be used with any b37 array data)

plink --bfile ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5 \
      --allow-no-sex \
      --exclude range $COMPLEX \
      --make-bed --out ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl

#LD prune

plink --bfile ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl \
      --allow-no-sex \
      --indep 50 5 1.5 \
      --out ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning
      
#Extract only the prune in SNPs

plink --bfile ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl \
      --allow-no-sex \
      --extract ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.prune.in \
      --make-bed --out ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning
      
#Check summary stats on final bed files used for PCA
plink --bfile ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning \
      --allow-no-sex \
      --freq \
      --missing \
      --hardy midp \
      --out ${OUTDIR}olink5000.pca.file.summary
      
# 220216 SNPs used

#  summary(snp.m$F_MISS)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000800 0.001600 0.002436 0.003000 0.027400 

#  summary(hwe$P)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000101 0.2169000 0.4690000 0.4744000 0.7428000 0.9887000 

#  summary(frq$MAF)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01000 0.01932 0.03195 0.06220 0.06435 0.50000 

#----------------------------------------------------------------------------------#
# Perform MDS clustering for merged files #
#----------------------------------------------------------------------------------#

#Run genome

plink --bfile ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning \
      --allow-no-sex \
      --genome \
      --out  ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome
      
#Run PCA with first 20PC

# cluster and reaf-genome may be unnecessary

plink --bfile ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning \
      --allow-no-sex \
      --read-genome ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.genome \
      --cluster \
      --pca 20 header \
      --out ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca

#Count SNPs in each file for the log
wc -l ${OUTDIR}olink5000.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning*
      