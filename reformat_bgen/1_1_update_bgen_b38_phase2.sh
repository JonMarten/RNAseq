#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --mem 10G
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/logs/genotype_processing/liftover_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH --job-name=liftover_bgen
#SBATCH -A PAUL-SL3-CPU

# Adapted from 05_make_small_bgen.sh. SNP filter lists are re-used from 02_generate_snp_ids.R

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load qctool

CHR=$SLURM_ARRAY_TASK_ID
GENPATH=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/genotypes/processing
INGEN=/home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/impute_${CHR}_interval.bgen
MIDGEN=${GENPATH}/b38_bgen/impute_${CHR}_interval_b38_no0_rnaSeqPhase1_2.bgen
OUTGEN=${GENPATH}/b38_bgen/filtered/impute_${CHR}_interval_b38_filtered_no0_rnaSeqPhase1_2.bgen
INSAMPLE=/home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples
MIDSAMPLE=${GENPATH}/b38_bgen/chr${CHR}_rnaSeqPhase1_2_sample.txt
MAP=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr${CHR}_b37_to_b38_map.txt
SAMPLEFILTER=${GENPATH}/rna_seq_phase1-2_affy_ids.txt
SNPFILTER=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/snp_inclusion_filters/c${CHR}_b38_filter_snps.txt

# Map to b38 and filter to retain only RNA seq samples
#qctool\
# -g $INGEN\
# -s $INSAMPLE\
# -incl-samples $SAMPLEFILTER\
# -map-id-data $MAP\
# -og $MIDGEN
#
## Create new sample file
#qctool\
# -g $MIDGEN\
# -os $MIDSAMPLE

# Filter SNPs to include only non-duplicated SNPs mapped to b38
qctool\
 -g $MIDGEN\
 -incl-snpids $SNPFILTER\
 -s $MIDSAMPLE\
 -og $OUTGEN