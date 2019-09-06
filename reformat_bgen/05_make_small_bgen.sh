#!/bin/bash
#SBATCH -p skylake-himem
#SBATCH --mem 10G
#SBATCH --time=12:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/make_small_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#SBATCH --job-name=make_small_bgen
#SBATCH -A PETERS-SL3-CPU

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load qctool

CHR=$SLURM_ARRAY_TASK_ID
GENPATH=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover
INGEN=/home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/impute_${CHR}_interval.bgen
MIDGEN=${GENPATH}/b38_bgen/impute_${CHR}_interval_b38_no0_rnaSeqPhase1.bgen
OUTGEN=${GENPATH}/b38_bgen/filtered/impute_${CHR}_interval_b38_filtered_no0_rnaSeqPhase1.bgen
INSAMPLE=/home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples
MIDSAMPLE=${GENPATH}/b38_bgen/chr${CHR}_rnaSeqPhase1_sample.txt
MAP=${GENPATH}/INTERVAL_chr${CHR}_b37_to_b38_map.txt
SAMPLEFILTER=${GENPATH}/rnaseq_affy_ids.txt
SNPFILTER=${GENPATH}/snp_inclusion_filters/c${CHR}_b38_filter_snps.txt

# Map to b38 and filter to retain only RNA seq samples
qctool\
 -g $INGEN\
 -s $INSAMPLE\
 -incl-samples $SAMPLEFILTER\
 -map-id-data $MAP\
 -og $MIDGEN

# Create new sample file
qctool\
 -g $MIDGEN\
 -os $MIDSAMPLE

qctool\
 -g $MIDGEN\
 -incl-snpids $SNPFILTER\
 -s $MIDSAMPLE\
 -og $OUTGEN