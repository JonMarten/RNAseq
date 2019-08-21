#!/bin/bash
#SBATCH --job-name=filter_bgen
#SBATCH -A PAUL-SL2-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 10G
#SBATCH --time=36:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/filter_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

# get start time
start=$(date +%s.%N)

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load qctool

# use unique identifiers to retain only non-duplicated SNPs mapped to b38
qctool\
 -g /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_no0.bgen\
 -incl-snpids /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/snp_inclusion_filters/c${SLURM_ARRAY_TASK_ID}_b38_filter_snps.txt\
 -s /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples\
 -og /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered_no0.bgen
 
end=$(date +%s.%N)

runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"