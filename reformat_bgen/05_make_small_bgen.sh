#SBATCH -p skylake-himem
#SBATCH --mem 10G
#SBATCH --time=12:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/filter_bgen_tinytest_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk
#!/bin/bash
#SBATCH --job-name=make_small_bgen
#SBATCH -A PETERS-SL3-CPU

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load qctool

#qctool\
# -g /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen\
# -s /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples\
# -incl-samples /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/rnaseq_affy_ids.txt\
# -map-id-data /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr${SLURM_ARRAY_TASK_ID}_b37_to_b38_map.txt\
# -og /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_no0_rnaSeqPhase1.bgen

qctool -g impute_${SLURM_ARRAY_TASK_ID}_interval_b38_no0_rnaSeqPhase1.bgen -os chr${SLURM_ARRAY_TASK_ID}_rnaSeqPhase1_sample.txt

qctool\
-g /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_no0_rnaSeqPhase1.bgen\
-incl-snpids /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/snp_inclusion_filters/c${SLURM_ARRAY_TASK_ID}_b38_filter_snps.txt\
-s /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/chr${SLURM_ARRAY_TASK_ID}_rnaSeqPhase1_sample.txt\
-og /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered_no0_rnaSeqPhase1.bgen

