#!/bin/bash
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 20G
#SBATCH --job-name=plink_convert_b38
#SBATCH --time=12:0:0
#SBATCH --output=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/plink_bfile/logs/plink_convert_b38_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk


#! This line enables the module command:
. /etc/profile.d/modules.sh     

module purge                  
module load rhel7/default-peta4
module load plink
module load qctool

# get start time
start=$(date +%s.%N)

cd /rds/user/jm2294/hpc-work/projects/RNAseq/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover

qctool -g /rds/user/jm2294/hpc-work/projects/RNAseq/GENETIC_DATA/b37_b38_liftover/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.bgen\
 -og /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/plink_bfile/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.vcf

plink2\
 --vcf /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/plink_bfile/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered.vcf\
 --make-bed\
 --out /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/plink_bfile/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_filtered

# get end time
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"