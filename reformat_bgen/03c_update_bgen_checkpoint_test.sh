#!/bin/bash
#SBATCH -J checkpoint_update_bgen
#SBATCH -A PETERS-SL3-CPU
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/checkpoint_update_bgen_%A_%a.log
#SBATCH --time=12:00:00
#SBATCH -p skylake-himem
#SBATCH --mem 20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load dmtcp/2.6.0-intel-17.0.4
module load qctool
ulimit -s 8192

RESTARTSCRIPT="dmtcp_restart_script.sh"
export DMTCP_QUIET=2

jobDir=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/checkpointing_files/checkpoint_chr$SLURM_ARRAY_TASK_ID
mkdir $jobDir
cd $jobDir

runcmd='qctool -g /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen -s /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples -map-id-data /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr${SLURM_ARRAY_TASK_ID}_b37_to_b38_map.txt -og /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_checkpoint.bgen'
qctool -g /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/impute_${SLURM_ARRAY_TASK_ID}_interval.bgen -s /home/jm2294/rds/rds-jmmh2-pre_qc_data/interval/affy_ukbiobank_array/raw_data/genetics/imputed/interval.samples -map-id-data /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/INTERVAL_chr${SLURM_ARRAY_TASK_ID}_b37_to_b38_map.txt -og /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/impute_${SLURM_ARRAY_TASK_ID}_interval_b38_checkpoint.bgen
tint=39000

echo "Start coordinator"
date
eval "dmtcp_coordinator --daemon --coord-logfile dmtcp_log.txt --exit-after-ckpt --exit-on-last -i "$tin$ sleep 2
sleep 2
cport=$(<cport.txt)
echo "$cport"

if [ -f "$RESTARTSCRIPT" ]
then
    echo "Resume the application"
    CMD="dmtcp_restart -p "$cport" -i "$tint" ckpt*.dmtcp"
    echo $CMD
    eval $CMD
else
    echo "Start the application"
    CMD="dmtcp_launch --rm --infiniband --no-gzip -h localhost -p "$cport" "$runcmd
    echo $CMD
    eval $CMD
fi

echo "Stopped program execution"
date