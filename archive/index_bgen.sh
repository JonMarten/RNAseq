#!/bin/bash
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake
#SBATCH --mem 10G
#SBATCH --job-name=index_bgen
#SBATCH --time=5:0:0
#SBATCH --output=/home/jm2294/rds/hpc-work/projects/RNAseq/GENETIC_DATA/b37_b38_liftover500kb_window/index_bgen_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4

module load bgen

# get start time
start=$(date +%s.%N)

# Create bgi index file for bgen files
cd /home/jm2294/rds/hpc-work/projects/RNAseq/GENETIC_DATA/b37_b38_liftover
CHR=$SLURM_ARRAY_TASK_ID
bgenix -g impute_${CHR}_interval_b38_filtered.bgen -index 

end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"