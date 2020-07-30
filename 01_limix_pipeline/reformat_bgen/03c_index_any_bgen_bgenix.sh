#!/bin/bash
#SBATCH --job-name=index_any_bgen_bgenix
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 12G
#SBATCH --time=12:0:0
#SBATCH --output=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/logs/index_any_bgen_bgenix_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge
module load rhel7/default-peta4
module load bgen

# NOTE: This script must be run with the path to the bgen file as an argument
if (( $# != 1 ))
then
  echo "Please supply the path to the bgen to be indexed as an argment, eg [sbatch script.sh /path/to/file.bgen]"
  exit 1
fi

BGEN=$1

# index bgen files
bgenix -g $BGEN -index -clobber