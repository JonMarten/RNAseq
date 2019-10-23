#!/bin/bash
#SBATCH -A PAUL-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --job-name=eqtl_merge_results
#SBATCH --time=1:0:0
#SBATCH --mem 5G
#SBATCH --output=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/logs/eqtl_merge_results_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

# NOTE: This script must be run with the path to the limix results as an argument
if (( $# != 1 ))
then
  echo "Please supply the path to the results to be processed as an argment, eg [sbatch script.sh /path/to/folder]"
  exit 1
fi

source activate limix_qtl

# Set directory of results to merge
RESULT_DIR=$1
OUTPUT_DIR=$(echo ${RESULT_DIR}/processed)
mkdir -p $OUTPUT_DIR

# Start runtime
start=$(date +%s.%N)

# Merge limix pipeline output files together into plain text format

python /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py\
 -id $RESULT_DIR\
 -od ${OUTPUT_DIR}/processed_

conda deactivate 

# Print runtime
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")
echo "Runtime was $runtime"
