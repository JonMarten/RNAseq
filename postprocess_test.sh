#!/bin/bash
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --job-name=eqtl_merge_results
#SBATCH --time=1:0:0
#SBATCH --mem 5G
#SBATCH --output=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/logs/eqtl_merge_results_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

source activate limix_qtl

# Set directory of results to merge
RESULT_DIR=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/results/eqtl_test_13633258
OUTPUT_DIR=$(echo ${RESULT_DIR}/processed)
mkdir -p $OUTPUT_DIR

# Start runtime
start=$(date +%s.%N)

# Merge limix pipeline output files together into plain text format

python /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py\
 -id $RESULT_DIR\
 -od $OUTPUT_DIR

conda deactivate 

# Print runtime
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")
echo "Runtime was $runtime"
