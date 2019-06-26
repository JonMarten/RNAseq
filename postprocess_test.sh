#!/bin/bash
#SBATCH --job-name=eqtl_chunktest_newpipe
#SBATCH --time=1:0:0
#SBATCH --cpus-per-task=1
#SBATCH --partition=short 
#SBATCH --output=/home/jm2294/projects/RNAseq/test_run_chunks/output_new_pipeline/test/eqtl_chunktest_reformat_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

source activate limix_qtl

start=$(date +%s.%N)


python /home/jm2294/projects/RNAseq/hipsci_pipeline_19_06_18/post-processing_QTL/minimal_postprocess.py\
 -id /home/jm2294/projects/RNAseq/test_run_chunks/output_new_pipeline/test\
 -od /home/jm2294/projects/RNAseq/test_run_chunks/output_new_pipeline/test 

conda deactivate 

end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"
