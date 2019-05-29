#!/bin/bash
#SBATCH --job-name=eqtl_chunktest
#SBATCH --time=48:0:0
#SBATCH --cpus-per-task=5 
#SBATCH --partition=long 
#SBATCH --output=eqtl_chunktest_%A_%a.log

# get start time
start=$(date +%s.%N)

# Get genomic positions for chunk
CHUNK=$(head /home/jm2294/projects/RNAseq/test_run_chunks/chunklist.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1)
#CHUNK=$(head /home/jm2294/projects/RNAseq/test_run_chunks/chunklist.txt -n 529 | tail -n 1)  ## TEST LINE FOR INTERACTIVE RUN

CHR=$(echo $CHUNK | cut -d ' ' -f1)
START=$(echo $CHUNK | cut -d ' ' -f2)
END=$(echo $CHUNK | cut -d ' ' -f3)

# Load limix environment
source activate limix_qtl

# Specify file paths
GENPATH=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq
PHEPATH=/home/jm2294/projects/RNAseq/test_run
OUTPATH=/home/jm2294/projects/RNAseq/test_run_chunks/output

# Specify files
GENFILE=${GENPATH}/impute_${CHR}_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids
ANFILE=${PHEPATH}/Feature_Annotation_Ensembl_gene_ids_autosomes.txt
PHEFILE=${PHEPATH}/phenotype_5281-fc-genecounts.txt
SAMPLEMAPFILE=${PHEPATH}/sample_mapping_file_gt_to_phe.txt
COVFILE=${PHEPATH}/INTERVAL_RNA_batch1_2_covariates_sex_age.txt
GR=$(echo ${CHR}:${START}-${END})

# Echo config for log file
echo Running Limix
echo "************** Parameters **************"
echo Chunk: $SLURM_ARRAY_TASK_ID
echo Genomic Region: chr$GR
echo Genotype File: $GENFILE
echo Phenotype File: $PHEFILE
echo Annotation File: $ANFILE
echo Sample Map File: $SAMPLEMAPFILE
echo Covariate File: $COVFILE
echo "****************************************"

# Run QTL mapping
python -u /home/jm2294/projects/RNAseq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen $GENFILE\
 -af $ANFILE\
 -pf $PHEFILE\
 -od $OUTPATH\
 --sample_mapping_file $SAMPLEMAPFILE\
 -cf $COVFILE\
 -c\
 -np 1000\
 -maf 0.001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w 500000\
 --block_size 1500\
 -gr $GR

python -u /home/jm2294/projects/RNAseq/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py\
 -id $OUTPATH\
 -od $OUTPATH 

conda deactivate 
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"
