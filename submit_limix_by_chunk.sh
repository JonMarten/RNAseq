#!/bin/bash
#SBATCH --job-name=eqtl_chunktest
#SBATCH --time=48:0:0 #or –t, max wall time, in min
#SBATCH --cpus-per-task=1 #or –c, max CPUs for the job
#SBATCH --partition=long #or –p, partition to run job on
#SBATCH --output=eqtl_chunktest.out # or –o, file name to log output to
#SBATCH --error=eqtl_chunktest.err # or –e, file name to log error to

# get start time
start=$(date +%s.%N)

# Get genomic positions for chunk
CHUNK=$(head /home/jm2294/projects/RNAseq/test_run_chunks/chunklist.txt -n $SLURM_ARRAY_JOB_ID | tail -n 1)
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
echo Chunk: chr$GR
echo Genotype: $GENFILE
echo Phenotype: $PHEFILE
echo Annotation: $ANFILE
echo Sample Map: $SAMPLEMAPFILE
echo Covariates: $COVFILE

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
 -od $OUTHPATH 

conda deactivate 
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"