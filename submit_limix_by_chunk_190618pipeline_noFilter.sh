#!/bin/bash
#SBATCH --job-name=eqtl_chunktest_newpipe_nofilter
#SBATCH --time=24:0:0
#SBATCH --cpus-per-task=10
#SBATCH --partition=medium 
#SBATCH --output=/home/jm2294/projects/RNAseq/test_run_chunks/output_newpipeline_nofilter/eqtl_chunktest_newpipe_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk


# get start time
start=$(date +%s.%N)

# Get genomic positions for chunk
CHUNK=$(head /home/jm2294/projects/RNAseq/test_run_chunks/chunklist_b38.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1)
#CHUNK=$(head /home/jm2294/projects/RNAseq/test_run_chunks/chunklist.txt -n 529 | tail -n 1)  ## TEST LINE FOR INTERACTIVE RUN
	
CHR=$(echo $CHUNK | cut -d ' ' -f1)
START=$(echo $CHUNK | cut -d ' ' -f2)
END=$(echo $CHUNK | cut -d ' ' -f3)

# Load limix environment
source activate limix_qtl

# Specify file paths
GENPATH=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq
PHEPATH=/home/jm2294/projects/RNAseq
OUTPATH=/home/jm2294/projects/RNAseq/test_run_chunks/output_newpipeline_nofilter

# Specify files. NOTE THAT GENFILE DOES NOT NEED .bgen SUFFIX
#GENFILE=${GENPATH}/impute_${CHR}_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids
GENFILE=${GENPATH}/b37_b38_liftover/impute_${CHR}_interval_b38_filtered_on_rsids
ANFILE=${PHEPATH}/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt
PHEFILE=${PHEPATH}/test_run/phenotype_5281-fc-genecounts.txt
SAMPLEMAPFILE=${PHEPATH}/test_run/sample_mapping_file_gt_to_phe.txt
COVFILE=${PHEPATH}/test_run/INTERVAL_RNA_batch1_2_covariates_sex_age.txt
GR=$(echo ${CHR}:${START}-${END})
BLOCKSIZE=3000
WINDOW=1000000
PERMUTATIONS=100

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
echo Block Size: $BLOCKSIZE
echo Window: $WINDOW
echo Permutations: $PERMUTATIONS
echo "****************************************"

# Run QTL mapping
python -u /home/jm2294/projects/RNAseq/hipsci_pipeline_19_06_18/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen $GENFILE\
 -af $ANFILE\
 -pf $PHEFILE\
 -od $OUTPATH\
 --sample_mapping_file $SAMPLEMAPFILE\
 -c\
 -np $PERMUTATIONS\
 -maf 0.001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w $WINDOW\
 --block_size $BLOCKSIZE\
 -gr $GR\
 -cf $COVFILE
# --variant_filter /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/b38_biallelic_snps_only_no_indels.txt
# --variant_filter /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/b38_indels_only.txt

#python -u /home/jm2294/projects/RNAseq/hipsci_pipeline/post-processing_QTL/minimal_postprocess.py\
# -id $OUTPATH\
# -od $OUTPATH 

conda deactivate 
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"
