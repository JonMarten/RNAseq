#!/bin/bash
#SBATCH -A PETERS-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 120G
#SBATCH --job-name=eqtl_test
#SBATCH --time=12:0:0
#SBATCH --output=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/logs/eqtl_test_fullrun_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

# get start time
#start=$(date +%s)

# Get genomic positions for chunk
CHUNK=$(head /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/chunklist_b38_50genes.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1)
CHR=$(echo $CHUNK | cut -d ' ' -f1)
START=$(echo $CHUNK | cut -d ' ' -f2)
END=$(echo $CHUNK | cut -d ' ' -f3)

# Load limix environment
source activate limix_qtl

# Specify file paths. Current config creates a new output directory for every new job submitted
GENPATH=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered
PHEPATH=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing
OUTPATH=/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/analysis/00_testing/results/eqtl_test_fullrun_${SLURM_ARRAY_JOB_ID}

# Create output directory if it doesn't exist
mkdir -p $OUTPATH

# Specify files. NOTE THAT GENFILE DOES NOT NEED .bgen SUFFIX
GENFILE=${GENPATH}/impute_${CHR}_interval_b38_filtered
ANFILE=${PHEPATH}/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt
PHEFILE=${PHEPATH}/phenotype/phenotype_5281-fc-genecounts.txt
SAMPLEMAPFILE=${PHEPATH}/phenotype/sample_mapping_file_gt_to_phe.txt
COVFILE=${PHEPATH}/covariates/INTERVAL_RNA_batch1_2_covariates_sex_age.txt
GR=$(echo ${CHR}:${START}-${END})
BLOCKSIZE=3000
WINDOW=500000
PERMUTATIONS=100
MAF=0.01

# Echo config for log file
echo Running Limix
echo "************** Parameters **************"
echo Chunk: $SLURM_ARRAY_TASK_ID
echo Genomic Region: chr$GR
echo Genotype File: ${GENFILE}.bgen
echo Phenotype File: $PHEFILE
echo Annotation File: $ANFILE
echo Sample Map File: $SAMPLEMAPFILE
echo Covariate File: $COVFILE
echo Block Size: $BLOCKSIZE
echo Window: $WINDOW
echo Permutations: $PERMUTATIONS
echo MAF: $MAF
echo "****************************************"

# Run QTL mapping
python -u /home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen $GENFILE\
 -af $ANFILE\
 -pf $PHEFILE\
 -od $OUTPATH\
 --sample_mapping_file $SAMPLEMAPFILE\
 -c\
 -np $PERMUTATIONS\
 -maf $MAF\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w $WINDOW\
 --block_size $BLOCKSIZE\
 -gr $GR\
 -cf $COVFILE

conda deactivate 
 
#end=$(date +%s)    
#runtime=$(python -c "print(round((${end} - ${start})/3600,3))")

#echo "Runtime was $runtime hours"
