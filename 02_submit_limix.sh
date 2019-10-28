#!/bin/bash
#SBATCH -A PAUL-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 150G
#SBATCH --job-name=eqtl_phase1_cis_18373genes_20PEER
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/logs/eqtl_phase1_cis_18373genes_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

# get start time
start=$(date +%s.%N)

# Specify number of genes per chunk. 15 or 5 at the moment. 15 times out on SL3 nodes for some chunks.
CHUNKSIZE=5

# Get genomic positions for chunk
CHUNK=$(head /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/annotation_file/chunklist_b38_${CHUNKSIZE}genes_filtered.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1)
CHR=$(echo $CHUNK | cut -d ' ' -f1)
START=$(echo $CHUNK | cut -d ' ' -f2)
END=$(echo $CHUNK | cut -d ' ' -f3)

# Load limix environment
source activate limix_qtl

# Specify input file paths. Current config creates a new output directory for every new job submitted
GENPATH=/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered
PHEPATH=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping


COVS=age_sex_rin_batch_PC10_PEER20

# Specify files. NOTE THAT GENFILE DOES NOT NEED .bgen SUFFIX
GENFILE=${GENPATH}/impute_${CHR}_interval_b38_filtered_no0_rnaSeqPhase1
ANFILE=${PHEPATH}/annotation_file/Feature_Annotation_Ensembl_gene_ids_autosomes_b38.txt
PHEFILE=${PHEPATH}/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.txt
SAMPLEMAPFILE=${PHEPATH}/phenotype/sample_mapping_file_gt_to_phe_phase1.txt
COVFILE=${PHEPATH}/covariates/INTERVAL_RNAseq_phase1_${COVS}.txt
GR=$(echo ${CHR}:${START}-${END})
BLOCKSIZE=2000	
WINDOW=500000
PERMUTATIONS=500
MAF=0.005

# Create output directory if it doesn't exist
OUTPATH=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_${COVS}_${CHUNKSIZE}GenesPerChunk
mkdir -p $OUTPATH

# Echo config for log file
echo Running Limix
echo "************** Parameters **************"
echo Genes per chunk: $CHUNKSIZE
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
 
end=$(date +%s.%N)
 
runtime=$(python -c "print(${end} - ${start})")

#echo "Runtime was $runtime hours"
