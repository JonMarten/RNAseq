source activate limix_qtl

snakemake 
	--snakefile ./snakemake 
	--jobs 50 
	--latency-wait 30 
	--cluster-config ../cluster.json 
	--cluster 'bsub -q {cluster.queue}'
	-n {cluster.n} 
	-R "rusage[mem={cluster.memory}]" 
	-M {cluster.memory} 
	-o ./cis_eqlts.20181120.o 
	-e ./cis_eqlts.20181120.e
	
	
GENPATH=/home/jm2294/projects/RNAseq/hipsci_pipeline/geuvadis_CEU_test_data/Genotypes	
PHEPATH=/home/jm2294/projects/RNAseq/hipsci_pipeline/geuvadis_CEU_test_data/Expression
OUTPATH=/home/jm2294/projects/RNAseq/hipsci_pipeline/geuvadis_CEU_test_data/output_test
	
python /home/jm2294/projects/RNAseq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen ${GENPATH}/Geuvadis\
 -af ${PHEPATH}/Geuvadis_CEU_Annot.txt\
 -pf ${PHEPATH}/Geuvadis_CEU_YRI_Expr.txt.gz\
 -od ${OUTPATH}\
 -cf ${PHEPATH}/Geuvadis_CEU_YRI_covariates.txt\
 -c\
 --kinship_file ${GENPATH}/Geuvadis_chr1_kinship.normalized.txt\
 -np 10000\
 -maf 0.000000001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w 250000\
 --block_size 1500\
 --sample_mapping_file /home/jm2294/projects/RNAseq/hipsci_pipeline/geuvadis_CEU_test_data/Geuvadis_CEU_gte.txt
 
 -gr {chunkFull}\
 
 python /home/jm2294/projects/RNAseq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --plink ${GENPATH}/Geuvadis\
 -af ${PHEPATH}/Geuvadis_CEU_YRI_Expr.txt\
 -pf ${PHEPATH}/Geuvadis_CEU_Expr.txt\
 -od ${OUTPATH}\
 -c\
# -cf ${PHEPATH}/Geuvadis_CEU_YRI_covariates.txt\
 --kinship_file ${GENPATH}/Geuvadis_chr1_kinship.normalized.txt\
 -np 10000\
 -maf 0.000000001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w 250000\
 --block_size 1500\
 --sample_mapping_file /home/jm2294/projects/RNAseq/hipsci_pipeline/geuvadis_CEU_test_data/Geuvadis_CEU_gte.txt

### test on own data. 
source activate limix_qtl
GENPATH=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq
PHEPATH=/home/jm2294/projects/RNAseq/test_run
OUTPATH=/home/jm2294/projects/RNAseq/test_run/output
	
python /home/jm2294/projects/RNAseq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen ${GENPATH}/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile\
 -af ${PHEPATH}/Feature_Annotation_Ensembl_gene_ids_autosomes.txt\
 -pf ${PHEPATH}/phenotype_5281-fc-genecounts.txt\
 -od ${OUTPATH}\
 -t\
 --sample_mapping_file ${PHEPATH}/sample_mapping_file_gt_to_phe.txt\
 #-cf ${PHEPATH}/Geuvadis_CEU_YRI_covariates.txt\
 #--kinship_file ${GENPATH}/Geuvadis_chr1_kinship.normalized.txt\
 -np 10000\
 -maf 0.000000001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w 250000\
 --block_size 1500\
 
source activate limix_qtl
GENPATH=/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq
PHEPATH=/home/jm2294/projects/RNAseq/test_run
OUTPATH=/home/jm2294/projects/RNAseq/test_run/output
	
python /home/jm2294/projects/RNAseq/hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py\
 --bgen ${GENPATH}/impute_22_23500000-24500000_interval_RNAseq_batch1_withsamples_testfile_uniqueRSids\
 -af ${PHEPATH}/Feature_Annotation_Ensembl_gene_ids_autosomes.txt\
 -pf ${PHEPATH}/phenotype_5281-fc-genecounts.txt\
 -od ${OUTPATH}\
 --sample_mapping_file ${PHEPATH}/sample_mapping_file_gt_to_phe.txt\
 -t\
 -np 10000\
 -maf 0.000000001\
 -hwe 0.00001\
 -cr 0.95\
 -gm standardize\
 -w 250000\
 --block_size 50\
 #-c\
 #-cf ${PHEPATH}/Geuvadis_CEU_YRI_covariates.txt\
 #--kinship_file ${GENPATH}/Geuvadis_chr1_kinship.normalized.txt\