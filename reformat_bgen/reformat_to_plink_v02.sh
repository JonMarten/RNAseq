module load ceuadmin/plink/2.0_09_09_18

cd /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/bgen_b38_filtered

plink2 \
 --bgen impute_22_interval_b38_filtered.bgen\
 --make-bed\
 --sort-vars\
 --out plink/impute_22_interval_b38_filtered