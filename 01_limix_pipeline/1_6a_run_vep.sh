#!/bin/bash
#SBATCH -A PAUL-SL3-CPU
#SBATCH -p skylake-himem
#SBATCH --mem 150G
#SBATCH --job-name=vep
#SBATCH --time=12:0:0
#SBATCH --output=/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/logs/run_vep_%A.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

. /etc/profile.d/modules.sh     
module purge
module load ceuadmin/tabix/0.2.6

#/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20/processed/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_vepQuery_rsids
INFILE=${1}.txt
OUTFILE=${1}_vepoutput.txt

vep\
 --force_overwrite\
 --cache\
 --dir_cache /rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/vep\
 --everything\
 -i $INFILE
 -o $OUTFILE
