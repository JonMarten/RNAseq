# Get SNP info from master bgen files
#SBATCH -p long
#SBATCH --job-name=update_bgen
#SBATCH -t 36:0:0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm2294@medschl.cam.ac.uk

# get start time
start=$(date +%s.%N)

module load gcc/5.2.0 
module load qctool2/rc4-6.8

module load gcc/5.2.0 
module load qctool2/rc4-6.8

qctool\
 -g /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_22_interval_b38_filtered.bgen\
 -s /home/jm2294/GENETIC_DATA/INTERVAL/master/impute_22_interval.bgen.sample\
 -excl-rsids /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/c22_b38_filter_snps_by_rsids.txt\
 -og /home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/impute_22_interval_b38_filtered_on_rsids.bgen
 
end=$(date +%s.%N)    
runtime=$(python -c "print(${end} - ${start})")

echo "Runtime was $runtime"