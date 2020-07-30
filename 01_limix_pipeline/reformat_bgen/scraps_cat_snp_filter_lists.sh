for i in {1..22}
do
	echo c${i}_b38_filter_snps.txt
	cat c${i}_b38_filter_snps.txt >> all_b38_filter_snps.txt
done