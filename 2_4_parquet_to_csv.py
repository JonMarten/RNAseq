import pandas as pd
#df = pd.read_parquet("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results/tensorqtl_trans_test.trans_qtl_pairs.parquet")
#df.to_csv("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results/tensorqtl_trans_test.trans_qtl_pairs.csv", index=False)

dir = "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results/"
# cis files
for i in range(1,23):
  df = pd.read_parquet(dir + "tensorqtl_cis_MAF0.005_cisnominal.cis_qtl_pairs." + str(i) + ".parquet")
  df.to_csv(dir + "tensorqtl_cis_test.cis_qtl_pairs." + str(i) + ".csv", index=False)
