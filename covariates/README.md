#### Derivation of coviarate file
Covariate data is provided for the whole cohort by the data management team; the current data relase is stored in `V:\INTERVAL\INTERVAL_MasterMerge\p1103\20190531\Data`. Phenotypes are stored in `INTERVALdata_31MAY2019.csv` and `INTERVALdata_P3_31MAY2019.csv`. `omicsMap.csv` and `omicsMap_P3.csv` map the randomly-generated 7 digit individual identifiers used by the data management team to the identifiers used on the various omic datasets, including the RNA-seq data. 

The INTERVAL study is split into phases which do not directly affect the RNA seq data but do affect how phenotype data is stored. **Each individual has only one timepoint at which RNA was measured.** This is preferentially two years after recruitment but may be later if no blood was available at this timepoint.

[This script](covariates/02_make_id_mapper_file.R) processes the omicsmap files to create a single file, `rna_id_mapper.csv`. This maps the phenotype '**identifier**' to the INT_RNA identifier used for the RNA seq files. It also lists '**batch**' (corresponding to RNA seq data release batch from the Sanger) and '**phase**', corresponding to the phase of the sample used for sequencing. This is obtained from the omicsmap file - the column containing the INT_RNA identifier specifies the phase. Where an id is included in multiple columns it is assigned one on the following order of preference: 24m > 48m > p3. Phase information is used to calculate the correct age for the RNA seq blood sample.

#### Calculating Ages
Ages are calculated using date of birth (assuming the 15th of the month) and the appointment date corresponding to blood sample, as obtained above. The plots below show how age varies with phase and batch, both absolute age and time between initial appointment and the blood sample used for RNA. Note that all batch 1 samples (the pilot study) are in the 24 month phase. Currently, phentoype data is only mapped to INT_RNA identifiers for batches 1 and 2.

![ageplot](https://github.com/JonMarten/RNAseq/blob/master/covariates/bothplot.png?raw=true)
