# Phenotypes in the INTERVAL study
The INTERVAL study is split into phases which do not directly affect the RNA seq data but do affect how phenotype data is stored. **Each individual has only one timepoint at which RNA was measured.** This is preferentially two years after recruitment but may be later if no blood was available at this timepoint.

Non-omic phenotye data is provided for the whole cohort by the data management team; the current data relase is stored in `V:\INTERVAL\INTERVAL_MasterMerge\p1103\20190531\Data`. Phenotypes are stored in `INTERVALdata_31MAY2019.csv` and `INTERVALdata_P3_31MAY2019.csv`. These datasets use randomly-generated 7 digit identifiers for each individual - this is intended to prevent merging of phenotype data between different projects. The files `omicsMap.csv` and `omicsMap_P3.csv` map the random identifiers  identifiers used on the various omic datasets. 

# Generating the coviarate file
![Covariate generation flowchart](https://github.com/JonMarten/RNAseq/blob/master/covariates/cov_flowchart.png)

Three R scripts are used to generate the master covariate file. [The first](01_get_batch_ids.R) scans the RNA seq gene count files downloaded from the Sanger to get a list of `INT_RNA` ids for individuals with RNA seq data released, along with the batch in which they were sequenced. 

[The second](02_make_id_mapper_file.R) script combines this list with the two omicsMap files to create a file that maps the data management identifier to the `INT_RNA` id for all individuals with RNA seq data. This file also lists the '**phase**', which corresponds to the blood sample used for RNA seq. where an individual blood sample is assigned to multiple phases, they are assigned one on the following order of preference: 24m > 48m > p3.

[The third](03_make_master_RNA_covariate_file.R) script reads in the separate phenotype data files and selects the cell count data that corresponds to the blood sample used for RNA seq. It also calculates age at time of blood draw (assuming DOB as the 15th of the month, as more precise dates are not available). It then outputs a single column for each phenotype, appended with the suffix `___RNA`.

## Calculating Ages
The plots below show how age varies with phase and batch, both absolute age and time between initial appointment and the blood sample used for RNA. Note that all batch 1 samples (the pilot study) are in the 24 month phase. 
![ageplot](https://github.com/JonMarten/RNAseq/blob/master/covariates/bothplot.png?raw=true)
