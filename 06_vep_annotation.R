# Generate query file for VEP web interface and merge into results file
# I'll be damned if I can figure out how to get the VEP to run on CSD3 without crashing, so it has to be done in this clumsy way. 

library(data.table)
library(stringr)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_5GenesPerChunk/processed")
prefix <- "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20"

eSNPs <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs.txt", data.table = F)

vep1 <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_rsidsForVEP_part1_results.txt", data.table = F)
vep2 <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_rsidsForVEP_part2_results.txt", data.table = F)

vep <- rbind(vep1, vep2)
names(vep)[1] <- "rsid"
m <- left_join(eSNPs, vep, by = "rsid")

# Merge in VEP results          

vep <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_eSNPs.txt_vepQuery_rsids_vepOfflineResults.txt", data.table = F)
names(vep)[1] <- "rsid"         
vepTrim <- vep[duplicated(vep$rsid),]
m <- left_join(eSNPs, vepTrim, by = "rsid")

m %>% 
  group_by(rsid) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  arrange(desc(count)) %>%
  data.frame() %>%
  head(30)

# Remove duplicate SNP/feature pairings (these arise due to overlapping windows)
# Identical rows are removed above with distinct(). Non-identical rows have different empirical pvalues.
# The row with the larger p-value is retained to be conservative.
c22 <- c22 %>%
  mutate(feat_snp_id = paste0(feature_id, "__", snp_id))
fsCount <- c22 %>% 
  group_by(feat_snp_id) %>%
  summarise(count = n(),
            minp = min(empirical_feature_p_value),
            maxp = max(empirical_feature_p_value)) %>%
  filter(count > 1) %>%
  arrange(desc(count))

c22_nodupe <- c22 %>%
  filter(!feat_snp_id %in% fsCount$feat_snp_id)

c22_dupe <- c22 %>%
  filter(feat_snp_id %in% fsCount$feat_snp_id) %>%
  mutate(matchp = fsCount$maxp[match(feat_snp_id, fsCount$feat_snp_id)]) %>%
  filter(empirical_feature_p_value == matchp) %>%
  dplyr::select(-matchp)

c22_filt <- rbind(c22_nodupe, c22_dupe) %>%
  arrange(feature_start, snp_position)

merge <- left_join(c22_filt, distinct(vep), by = "rsid") %>% distinct
fwrite(merge, "results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_vepannotated.csv")
  
## Annotate with intron number

introns <- str_split_fixed(merge$INTRON, "/", 2) %>%
  data.frame(stringsAsFactors = F)
names(introns) <- c("variant_intron_no", "total_no_introns")
introns <- introns %>%
  mutate(variant_intron_no = as.numeric(variant_intron_no),
         total_no_introns = as.numeric(total_no_introns))

merge <- cbind(merge, introns)      

# Get non-annotated SNPs
merge %>% filter(is.na(IMPACT))
         
         
         
         