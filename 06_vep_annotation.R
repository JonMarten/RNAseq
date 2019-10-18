# Generate query file for VEP web interface and merge into results file

library(data.table)
library(stringr)
library(dplyr)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/")

jobids <- c("cis_eqtls_18373genes_age_sex_rin_batch_PC10", "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20")
i <- 1
resdir <- paste0("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/", jobids[i],"/processed")
setwd(resdir)

c22 <- fread(paste0(jobids[i],"_eSNPs.txt"), data.table=F) %>%
  distinct()
# vep lookup including non-rsid snps?
c22$alleles <- str_extract(c22$snp_id, "[A-Z]_[A-Z]")
c22 <- c22 %>%
  mutate(alleles = str_extract(snp_id, "[A-Z]+_[A-Z]+"))
almat <-  str_split_fixed(c22$alleles, "_", 2)  
c22$ref_allele <- ifelse(almat[,1] == c22$assessed_allele, almat[,2], almat[,1])


c22b <- c22 %>% mutate(type = ifelse(nchar(ref_allele) > nchar(assessed_allele), 
                             "DEL", 
                             ifelse(nchar(assessed_allele) > nchar(ref_allele), 
                                    "INS", 
                                    ifelse(nchar(assessed_allele) == nchar(ref_allele), 
                                           "SNP", 
                                           NA
                                           )
                                    )
                             )
                       )




vepQuery <- c22 %>%
  select(chromosome = snp_chromosome, start = snp_position, ref_allele, assessed_allele, strand  = feature_strand, identifier = rsid) %>%
  mutate(end = start, allele = paste0(ref_allele, "/", assessed_allele)) %>%
  select(-ref_allele, -assessed_allele) %>%
  select(chromosome, start, end, allele, strand, identifier) %>%
  arrange(chromosome, start) %>%
  filter(!is.na(chromosome) & !is.na(start) & !is.na(end))

fwrite(vepQuery, file = "results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_vepQuery.txt", sep = "\t", col.names = F)

rsidList <- vepQuery %>% pull(identifier)
write.table(rsidList, file = "results_merged_chr22_window500000_Perm100_MAF0.01_eSNPs_vepQuery_rsids.txt", sep = "\t", col.names = F, row.name = F, quote = F)


# Merge in VEP results          

vep <- fread("vep.txt", data.table = F)
names(vep)[1] <- "rsid"         

vep %>% 
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
         
         
         
         