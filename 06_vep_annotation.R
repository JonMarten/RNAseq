# Generate query file for VEP web interface and merge into results file
# I'll be damned if I can figure out how to get the VEP to run on CSD3 without crashing, so it has to be done in this clumsy way. 

library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(UpSetR)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_5GenesPerChunk/processed")
prefix <- "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20"

eSNPs <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs.txt", data.table = F)

vep1 <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_rsidsForVEP_part1_results.txt", data.table = F)
vep2 <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_rsidsForVEP_part2_results.txt", data.table = F)

vep <- rbind(vep1, vep2)
names(vep)[1] <- "rsid"
#m <- left_join(eSNPs, vep, by = "rsid")

vepTrim2 <- vep %>%
 # mutate(IMPACT = factor(IMPACT, levels = c("HIGH", "MODERATE","LOW","MODIFIER","-"))) %>%
  group_by(rsid) %>%
  arrange(Feature) %>%
  filter(row_number() == 1) %>%
  data.frame


# Merge in VEP results          

m <- left_join(eSNPs, vepTrim2, by = "rsid")
fwrite(m, file = "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_VEPannotated.csv")
# m <- fread(data.table = F, "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_VEPannotated.csv")

# Significant eSNPs per eGene
eSNPs <- eSNPs %>%
  mutate(geneLength = feature_end - feature_start,
         TSSdist = snp_position - feature_start,
         positionAlongGene = (snp_position - feature_start)/geneLength)

eSNPs$HLA <- grepl("HLA", eSNPs$gene_name)
  
snpCount <- eSNPs %>%
  group_by(feature_id) %>%
  summarise(geneLength = unique(geneLength),
            nSNPs_significant = n(),
            nSNPs_perlength = n() / unique(geneLength),
            HLA = unique(HLA)) %>%
  data.frame()

ggplot(snpCount, aes(x = geneLength, y = nSNPs_significant, colour = HLA)) + 
  geom_point() +
  labs(x = "Gene length", y = "Number of significant eSNPs") 

# Density plot of number of significant eSNPs per gene
g1 <- ggplot(snpCount, aes(x = nSNPs_significant)) +
  geom_density()

# Density plot of number of significant eSNPs per gene normalised by gene length
g2 <- ggplot(filter(snpCount, nSNPs_perlength < 1), aes(x = nSNPs_perlength)) +
  geom_density()

plot_grid(g1, g2, ncol = 1)

# Density plot of number of TSS distance
g3 <- ggplot(eSNPs, aes(x = TSSdist)) +
  geom_density()

# Density plot of number of position of SNP along gene
g4 <- ggplot(eSNPs, aes(x = positionAlongGene)) +
  geom_density()

ggplot(eSNPs, aes(x = BIOTYPE)) + geom_bar()

plot_grid(g3, g4, ncol = 1)

b <- m %>% 
  group_by(Consequence) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(Consequence = factor(Consequence, levels = unique(Consequence))) 


  ggplot(b[-1,], aes(x = Consequence, y = n, fill = Consequence)) + 
    geom_bar(stat="identity") +
    coord_flip() +
    theme(legend.position = "none")
  
  tab <- m %>%
    group_by(Consequence) %>%
    summarise(count = n()) %>%
    data.frame() %>%
    mutate(Consequence = as.factor(Consequence)) %>% 
    arrange(desc(count))
  
  tab$Consequence <- gsub(pattern = ",", "&", tab$Consequence)
  
  # Make upset plot for consequences 
  a <- c(tab$count)
  names(a) <- tab$Consequence
 
   upset(fromExpression(a), 
        nsets = 15,
        order.by = "freq",
        mb.ratio = c(0.60, 0.40),
        text.scale = 1.5)
  
  
  


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
         
         
         
         