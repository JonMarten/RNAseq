library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
#library(hexbin)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/results/cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_5GenesPerChunk/processed")

prefix <- "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20"

eSNPs <- fread("cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs.txt", data.table = F)
eQTLgen <- fread("zcat /rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/eQTLgen/cis-eQTLs_full_20180905.txt.gz", data.table = F)

# merge with eQTLgen
eQTLgen <- eQTLgen %>%
  mutate(matchID = paste0(Gene, ":", SNP)) %>%
  select(matchID, Zscore.eqtlgen = Zscore, assessed_allele.eqtlgen = AssessedAllele, ref_allele.eqtlgen = OtherAllele, Ncohorts.eqtlgen = NrCohorts, N.eqtlgen = NrSamples, P.eqtlgen = Pvalue, FDR.eqtlgen = FDR)

eSNPs <- eSNPs %>%
  mutate(matchID = paste0(feature_id, ":", rsid))
merge <- left_join(eSNPs, eQTLgen)

# Flip alleles to match eQTLgen and check direction mismatch
flips <- merge$assessed_allele != merge$assessed_allele.eqtlgen & !is.na(merge$assessed_allele.eqtlgen)
merge <- merge %>%
  mutate(zscore = beta / beta_se) %>%
  mutate(Zscore.eqtlgen.flipped = ifelse(flips, -Zscore.eqtlgen, Zscore.eqtlgen),
         flipped = ifelse(flips, 1, 0)) %>%
  mutate(signMismatch = ifelse(sign(zscore) != sign(Zscore.eqtlgen.flipped), 1, 0))

# Write out merged results
fwrite(merge, file = "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_eQTLgenReplication.txt")

# check which of our eGenes are significant in eqtlgen
eGenes.rep <- merge %>% 
  filter(!is.na(FDR.eqtlgen)) %>%
  group_by(feature_id) %>%
  summarise(gene_name = paste(unique(gene_name), collapse = ","), minP.eQTLgen = min(FDR.eqtlgen, na.rm = T)) %>% 
  arrange(minP.eQTLgen) %>%
  data.frame() %>%
  mutate(significant = ifelse(minP.eQTLgen < 0.05 / length(n()), 1, 0))

# Check which of our eSNPs are signficant in eqtlgen
eSNP_replication <- merge %>%
  filter(!is.na(FDR.eqtlgen)) %>%
  mutate(P.BH.eqtlgen = p.adjust(P.eqtlgen, method = "BH")) %>%
  mutate(sigrep = ifelse(P.BH.eqtlgen < 0.05 & signMismatch == 0, 1, 0))
# correlate Z-scores
cor(eSNP_replication$zscore, eSNP_replication$Zscore.eqtlgen.flipped)
cor(eSNP_replication$zscore[eSNP_replication$sigrep == 1], eSNP_replication$Zscore.eqtlgen.flipped[eSNP_replication$sigrep == 1])

# Write out merged eQTLgen replication
fwrite(eSNP_replication, file = "cis_eqtls_18373genes_age_sex_rin_batch_PC10_PEER20_eSNPs_eQTLgenReplication_reponly.txt")

# Pull genes with mismatched SNP Zscores
mismatchGenes <- eSNP_replication %>%
  filter(signMismatch == 1) %>%
  pull(feature_id) %>%
  unique()

# Classify SNPs as strand-ambiguous
eSNP_replication$feature_length <- eSNP_replication$feature_end - eSNP_replication$feature_start
eSNP_replication <- eSNP_replication %>%
  mutate(ambig = ifelse(assessed_allele == "A" & ref_allele == "T" |
                          assessed_allele == "T" & ref_allele == "A" |
                          assessed_allele == "C" & ref_allele == "G" |
                          assessed_allele == "G" & ref_allele == "C", 1, 0)) %>%
  mutate(veryambig = ifelse(ambig == 1 & maf > 0.45, 1, 0))

# Sort eSNP list 
eSNP_replication <- eSNP_replication %>% 
  arrange(feature_chromosome, feature_start, snp_position)

# get list of features with more than 2 significant SNPs mismatched  
mismatches <- eSNP_replication %>%
  filter(P.BH.eqtlgen < 0.05 & signMismatch  == 1) %>%
  group_by(feature_id) %>%
  summarise(n_mismatched = n(), chr = unique(feature_chromosome)) %>%
  data.frame 

flexmismatchGenes <- mismatches %>%
  filter(n_mismatched > 2) %>%
  pull(feature_id)

#Plot Z scores

g <- ggplot(eSNP_replication, aes(x = zscore, y = Zscore.eqtlgen.flipped)) + 
  geom_bin2d(bins = 100) +
  theme(aspect.ratio = 1, legend.position = "right") +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  scale_fill_continuous(type = "viridis") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") 

outname_plot <- paste0(prefix, "_eSNPs_eqtlgenReplication_zscoreplot.png")
save_plot(file = outname_plot, g, base_height = 10, base_width = 10)

# plot single gene
onegene <- ggplot(filter(eSNP_replication, feature_id == "ENSG00000229391"), 
             aes(x = zscore, y = eQTLgen_Zscore_flipped, color = Significant)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 0.8, alpha = 1, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
  theme_minimal() +
  scale_colour_manual(values = c("gray75", "orangered", "chartreuse3")) +
  facet_wrap(. ~ gene_name)  


# plot mismatches coloured by significance 
g2 <- ggplot(filter(eSNP_replication, feature_id %in% flexmismatchGenes), 
             aes(x = zscore, y = eQTLgen_Zscore_flipped, color = Significant)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 0.8, alpha = 1, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
  theme_minimal() +
  scale_colour_manual(values = c("gray75", "orangered", "chartreuse3")) +
  facet_wrap(. ~ gene_name)  

outname_plot2 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch.png")
save_plot(file = outname_plot2, g2, base_height = 7, base_width = 12)

# Plot mistmatches coloured by strand-ambiguity
g3 <- ggplot(filter(eSNP_replication, feature_id %in% flexmismatchGenes), 
             aes(x = zscore, y = eQTLgen_Zscore_flipped, color = as.factor(ambig))) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 0.8, alpha = 1, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
  theme_minimal() +
  scale_colour_manual(values = c("gray75", "orangered", "chartreuse3")) +
  facet_wrap(. ~ gene_name)  

outname_plot3 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch_ambig.png")
save_plot(file = outname_plot3, g3, base_height = 7, base_width = 12)

# Plot mismatches coloured by eQTLgen number of cohorts size
g4<- ggplot(filter(eSNP_replication, feature_id %in% flexmismatchGenes), 
            aes(x = zscore, y = eQTLgen_Zscore_flipped, color = eQTLgen_NrCohorts)) + 
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  geom_point(cex = 0.8, alpha = 1, pch = 20) +
  theme(aspect.ratio = 1, legend.position = "right") +
  labs(x = "Limix Z-score", y = "eQTLgen Z-score") +
  theme_minimal() +
  scale_colour_viridis_c() +
  facet_wrap(. ~ gene_name)  

outname_plot4 <- paste0(jobdf$job_id[i],"/processed/",jobdf$job_id[i], "_eSNPs_eqtlgenReplication_zscoreplot_mismatch_eQTLnCohorts.png")
save_plot(file = outname_plot4, g4, base_height = 7, base_width = 12)


fwrite(jobdf, file = "peer_comparison_summary.csv")

