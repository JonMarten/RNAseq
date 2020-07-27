# Get RNA and protein levels for GDF15
# RNA/protein levels for soma/olink
# Stratify by SNP genotype
setwd("/rds/user/jm2294/rds-jmmh2-projects/interval_rna_seq/GDF15")

lapply(c("dplyr", "ggplot2", "rbgen", "data.table", "cowplot", "tidyr", "ggsignif"), require, character.only = T)
theme_set(theme_cowplot())
snps <- c("rs16982345","rs1054221","rs1059369","rs189593084") 
bgen <- "/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/GENETIC_DATA/b37_b38_liftover/b38_bgen/filtered/impute_19_interval_b38_filtered.bgen"

gene <- "ENSG00000130513"
prot <- "GDF15"
chr <- 19

olink <- fread("/rds/project/jmmh2/rds-jmmh2-projects/sle_grs/OLINK/INTERVAL_Olink_inf1_cvd2_cvd3_merged.csv", data.table = F)
soma <- fread("/rds/project/jmmh2/rds-jmmh2-projects/sle_grs/Somalogic_all_proteins_agesexrntrans.csv", data.table = F)
gdf15.s <- soma %>% select(id, "GDF15.soma" = GDF15)
#gdf15.o <- olink %>% select(id, age, sex = sexPulse, "GDF15.olink" = GDF.15___Q99988) # no effect of age and sex so not adjusting
gdf15.o <- olink %>% select(id, "GDF15.olink" = GDF.15___Q99988)

gdf15.prot <- full_join(gdf15.o, gdf15.s) %>%
  mutate(id = as.character(id))

geno <- bgen.load(bgen, rsids = snps)

bgenToDose <- function(bgen){
  nsnps <- nrow(bgen$variants)
  nids <- length(bgen$samples)
  out <- matrix(nrow = nids, ncol = nsnps)
  rownames(out) <- bgen$samples
  colnames(out) <- rownames(bgen$variants)
  for(i in 1:nsnps){
    out[,i] <- bgen$data[i,,2] + 2*bgen$data[i,,3]
  }
  out2 <- list("variants" = bgen$variants, "dosage" = out)
  return(out2)
}

dose <- bgenToDose(geno)
dosage <- dose[[2]] %>% data.frame
dosage$id <- rownames(dosage)

doseProt <- full_join(dosage, gdf15.prot)

# Get RNA seq data
rna <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/gene_expr_filtered_TMM_rpkms_genesInvRankTransf_2747Samp_18373genes.csv", data.table = F)
rnacov <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates/INTERVAL_RNAseq_phase1_age_sex_rin_batch_PC10_PEER20.txt", data.table = F)

gdf15.rna <- rna %>% 
  filter(gene_id == gene) %>%
  select(-c(gene_id:gene_length)) %>%
  t %>% 
  data.frame
names(gdf15.rna) <- "GDF15.rna"
gdf15.rna$phenotype_individual_id <- rownames(gdf15.rna)

namemap <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)
gdf15.rna <- full_join(namemap, gdf15.rna) %>%
  filter(!is.na(GDF15.rna))
names(gdf15.rna)[1] <- "id"
gdf15.rna$id <- as.character(gdf15.rna$id)

# Adjust RNA for age + sex + RIN + batch
rnacov <- rename(rnacov, phenotype_individual_id = sample_id)

gdf15.rna <- full_join(gdf15.rna, rnacov)
rna.lm <- lm(data = gdf15.rna, GDF15.rna ~ RIN + age_RNA + sex + batch1 + batch2 + batch3 + batch4 + batch5 + batch6+ batch7)

gdf15.rna <- gdf15.rna %>% 
  mutate(GDF15.rna.res = rna.lm$residuals) %>%
  select(id, phenotype_individual_id, GDF15.rna.res)

gdf15 <- full_join(gdf15.rna, doseProt)

gdf15 <- gdf15 %>%
  filter(!is.na(GDF15.rna.res))

gdf15 <- gdf15 %>%
  mutate(rs189593084.hardcall = ifelse(rs189593084 >= 1.5, "AA", 
                                       ifelse(rs189593084 < 0.5, "CC", "CA")),
         rs1059369.hardcall = ifelse(rs1059369 >= 1.5, "AA", 
                                       ifelse(rs1059369 < 0.5, "TT", "TA")),
         rs1054221.hardcall = ifelse(rs1054221 >= 1.5, "CC", 
                                       ifelse(rs1054221 < 0.5, "TT", "CT")),
         rs16982345.hardcall = ifelse(rs16982345 >= 1.5, "AA", 
                                       ifelse(rs16982345 < 0.5, "GG", "GA"))) 

g1 <- ggplot(gdf15, aes(x = rs189593084.hardcall, y = rs189593084)) + geom_boxplot()
g2 <- ggplot(gdf15, aes(x = rs1059369.hardcall, y = rs1059369)) + geom_boxplot()
g3 <- ggplot(gdf15, aes(x = rs1054221.hardcall, y = rs1054221)) + geom_boxplot()
g4 <- ggplot(gdf15, aes(x = rs16982345.hardcall, y = rs16982345)) + geom_boxplot()

save_plot(plot_grid(g1,g2,g3,g4), filename = "hardcalls.png", base_height = 10, base_width = 10)

gdf15 <- gdf15 %>%
  select("geno.id" = id,
         "rna.id" = phenotype_individual_id,
         rs189593084:rs16982345,
         rs189593084.hardcall:rs16982345.hardcall,
         GDF15.rna.res,
         GDF15.olink,
         GDF15.soma)

fwrite(gdf15, file = "GDF15_olink_soma_rna_genotype.csv")

# Make plots

gdf15$rs189593084.hardcall <- factor(gdf15$rs189593084.hardcall, levels = c("AA","CA","CC"))

o1 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.olink, colour = rs189593084.hardcall)) +
  geom_point(alpha = 0.5, pch = 20) + 
  facet_wrap(.~rs189593084.hardcall, drop = F) + 
  coord_fixed()
o2 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.olink, colour =  rs1059369.hardcall)) + 
  geom_point(alpha = 0.5, pch = 20)+ 
  facet_wrap(.~rs1059369.hardcall)+ 
  coord_fixed()
o3 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.olink, colour = rs1054221.hardcall)) + 
  geom_point(alpha = 0.5, pch = 20)+ 
  facet_wrap(.~rs1054221.hardcall)+ 
  coord_fixed()
o4 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.olink, colour = rs16982345.hardcall)) + 
  geom_point(alpha = 0.5, pch = 20)+ 
  facet_wrap(.~rs16982345.hardcall)+ 
  coord_fixed()

save_plot(plot_grid(o1,o2,o3,o4), filename = "gdf15_olink_scatterplot.png", base_height = 5, base_width = 15)

s1 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.soma, colour = rs189593084.hardcall)) +
  geom_point(alpha = 0.5, pch = 20) + 
  facet_wrap(.~rs189593084.hardcall, drop = F) + 
  coord_fixed()
s2 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.soma, colour =  rs1059369.hardcall)) + 
  geom_point(alpha = 0.5, pch = 20)+ 
  facet_wrap(.~rs1059369.hardcall)+ 
  coord_fixed()
s3 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.soma, colour = rs1054221.hardcall)) + 
  geom_point(alpha = 0.5, pch = 20)+ 
  facet_wrap(.~rs1054221.hardcall)+ 
  coord_fixed()
s4 <- ggplot(gdf15, aes(x = GDF15.rna.res, y = GDF15.soma, colour = rs16982345.hardcall)) + 
  geom_point(alpha = 0.5, pch = 20)+ 
  facet_wrap(.~rs16982345.hardcall)+ 
  coord_fixed()

save_plot(plot_grid(s1,s2,s3,s4), filename = "gdf15_soma_scatterplot.png", base_height = 5, base_width = 15)

# Venn Diagram
library(VennDiagram)
sets <- list("RNA" = which(!is.na(gdf15$GDF15.rna)), 
             "Olink" = which(!is.na(gdf15$GDF15.olink)), 
             "SomaLogic" = which(!is.na(gdf15$GDF15.soma)))
venn.diagram(sets, filename = "venn.png")

# genotype box plots
gdf15.tall <- gdf15 %>% 
  mutate(rs189593084.hardcall = as.character(rs189593084.hardcall)) %>%
  select(-c(rs189593084:rs16982345)) %>%
  gather(-geno.id, -rna.id, -c(GDF15.rna.res:GDF15.soma), key = "SNP", value = "genotype") %>%
  mutate(SNP = gsub(pattern = ".hardcall", replacement = "", SNP)) %>%
  gather(-geno.id, -rna.id, -SNP, -genotype, key = "phenotype", value = "value") %>%
  mutate(phenotype = gsub("GDF15.","",phenotype) %>% toupper) %>%
  mutate(phenotype = factor(phenotype, levels = c("RNA.RES","OLINK","SOMA")))

boxplots <- ggplot(gdf15.tall, aes(x = genotype, y = value)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0)
  facet_wrap(. ~ phenotype + SNP, scales = "free") +
  geom_smooth(method = "lm", se = F, aes(group = 1))
save_plot(boxplots, filename = "gdf15_boxplots.png", base_height = 10, base_width = 15)


dat <- gdf15.tall %>% 
  filter(!is.na(value))

# linear models
lm(data = gdf15, GDF15.rna.res ~ rs1054221.hardcall) %>% summary
lm(data = gdf15, GDF15.olink ~ rs1054221.hardcall) %>% summary
lm(data = gdf15, GDF15.soma ~ rs1054221.hardcall) %>% summary
lm(data = gdf15, GDF15.rna.res ~ rs1059369.hardcall) %>% summary
lm(data = gdf15, GDF15.olink ~ rs1059369.hardcall) %>% summary
lm(data = gdf15, GDF15.soma ~ rs1059369.hardcall) %>% summary
lm(data = gdf15, GDF15.rna.res ~ rs16982345.hardcall) %>% summary
lm(data = gdf15, GDF15.olink ~ rs16982345.hardcall) %>% summary
lm(data = gdf15, GDF15.soma ~ rs16982345.hardcall) %>% summary
lm(data = gdf15, GDF15.rna.res ~ rs189593084.hardcall) %>% summary
lm(data = gdf15, GDF15.olink ~ rs189593084.hardcall) %>% summary
lm(data = gdf15, GDF15.soma ~ rs189593084.hardcall) %>% summary


lm(data = gdf15, GDF15.rna.res ~ rs1054221) %>% summary
lm(data = gdf15, GDF15.olink ~ rs1054221) %>% summary
lm(data = gdf15, GDF15.soma ~ rs1054221) %>% summary
lm(data = gdf15, GDF15.rna.res ~ rs1059369) %>% summary
lm(data = gdf15, GDF15.olink ~ rs1059369) %>% summary
lm(data = gdf15, GDF15.soma ~ rs1059369) %>% summary
lm(data = gdf15, GDF15.rna.res ~ rs16982345) %>% summary
lm(data = gdf15, GDF15.olink ~ rs16982345) %>% summary
lm(data = gdf15, GDF15.soma ~ rs16982345) %>% summary
lm(data = gdf15, GDF15.rna.res ~ rs189593084) %>% summary
lm(data = gdf15, GDF15.olink ~ rs189593084) %>% summary
lm(data = gdf15, GDF15.soma ~ rs189593084) %>% summary

############ 1
rb1 <- ggplot(gdf15, aes(x = gdf15$rs1059369.hardcall, y = GDF15.rna.res)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(cex = 0.3, comparisons = list(c("AA","TA"),c("TA","TT"),c("AA","TT")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

ob1 <- ggplot(gdf15, aes(x = gdf15$rs1059369.hardcall, y = GDF15.olink)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons = list(c("AA","TA"),c("TA","TT"),c("AA","TT")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

sb1 <- ggplot(gdf15, aes(x = gdf15$rs1059369.hardcall, y = GDF15.soma))  + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons = list(c("AA","TA"),c("TA","TT"),c("AA","TT")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
############ 2
rb2 <- ggplot(gdf15, aes(x = gdf15$rs1054221.hardcall, y = GDF15.rna.res)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0)+
  geom_signif(comparisons =  list(c("CC","CT"),c("CT","TT"),c("CC","TT")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

ob2 <- ggplot(gdf15, aes(x = gdf15$rs1054221.hardcall, y = GDF15.olink)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("CC","CT"),c("CT","TT"),c("CC","TT")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

sb2 <- ggplot(gdf15, aes(x = gdf15$rs1054221.hardcall, y = GDF15.soma))  + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("CC","CT"),c("CT","TT"),c("CC","TT")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

############ 3
rb3 <- ggplot(gdf15, aes(x = gdf15$rs16982345.hardcall, y = GDF15.rna.res)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("AA","GA"),c("GA","GG"),c("AA","GG")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

ob3 <- ggplot(gdf15, aes(x = gdf15$rs16982345.hardcall, y = GDF15.olink)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("AA","GA"),c("GA","GG"),c("AA","GG")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

sb3 <- ggplot(gdf15, aes(x = gdf15$rs16982345.hardcall, y = GDF15.soma))  + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("AA","GA"),c("GA","GG"),c("AA","GG")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))
#### 4
rb4 <- ggplot(gdf15, aes(x = gdf15$rs189593084.hardcall, y = GDF15.rna.res)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("CA","CC")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

ob4 <- ggplot(gdf15, aes(x = gdf15$rs189593084.hardcall, y = GDF15.olink)) + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("CA","CC")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))

sb4 <- ggplot(gdf15, aes(x = gdf15$rs189593084.hardcall, y = GDF15.soma))  + 
  geom_jitter(width = 0.2, alpha = 0.5, colour = "gray80") + 
  geom_boxplot(alpha = 0) +
  geom_signif(comparisons =  list(c("CA","CC")), map_signif_level = F, y_position = c(4,4.4,4.8), tip_length = 0.01) +
  geom_smooth(method = "lm", se = F, aes(group=1))


ugly <- plot_grid(rb1, rb2, rb3, rb4, ob1, ob2, ob3, ob4, sb1, sb2, sb3, sb4)
save_plot(ugly, filename = "gdf15_boxplots_ugly.png", base_height = 10, base_width = 15)







ggplot(gdf15, aes(x = rs189593084.hardcall, y = GDF15.rna.res)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("CA","CC")), map_signif_level = T) +
  geom_smooth(method = "lm", se = F, aes(group=1))

