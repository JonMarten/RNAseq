library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
theme_set(theme_cowplot())

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/mediation")


#soma <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/somalogic_proteomics/gwasqc/somalogic_qcgwas_4000.csv", data.table = F)
#soma_meta <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/somalogic_proteomics/Documents/Somamer_info.tsv", data.table = F)
#soma_meta <- soma_meta %>% 
#  mutate(colname = gsub("\\.", "", tolower(SOMAMER_ID)))
#
#uniprots <- soma_meta$UniProt
#write.table(quote = F, col.names = F, row.names = F, uniprots, file = "soma_uniprotID.txt")
## this was pasted into the converter at https://www.uniprot.org/uploadlists/ and the results pasted back into the file.

#uniprots <- fread(data.table = F, "soma_uniprotID.txt")
#unmapped <- fread(data.table = F, "soma_uniprot_unmapped.txt")

#eneexp <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/gene_expr_filtered_TMM_rpkms_genesI#vRankTransf_2747Samp_18373genes.csv", data.table = F)

#eneexp_unf <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/globus_phase2_recalled/results/combined/5591-star-fc-genecounts.txt", #ata.table = F)


# olink
#olink_panel <- fread("/rds/project/jmmh2/rds-jmmh2-projects/sle_grs/OLINK/olink_protein_panel_manifest.csv", data.table = F)
#library(stringr)
#onames <- str_split_fixed(olink_panel$protein, "___", 2)
#write.table(quote = F, col.names = F, row.names = F, onames[,2], file = "olink_uniprotID.txt")
#
#o_uniprots <- fread(data.table = F, "olink_uniprotID.txt")
#o_unmapped <- fread(data.table = F, "olink_uniprot_unmapped.txt")
#
## merge
#
#allENS <- union(union(uniprots$ENSEMBL, geneexp_unf$ENSEMBL_ID), o_uniprots$ENSEMBL)
#
#allENS_df <- data.frame("ENSEMBL_ID" = allENS,
#                        "RNAseq" = ifelse(allENS %in% geneexp_unf$ENSEMBL_ID, 1, 0),
#                        "RNAseq_filtered" = ifelse(allENS %in% geneexp$gene_id, 1, 0),
#                        "SOMA" = ifelse(allENS %in% uniprots$ENSEMBL, 1, 0),
#                        "OLINK" = ifelse(allENS %in% o_uniprots$ENSEMBL, 1, 0))
#
#vennlist = list("RNAseq" = allENS_df$ENSEMBL_ID[allENS_df$RNAseq == 1],
#                "RNAseq_filtered" = allENS_df$ENSEMBL_ID[allENS_df$RNAseq_filtered == 1],
#                "SOMA" = allENS_df$ENSEMBL_ID[allENS_df$SOMA == 1],
#                "Olink" = allENS_df$ENSEMBL_ID[allENS_df$OLINK == 1])
#
#
#library(VennDiagram)
#venn.diagram(vennlist, file = "venn.png")

# protein data
soma <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/somalogic_proteomics/gwasqc/somalogic_qcgwas_4000.csv", data.table = F)

# Merge olink data
olinkNeu <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_neu_qcgwas.csv", data.table = F)
names(olinkNeu)[-1] <- paste0("neu_", names(olinkNeu)[-1])
olinkCVD2 <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_qcgwas_cvd2.csv", data.table = F)
names(olinkCVD2)[-1] <- paste0("cvd2_", names(olinkCVD2)[-1])
olinkCVD3 <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_qcgwas_cvd3.csv", data.table = F)
names(olinkCVD3)[-1] <- paste0("cvd3_", names(olinkCVD3)[-1])
olinkINF <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/phenotype/olink_proteomics/gwasqc/olink_qcgwas_inf.csv", data.table = F)
names(olinkINF)[-1] <- paste0("inf_", names(olinkINF)[-1])

olink_merge <- full_join(by = "aliquot_id",olinkNeu,
                         full_join(by = "aliquot_id",olinkCVD2, 
                                   full_join(by = "aliquot_id",olinkCVD3, olinkINF)
                                   )
                         )
                                  
# scott's curated lists
scottSoma <- fread(data.table = F, "/rds/project/jmmh2/rds-jmmh2-projects/polygenic/internal/interval_grs_scan/analyses/processed_traits/somalogic_proteins/trait_info.tsv")
scottOlink <- fread(data.table = F, "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/geno_files/updated_olink_info.txt")
write.table(quote = F, col.names = F, row.names = F, scottOlink$UniProt, file = "olink_uniprotID_scott.txt")
# as above, manually convery to ENSEMBL ids. Required 4 to be looked up by hand:  c("P01137","Q8WWJ7","Q13541","Q9HAN9"). Backed up as olink_uniprotID_scott_backup.txt

scottOlinkNames <- fread("olink_uniprotID_scott_backup.txt", data.table = F)
names(scottOlinkNames) <- c("UniProt", "ENSEMBL")
scottOlink <- left_join(scottOlink, scottOlinkNames)

geneexp <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/gene_expr_filtered_TMM_rpkms_genesInvRankTransf_2747Samp_18373genes.csv", data.table = F)
geneexp_unf <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/globus_phase2_recalled/results/combined/5591-star-fc-genecounts.txt", data.table = F)

allENS <- union(union(union(scottSoma$Ensembl.Gene.ID, 
                      geneexp_unf$ENSEMBL_ID), 
                scottOlinkNames$To),geneexp$gene_id)

allENS_df <- data.frame("ENSEMBL_ID" = allENS,
                        "RNAseq" = ifelse(allENS %in% geneexp_unf$ENSEMBL_ID, 1, 0),
                        "RNAseq_filtered" = ifelse(allENS %in% geneexp$gene_id, 1, 0),
                        "SOMA" = ifelse(allENS %in% scottSoma$Ensembl.Gene.ID, 1, 0),
                        "OLINK" = ifelse(allENS %in% scottOlinkNames$ENSEMBL, 1, 0))

vennlist = list("RNAseq_phase2" = allENS_df$ENSEMBL_ID[allENS_df$RNAseq == 1],
                "RNAseq_filtered_phase1" = allENS_df$ENSEMBL_ID[allENS_df$RNAseq_filtered == 1],
                "SOMA" = allENS_df$ENSEMBL_ID[allENS_df$SOMA == 1],
                "Olink" = allENS_df$ENSEMBL_ID[allENS_df$OLINK == 1])

allENS_df <- allENS_df %>%
  mutate(SOMA_colname = ifelse(SOMA == 1, scottSoma$variable[match(ENSEMBL_ID, scottSoma$Ensembl.Gene.ID)], NA),
         OLINK_colname = ifelse(OLINK == 1, scottOlink$variable[match(ENSEMBL_ID, scottOlink$ENSEMBL)], NA))

library(VennDiagram)
venn.diagram(vennlist, file = "venn_scott.png")

# Rename all IDs to match Affy ID
omicsMap <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/processed/INTERVAL_omics_table_14MAY2020.csv", data.table = F)

geneexpT <- geneexp %>%
  select(-(gene_id:gene_length)) %>%
  t %>%
  data.frame
names(geneexpT) <- geneexp$gene_id
geneexpT <- geneexpT %>%
  mutate(affyID = omicsMap$affymetrix_ID[match(rownames(geneexpT), omicsMap$RNA_ID)]) %>%
  select(affyID, everything())

soma <- soma %>%
  mutate(affyID = omicsMap$affymetrix_ID[match(aliquot_id, omicsMap$soma_ID)]) %>%
  select(affyID, everything())

olink_merge <- olink_merge %>%
  mutate(affyID = omicsMap$affymetrix_ID[match(aliquot_id, omicsMap$olink_ID)]) %>%
  select(affyID, everything())

## Extract protein and gene expression for a named gene
pullProt <- function(gene){
  generow <-  allENS_df %>% filter(ENSEMBL_ID == gene)
  out <- geneexpT %>% select(affyID, eval(gene))
  if(generow$SOMA == 1) {
    somacol <- generow$SOMA_colname
    somaProt <- soma %>% select(affyID, eval(somacol))
    names(somaProt)[2] <- paste0("soma__", names(somaProt)[2])
    out <- full_join(out, somaProt, by = "affyID")
  }
  if(generow$OLINK == 1) {
    olinkcol <- generow$OLINK_colname
    olinkProt <- olink_merge %>% select(affyID, eval(olinkcol))
    names(olinkProt)[2] <- paste0("olink__", names(olinkProt)[2])
    out <- full_join(out, olinkProt, by = "affyID")
  }
  return(out)
}

ltbr <- pullProt("ENSG00000111321")
ggplot(ltbr, aes(x = ENSG00000111321, y = olink__cvd3_ltbr___p36941))

