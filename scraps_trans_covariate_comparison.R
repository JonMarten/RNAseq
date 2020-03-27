setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/03_tensorqtl/results/python_module_method")
library(dplyr)
library(data.table)

resBase <- fread(data.table = F, "tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_merged_annotated_transSNPs_eSNPs_sentinels.csv")
resPeer <- fread(data.table = F, "tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_PEER20_merged_annotated_transSNPs_eSNPs_sentinels.csv")
resCell <- fread(data.table = F, "tensorqtl_trans_MAF0.005_age_sex_rin_batch_readDepth_PC10_NeutPCT_LympPCT_MonoPCT_EoPCT_BasoPCT_merged_annotated_transSNPs_eSNPs_sentinels.csv")

library(VennDiagram)

egBase <- unique(resBase$gene_name)
egCell <- unique(resCell$gene_name)
egPeer  <- unique(resPeer$gene_name)

venn.diagram(list(egBase, egCell, egPeer), category.names = c("Base", "Cell Percentages", "PEER"), filename = "transQTL_eGenes_venndiagram.png")

notPeerBase <- egBase[which(!egBase %in% egPeer)]
notPeerCell <- egCell[which(!egCell %in% egPeer)]
notPeerAll <- union(notPeerBase, notPeerCell)
peerOnly <- egPeer[which(!(egPeer %in% egCell | egPeer %in% egBase))]

basesig <- union(union(egBase, egCell), egPeer)

#enrichment
baseset <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/phenotype/INTERVAL_RNAseq_phase1_filteredSamplesGenes_TMMNormalised_FPKM_Counts_foranalysis.txt", data.table = F)
basegenes <- baseset[,1]
write.table(basegenes,"enrichment_base_set.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(notPeerBase, "enrichment_baseNotPeer.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(notPeerCell, "enrichment_cellNotPeer.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(notPeerAll, "enrichment_anyNotPeer.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(peerOnly, "enrichment_PeerOnly.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(basesig, "enrichment_baseEgenesOnly.txt", sep = "\t", quote = F, row.names = F, col.names = F)