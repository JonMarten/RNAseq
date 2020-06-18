library(data.table)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(viridis)

setwd("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/01_cis_eqtl_mapping/covariates")

  covs <- fread("INTERVAL_RNA_batch1-8_covariates_release_2019_08_15.csv", data.table = F)
  peer <- fread("PEER_factors_100Factors_plusCovariates_1000Iterations_AllBatches_TransformGenes_18373Genes.csv", data.table = F)
  
  names(peer)[1] <- "sample_id"
  
  namemap <- fread("../phenotype/sample_mapping_file_gt_to_phe_phase1.txt", data.table = F)
  names(namemap)[2] <- "sample_id"
  
  covs <- inner_join(namemap, covs)
  
  allcovs <- inner_join(covs, select(peer, -age_RNA), by = "sample_id")
  factors <- c("sequencingBatch", "sex","centre","batch1","batch2","batch3","batch4","batch5","batch6","batch7","batch8","sex2")
  allcovs <- allcovs %>%
    select(-factors)
  
  allcovs.m <- allcovs %>%
    #select(-(genotype_individual_id:Asset_Name), -(Instrument:Well), -(Extr_Sample_ID:Box_Position), -(Extraction_Date:Extraction_Time), -(FluidX_Plate_Position:Freezer_Shelf), -Agilent_RINe, -(Normalization_Plate_ID:Normalization_Plate_Position), -(Rack_ID:Plate_Position), -intervalPhase, -ethnicity, -ABORH, -(attendanceDate___RNA:processTime___RNA),-intercept) %>%
    select(-(genotype_individual_id : processTime___RNA), -(RIN:intercept)) %>%
    rename(PEER01 = PEER1, 
           PEER02 = PEER2, 
           PEER03 = PEER3, 
           PEER04 = PEER4, 
           PEER05 = PEER5, 
           PEER06 = PEER6, 
           PEER07 = PEER7, 
           PEER08 = PEER8, 
           PEER09 = PEER9) %>%
    as.matrix() %>%
    na.exclude() 
  
  allcovs.m.cor <- allcovs.m %>% rcorr(type = "spearman")
  
    r <- data.frame(allcovs.m.cor$r) 
    row.names(r) <- names(r)
    r$pheno = names(r)
    rtall <- melt(r, id.vars = "pheno")
    rtall <- rtall %>%
      filter(grepl("PEER", pheno) & !grepl("PEER", variable) & !is.na(value) & !is.infinite(value))
    
  strongcorrs <- rtall %>% filter(pheno %in% paste0("PEER",1:20) & abs(value) > 0.00) %>% pull(variable) %>% as.character
  
  rtall <- rtall %>% 
    filter(variable %in% strongcorrs)
  
    bloodNames <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/blood_trait_names.csv", data.table = F)
    names(bloodNames)[1] <- "variable"
    rtall2 <- left_join(rtall, bloodNames, by = "variable")
    rtall2 <- rtall2 %>% 
      mutate(variable = ifelse(!is.na(humanName), humanName, ifelse(!is.na(Description), Description, variable))) %>%
      filter(is.na(Use) | Use == 1)
    
    r2 <- dcast(rtall2, pheno ~ variable)
    
    r2 <- r2 %>% 
      filter(pheno %in% paste0("PEER",c(paste0("0",1:9), as.character(10:20))))
    row.names(r2) <- r2$pheno
    r2 <- select(r2, -pheno)
    r2[which(is.na(r2), arr.ind = T)] <- 0
    hmp <- pheatmap(t(r2), 
             display_numbers = T, 
             #cluster_rows = F, 
             cluster_cols = F, 
             na_col = "gray50",
             color = colorRampPalette(brewer.pal(name = "RdBu", n = 11))(50))
    
  plot(hmp$tree_row)
  abline(h=0.5, col="red", lty=2, lwd=2)
  
  groups <- cutree(hmp$tree_row, h = 0.5) 
  
  rowAnno <- data.frame(row.names = names(groups), "group" = as.factor(groups))

  rowAnno$cor_sum <- apply(X = t(r2), MARGIN = 1, FUN = function(x){sum(abs(x))})
  #rowAnno$cor_max <- apply(X = t(r2), MARGIN = 1, FUN = function(x){max(abs(x))})
  rowAnno$class <- bloodNames$class[match(rownames(rowAnno),bloodNames$humanName)]
  
  #rowAnno %>% 
  #  mutate(pheno = rownames(rowAnno)) %>%
  #  group_by(group) %>%
  #  summarise(n = n(), 
  #            maxcorsum = max(cor_sum),
  #            keepPheno_corSum = pheno[which.max(cor_sum)],
  #            maxcor = max(cor_max),
  #            keepPheno_maxcor = pheno[which.max(cor_max)])
  
  
  #pal <- colorRampPalette(brewer.pal(name = "RdBu", n = 11))(50)
  pal <- colorRampPalette(c("#087ccf","#ffffff","#cf1508"))(50)
  pheatmap(t(r2), 
           display_numbers = T,
           cluster_rows = T, 
           cluster_cols = F, 
           cutree_rows = 6,
           na_col = "gray50",
           annotation_row = select(rowAnno, -group),
           color = pal,
           treeheight_row = 200,
           filename = "blood_cell_traits_PEER_corplot.png",
           width = 20,
           height = 10,
           cellwidth = 30,
           cellheight = 30)
  
  
  
# Filter to include only most informative traits
#keep <- c("RIN","LYMPH%","NEUT%","WBC#","EO%","MONO%","Agilent_Conc_ng_ul","IRF","RET%","height","BASO%","PCT","age_RNA","ReadDepth")
#r3 <- t(r2)
#r3 <- r3[which(rownames(r3) %in% keep), ]
##r3[which(abs(r3) < 0.1, arr.ind = T)] <- 0
#pheatmap(r3, 
#         display_numbers = T, 
#         cluster_rows = T, 
#         cluster_cols = F
#         na_col = "gray50",
#         color = colorRampPalette(brewer.pal(name = "RdBu", n = 11))(50))
