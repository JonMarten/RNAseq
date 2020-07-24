library(data.table)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(colorout)

setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates")

# Read in data, merge and convert to matrix
covs <- fread("processed/INTERVAL_RNA_batch1-12_master_covariates_release_2020_07_01.csv", data.table = F)
peer <- fread("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/peer_factors/peer_30Fact_100Iter/PEER_factors.txt", data.table = F)
names(peer)[1] <- "sample_id"
peer <- peer %>%
  select(sample_id, PEER1:PEER30)
  
allcovs <- inner_join(covs, peer) %>%
  select_if(is.numeric)

allcovs.m <- allcovs %>%
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

# calculate correlation matrix  
allcovs.m.cor <- allcovs.m %>% rcorr(type = "spearman")
r <- data.frame(allcovs.m.cor$r) 
row.names(r) <- names(r)
r$pheno = names(r)

# Filter to retain correlations above threshold thresh
corThresh <- 0.05 # SET THRESHOLD HERE
numPeers <- 30 # SET NUMBER OF PEER FACTORS HERE

rtall <- melt(r, id.vars = "pheno") %>%
  filter(grepl("PEER", pheno) & !grepl("PEER", variable) & !is.na(value) & !is.infinite(value))

strongcorrs <- rtall %>% 
  filter(pheno %in% paste0("PEER",1:numPeers) & abs(value) >= corThresh) %>% 
  pull(variable) %>% 
  unique %>%
  as.character

rtall <- rtall %>% 
  filter(variable %in% strongcorrs)

# Add in human readable names for blood cell traits  
bloodNames <- fread("/rds/project/jmmh2/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/covariates/blood_trait_names.csv", data.table = F)
names(bloodNames)[1] <- "variable"
rtall2 <- left_join(rtall, bloodNames, by = "variable")
rtall2 <- rtall2 %>%
  mutate(variable = ifelse(!is.na(humanName), humanName, ifelse(!is.na(Description), Description, variable))) #%>%
  #filter(is.na(Use) | Use == 1)
    
r2 <- reshape2::dcast(rtall2, pheno ~ variable, fun.aggregate = function(x) {max(x, na.rm = T)}) %>% 
    filter(pheno %in% paste0("PEER",c(paste0("0",1:9), as.character(10:numPeers))))
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
  rowAnno$class <- bloodNames$class[match(rownames(rowAnno),bloodNames$humanName)]
  
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
           height = 20,
           cellwidth = 30,
           cellheight = 30)
  
  
  

