# Reformat covariates for LIMIX pipeline
library(dplyr)
library(data.table)
setwd("U:/Projects/RNAseq/covariates")

dat <- fread("data_release_20190815/INTERVALdata_15AUG2019.csv", data.table = F)
ids <- fread("data_release_20190815/omicsMap.csv", data.table = F)
ids3 <- fread("data_release_20190815/omicsMap_P3.csv", data.table = F)
ids[which(ids == "", arr.ind = T)] <- NA
ids3[which(ids3 == "", arr.ind = T)] <- NA
allids <- full_join(ids, ids3)
#dat2 <- full_join(allids, dat)

# Read in ids from RNA seq data
rna_seq_ids <- fread("../globus/RNA_seq_ids.csv", data.table = F)

# filter to just RNA id columns and rows with an RNA id
RNAcols <- grep("RNA", names(allids))
rnaids <- allids %>% 
  select(identifier, RNAcols, attendanceDate_p3) %>%
  mutate(RNAseq_gwasQC_24m = as.character(RNAseq_gwasQC_24m),
         RNAseq_gwasQC_48m = as.character(RNAseq_gwasQC_48m),
         RNAseq_gwasQC_p3 = as.character(RNAseq_gwasQC_p3))

RNArows <- rowSums(is.na(select(rnaids, -identifier, -attendanceDate_p3))) != ncol(select(rnaids, -identifier, -attendanceDate_p3)) # list rows with at least one non-NA RNA id
rnaids <- rnaids %>% 
  filter(RNArows)

# Make single column for RNA identifier
RNA_any <- NA
for (i in 1:nrow(rnaids)) {
  idsb <- as.character(select(rnaids, -identifier, -attendanceDate_p3)[i,])
  RNA_any[i] <- paste(unique(na.exclude(idsb)), collapse = ",")
}
rnaids$RNA_any <- RNA_any

# Check overlap with RNA seq data identifiers
rnamiss <- which(!rna_seq_ids$RNA_id %in% rnaids$RNA_any)
phemiss <- which(!rnaids$RNA_any %in% rna_seq_ids$RNA_id )

if (length(rnamiss) == 0) {
  cat("\n No RNA seq ids missing in phenotype data.")
} else {
  cat("\n",
      length(rnamiss),
      " RNA seq ids missing in phenotype data:", 
      paste0(rna_seq_ids$RNA_id[rnamiss], collapse = ", "))
}

if (length(phemiss) == 0) {
  cat("\n No phenotype data ids missing in RNA seq.")
} else {
  cat("\n",
      length(phemiss),
      " phenotype data ids missing in RNA seq:", 
      paste0(rnaids$RNA_any[phemiss], collapse = ", "))
}

#
rnaidsPhase <- rnaids %>% select(identifier, RNAseq_RAW_24m, RNAseq_RAW_48m, RNAseq_RAW_p3, RNA_any,attendanceDate_p3) 
rnaidsPhase$RNAphase <- ifelse(!is.na(rnaidsPhase$RNAseq_RAW_24m), 
                     "24m",
                     ifelse(!is.na(rnaidsPhase$RNAseq_RAW_48m), 
                            "48m",
                            ifelse(!is.na(rnaidsPhase$RNAseq_RAW_p3), 
                                   "p3",
                                   "none")))

# Output mapper file
rna_id_mapper <- rnaidsPhase %>%
  select(identifier, RNA_id = RNA_any, phase = RNAphase, attendanceDate_p3) %>% 
  filter(!is.na(identifier))
rna_id_mapper <- full_join(rna_id_mapper, rna_seq_ids) #%>%
#  filter(!is.na(batch)) # Remove this filter to allow unnumbered batch 5 ids through for now
write.csv(rna_id_mapper, "rna_id_mapper.csv", quote = F, row.names = F)


# Check overlap between batches for different INTERVAL phases for Artika
all <- merge(dat,rnaids, all = T)
all <- merge(all, rna_id_mapper, all = T)
all %>% group_by(batch) %>%
  summarise(
    ids_RAW_24 = length(which(!is.na(RNAseq_RAW_24m))),
    ids_RAW_48 = length(which(!is.na(RNAseq_RAW_48m))),
    ids_RAW_p3 = length(which(!is.na(RNAseq_RAW_p3))))

all2 <- all %>% 
  select(RNAseq_RAW_24m, RNAseq_RAW_48m, RNAseq_RAW_p3, RNA_any, batch)
all2$num_time_points = apply(X = all2, MARGIN = 1, FUN = function(x){length(which(!is.na(x[1:3])))})

# Venn diagram showing overlap of phases
library(VennDiagram)

nonmiss <- list("24m" = which(!is.na(rnaids$RNAseq_RAW_24m)),
                "48m" = which(!is.na(rnaids$RNAseq_RAW_48m)),
                "p3" = which(!is.na(rnaids$RNAseq_RAW_p3)))
venn.diagram(nonmiss, filename = "venn.tiff")
                
