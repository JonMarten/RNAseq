# Make new alternate ID for bgen 1.2 files
setwd("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq")
library(dplyr)
library(data.table)
a <- fread("c22_filtered_snp_stats.txt", skip = 8, data.table=F)
a <- a %>%
  mutate(indel = ifelse(nchar(alleleA)>1 | nchar(alleleB)>1, 1, 0))

# Create CPTIDs for SNPs (sort alleles in alphabetical order)
snps <- a %>%
  filter(indel == 0) %>%
  mutate(alleleA_sort = ifelse(alleleA < alleleB, alleleA, alleleB),
         alleleB_sort = ifelse(alleleA < alleleB, alleleB, alleleA)) %>%
  mutate(cptid = paste0(chromosome, ":", position, "_", alleleA_sort, "_", alleleB_sort))

# Create CPTIDs for indels (sort alleles by smallest first)
indels <- a %>%
  filter(indel == 1) %>%
  mutate(alleleA_sort = ifelse(nchar(alleleA) <= nchar(alleleB), alleleA, alleleB),
         alleleB_sort = ifelse(nchar(alleleA) > nchar(alleleB), alleleA, alleleB)) %>%
  mutate(cptid = paste0(chromosome, ":", position, "_", alleleA_sort, "_", alleleB_sort))

# Check for duplicated ids on SNPs
snps %>% 
  group_by(cptid) %>%
  filter(n() > 1) %>%
  data.frame()

# Check for duplicated ids on indels
dupes <- indels %>% 
  group_by(cptid) %>%
  filter(n() > 1) %>%
  data.frame()

# Append rsid to indels with duplicated cptids
if(nrow(dupes > 0)){
  dupevec <- which(indels$cptid %in% dupes$cptid)
  indels$cptid[dupevec] <- paste0(indels$cptid[dupevec], "___",indels$rsid[dupevec])
}

# Merge SNPs and indels back together
b <- rbind(snps, indels) %>%
  arrange(chromosome, position)


##Create a marker name column
snps.index <- which(nchar(assoc[,alleleA])==1 & nchar(assoc[,alleleB])==1)
indels.index <- which(nchar(assoc[,alleleA])>1 | nchar(assoc[,alleleB])>1)

assoc[snps.index,name:=paste(chromosome,position,alleleA,alleleB,sep=":")]
assoc[intersect(snps.index,which(alleleA>alleleB)),name:=paste(chromosome,position,alleleB,alleleA,sep=":")]
assoc[indels.index,name:=paste(chromosome,position,alleleA,alleleB,sep=":")]
assoc[intersect(indels.index,which(nchar(alleleB)>nchar(alleleA))),name:=paste(chromosome,position,alleleB,alleleA,sep=":")]
assoc[,name:=paste0("chr",name)]

