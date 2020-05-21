library(data.table)
library(dplyr)
a <- fread("rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.rel", data.table = F)
b <- fread("rnaseq_phase1.typed.maf1pc.geno3pc.hwe1minus5.complex-region-excl.indep-50-5-1.5.pruning.genome.pca.rel.id", data.table = F)

a1 <- as.matrix(a)
colnames(a1) <- as.character(b$IID)
rownames(a1) <- as.character(b$IID)

a1[upper.tri(a1, diag = FALSE)] %>% sort %>% max

fwrite(a1, file = "kinship_matrix.csv", row.names = T)

# drop one of each realated pair
dropped <- character()
a2 <- a1
maxrel = 1
while(maxrel > 0.1) {
  maxrel <- a2[upper.tri(a2, diag = FALSE)] %>% sort %>% max
  cat("\nMaxrel: ", maxrel)
  rels <- which(arr.ind = T, a2 == maxrel)
  removerow <- rels[1,1]
  removeid <- rownames(rels)[1]
  cat("\nRemoving: ", removeid)
  a2 <- a2[-removerow, -removerow]
  dropped <- c(dropped, removeid)
  rm(rels, removerow, removeid)
}
# for some reason this code drops one ID with relatedness < 0.1