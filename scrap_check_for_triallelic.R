# scrap code to check if all chunks have one or more triallelic SNPs

snpstats <- fread("interval_b38_snpstats_filtered.txt", skip = 8, data.table = F)
dupeAlt <- snpstats %>% group_by(alternate_ids) %>% filter(n() > 1) %>% data.frame
dupeRSIDs <- snpstats %>% group_by(rsid) %>% filter(n() > 1) %>% data.frame

chunks <- fread("~/projects/RNAseq/test_run_chunks/chunklist.txt", data.table= F)
names(chunks) <- c("chr","start","end")

chunks22 <- chunks %>% filter(chr==22)

tria_in_chunk <- data.frame("chunk" = 529:540, "tria" = NA)
for(i in 1:nrow(chunks22)){
  vec <- which(dupeRSIDs$position > chunks22$start[i] & dupeRSIDs$position < chunks22$end[i])
  tria_in_chunk$tria[i] <- paste(c(min(vec), max(vec)), collapse = ":")
}