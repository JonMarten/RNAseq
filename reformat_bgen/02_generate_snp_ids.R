# Make new alternate ID for bgen 1.2 files
setwd("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq")
library(dplyr)
library(data.table)
library(stringr)

chr <- 22

# Read in SNP stats output from qctool and add column to check if variant is an indel
snpstats <- fread(paste0("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/snp_stats/impute_",chr,"_interval_snp_stats_unfiltered.txt"), skip = 8, data.table=F)
snpstats <- snpstats %>%
  mutate(indel = ifelse(nchar(alleleA)>1 | nchar(alleleB)>1, 1, 0)) %>%
  #filter(rsid != ".") %>%
  mutate(alleleA_sort = ifelse(alleleA <= alleleB, alleleA, alleleB),
         alleleB_sort = ifelse(alleleA <= alleleB, alleleB, alleleA)) %>%
  mutate(match_id = paste0(str_pad(chr,2,pad = "0"),"_", position, "_", alleleA, "_", alleleB, ":", rsid))


# read in b38 map
b38 <- fread(paste0("/home/jm2294/GENETIC_DATA/INTERVAL/RNAseq/b37_b38_liftover/INTERVAL_imputed_liftover_hg38_cleaned/INTERVAL_imputed_liftover_hg38_cleaned_chr",chr,".vcf"), data.table = F, skip = 5)
names(b38) <- c("CHROM.b38","POS.b38","match_id","REF.b38","ALT.b38")

# b38 <- b38 %>%
#   distinct %>%
#   filter(rsid != ".") %>%
#   mutate(REF.b38_sort = ifelse(REF.b38 < ALT.b38, REF.b38, ALT.b38),
#          ALT.b38_sort = ifelse(REF.b38 < ALT.b38, ALT.b38, REF.b38))  %>%
#   mutate(match_id = paste0(rsid,"_",REF.b38_sort,"_",ALT.b38_sort))

snpsmerge <- inner_join(snpstats, b38, by = "match_id")

# Make list of SNPs that don't have b38 positions
no_b38 <- snpstats$match_id[which(!snpstats$match_id %in% snpsmerge$match_id)]
no_snpstats <- b38$match_id[which(!b38$match_id %in% snpsmerge$match_id)]

cat("\n",length(no_b38), "variants do not match to b38.")

drops <- snpstats %>%
  filter(match_id %in% no_b38) %>%
  mutate(chrpos = paste0(chromosome,":",position)) %>%
  pull(chrpos)

# Create CPTIDs for SNPs (sort alleles in alphabetical order)
snps <- snpsmerge %>%
  filter(indel == 0) %>%
  mutate(cptid.b38 = paste0(CHROM.b38, ":", POS.b38, "_", REF.b38_sort, "_", ALT.b38_sort))

# Create CPTIDs for indels (sort alleles by smallest first)
indels <- snpsmerge %>%
  filter(indel == 1) %>%
  filter(!is.na(CHROM.b38)) %>%
  mutate(REF.b38_sort = ifelse(nchar(REF.b38) <= nchar(ALT.b38), REF.b38, ALT.b38),
         ALT.b38_sort = ifelse(nchar(REF.b38) > nchar(ALT.b38), REF.b38, ALT.b38)) %>%
  mutate(cptid.b38 = paste0(CHROM.b38, ":", POS.b38, "_", REF.b38_sort, "_", ALT.b38_sort))

# Merge SNPs and indels back together
b <- rbind(snps, indels) %>%
  arrange(chromosome, position)

dupes <- b %>% 
  group_by(cptid.b38) %>%
  filter(n() > 1) %>%
  data.frame()
if(nrow(dupes > 0)){
  dupevec <- which(b$cptid.b38 %in% dupes$cptid.b38)
  b$cptid.b38[dupevec] <- paste0(b$cptid.b38[dupevec], "___",b$rsid.x[dupevec])
}

# For qctool:
# -map-id-data	
# Update the chromosome, position, IDs and/or alleles of a set of SNPs with new values. 
# The argument must be a file with six named columns giving the original SNPID, rsid, chromosome, position and alleles, followed by another six columns containing the values to replace with. SNPs not in this file will be passed to the output file unchanged. This option only affects the identifying data, not genotypes themselves.
#map.txt is a file with 12 (named) columns: old SNPID, rsid, chromosome, position, alleleA, alleleB, and new SNPID, rsid chromosome, position, alleleA, alleleB. (from https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=OXSTATGEN;841f4240.1607). The "map" file given to -map-id-data must be a text file with twelve named columns, in the following order: the current SNPID, rsid, chromosome, position, first and second alleles, followed by the desired updated SNPID, rsid, chromosome, position and alleles. The first line is treated as column names (currently it doesn't matter what these are called.) Variants not in this file are not affected by the mapping, and will be output unchanged. (https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/examples/altering_id_data.html)
newMap <- b %>%
  mutate(rsid=rsid.x, alternate_id.b38=cptid.b38,rsid.38=rsid.x, chromosome.38=CHROM.b38, position.b38=POS.b38,alleleA.b38=alleleA, alleleB.b38=alleleB) %>%
  select(alternate_ids, rsid, chromosome, position, alleleA, alleleB,
         alternate_id.b38, rsid.38, chromosome.38, position.b38, alleleA.b38, alleleB.b38)
write.table(newMap, quote = F, sep = " ", col.names = T, row.names = F, file = "b37_b38_liftover/INTERVAL_chr22_b37_to_b38_map.txt")

# As not all SNPs are currently mappable in the b38 data due to ambiguous rsids (eg snps coded as "."), SNPs have to be filtered to retain only those which have unambiguous b38 positions.
#--incl-snps	Exclude all SNPs not in the given file(s) from the analysis. The format of this file is the same as that output by the -write-snp-excl-list option. It must have six columns interpreted as SNPID, rsid, chromosome, position, first and second alleles. 
snplist <- newMap %>%
  select(alternate_id.b38)
write.table(snplist, row.names = F, col.names = F, quote = F, file = "b37_b38_liftover/c22_b38_filter_snps.txt")


