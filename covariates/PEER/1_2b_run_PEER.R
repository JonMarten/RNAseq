# ----------------------------------------
# Running peer analysis 
# ----------------------------------------
# Gene expression data has 4143samples x 17964 genes
# Gene data TMM nornalised and each gene is inverse rank normalised

# Covariates data includes 
# we including age, sex, BMI and blood cell traits: 
#BASO_PCT___RNA, EO_PCT___RNA, LYMPH_PCT___RNA, MONO_PCT___RNA,NEUT_PCT___RNA
#WBC_10_9_L___RNA, IRF_PCT___RNA, RET_PCT___RNA, HCT_PCT___RNA, HGB_g_dL___RNA, MCH_pg___RNA
#MCHC_g_dL___RNA, MCV_fL___RNA, RBC_10_12_L___RNA, RDW_SD_fL___RNA, MPV_fL___RNA, PCT_PCT___RNA, PDW_fL___RNA, PLT_F_10_9_L___RNA

library(peer)
setwd("/home/jm2294/rds/rds-jmmh2-projects/interval_rna_seq/analysis/04_phase2_full_analysis/peer_factors")
gene_expr <- read.csv("peer_InputFiles/PEER_Proceesed/GeneExpr_PEER_TmmInvRankNormalised_swapsSwapped_mismatchRemoved.csv", row.names = 1)
covariates <-read.csv("peer_InputFiles/PEER_Proceesed/Covariates_for_PEER.csv", row.names = 1)

table(rownames(gene_expr) == rownames(covariates))

gene_expr <- as.matrix(gene_expr)
covariates <- as.matrix(covariates)

# Run PEER

# Set up model
model <- PEER()

# adds intercept term to PEER linear model fits

PEER_setPhenoMean(model, gene_expr)
dim(PEER_getPhenoMean(model))

# adds additional factor to account for the mean expression
PEER_setAdd_mean(model, TRUE)

# PEER will find factors independent of these covariates
PEER_setCovariates(model, covariates)
dim(PEER_getCovariates(model))

# For automatic factor selection the PEER Nature Protocols paper recommends 
# setting the number of factors to 25% of the sample size to a maximum of 100.
# We can then select the appropriate number of factors by looking at their variance.
nK <- min(30, round(nrow(gene_expr)/4))
nK
PEER_setNk(model, nK)
PEER_getNk(model)

# Increase to 10,000 if the next command gives a message about failing to converge
PEER_setNmax_iterations(model, 1000)  

# Run PEER factor analysis
PEER_update(model)

# Get the hidden factors that PEER identified
factors <- PEER_getX(model)
dim(factors)
rownames(factors) <- rownames(gene_expr)
colnames(factors) <- c(colnames(covariates), "intercept", paste0("PEER", 1:nK))
saveRDS(factors, "peer_30Fact_100Iter_sampleSwapsFixed/PEER_factors.rds")
write.table(factors,"peer_30Fact_100Iter_sampleSwapsFixed/PEER_factors.txt",quote=F,sep="\t")

# Get the variance explained by each covariate and factor - we will use this to choose
# the number of factors to use in downstream analyses
factor_variance <- PEER_getAlpha(model)
factor_variance <- 1/factor_variance # Convert from inverse variance to variance
rownames(factor_variance) <- c(colnames(covariates), "intercept", paste0("PEER", 1:nK))
factor_variance <- factor_variance[-(1:(ncol(covariates)+1)),] # ignore covariates and intercept 
saveRDS(factor_variance, "peer_30Fact_100Iter_sampleSwapsFixed/PEER_factors_varaince_explained.rds")
write.table(factor_variance,"peer_30Fact_100Iter_sampleSwapsFixed/PEER_factorvariance.txt",quote=F,sep="\t")

# Look at this plot: for downstream analysis you want to choose the first N factors 
# before the variance drops to 0
pdf("peer_30Fact_100Iter_sampleSwapsFixed/PEER_factor_kN_selection.pdf")
plot(factor_variance, type="l", ylab="Variance of PEER factor weights", xlab="", xaxt="n")
points(factor_variance, pch=19)
dev.off()


# Other information from the PEER analysis we don't need, but you might
# want to save anyway:
# The weights correspond to how much each factor/covariate affects each gene
weights <- PEER_getW(model)
colnames(weights) <- c(colnames(covariates), "intercept", paste0("PEER", 1:nK))
rownames(weights) <- colnames(gene_expr)
write.table(weights,"peer_30Fact_100Iter_sampleSwapsFixed/PEER_weights.txt",quote=F,sep="\t")
