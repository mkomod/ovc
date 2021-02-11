library(data.table)
library(survival)
library(SmCCNet)                       # Network CCA (Shi et al.)
library(rcca)                          # Rel CCA (mkomod/rcca)

CORES <- 8
PERMUTATIONS <- 1000
RUN_ALL <- FALSE

source("00-functions.R")
source("01-process_data.R")

# Variable selection ---
# Radiomics and RNA features i.e. we exc TCGA.clinical covariates
if (RUN_ALL) {
    covariates <- names(TCGA)[(ncol(TCGA.clinical) + 1):(ncol(TCGA))]
    covariate.information <- matrix(nrow=length(covariates), ncol=5)
    colnames(covariate.information) <- c("coef", "exp_coef", "se_coef", "z", 
					 "p_value")
    rownames(covariate.information) <- covariates
    for (i in seq_along(covariates)) {
	covariate <- covariates[i]
	model.formula <- as.formula(paste(
					  "Surv(Overall_survival_days, OS_event)",
					  covariate,
					  sep=" ~ ")) 
	cfit <- survival::coxph(model.formula, data=TCGA)
	covariate.information[i, ] <- as.vector(summary(cfit)$coeff)
    }
    save(covariate.information, file="../RData/covariate_information.RData")
} else {
    load("../RData/covariate_information.RData")
}
covariate.rad.info <- sort(covariate.information[
	rownames(covariate.information) %in% names(TCGA.radiomics), "p_value"])
covariate.rna.info <- sort(covariate.information[
	rownames(covariate.information) %in% names(TCGA.rna), "p_value"])

# Set a threshold and select variables
pval <- 0.05
covariates.selected <- which(covariate.information[ , "p_value"] <= pval)

# Datasets from selected variables ---
TCGA.rad.covariates <- names(covariate.rad.info[which(covariate.rad.info<=pval)])
TCGA.rna.covariates <- names(covariate.rna.info[which(covariate.rna.info<=pval)])
TCGA.rad.subset <- TCGA[ , ..TCGA.rad.covariates]
TCGA.rna.subset <- TCGA[ , ..TCGA.rna.covariates]

# Center and standardise the data
TCGA.rna.mat <- as.matrix(TCGA.rna.subset)
TCGA.rad.mat <- as.matrix(TCGA.rad.subset)
TCGA.rna.std <- scale(TCGA.rna.mat, T, T)
TCGA.rad.std <- scale(TCGA.rad.mat, T, T)
# save(list=c("TCGA.rna.std", "TCGA.rad.std"), file="../RData/tcga_data.RData")


# Correlations ---
TCGA.rna.cor <- cor(TCGA.rna.mat)
TCGA.rna.cov <- cov(TCGA.rna.mat)
TCGA.rad.cor <- cor(TCGA.rad.mat)
TCGA.rad.cov <- cov(TCGA.rad.mat)
TCGA.rna.rad.cor <- cor(TCGA.rna.mat, TCGA.rad.mat)
TCGA.rna.rad.cov <- cov(TCGA.rna.mat, TCGA.rad.mat)
colnames(TCGA.rna.cor) <- rownames(TCGA.rna.cor) <- NULL
colnames(TCGA.rad.cor) <- rownames(TCGA.rad.cor) <- NULL

# Principle Components (PCs) ---
S.rad <- t(TCGA.rad.std) %*% TCGA.rad.std
S.rna <- t(TCGA.rna.std) %*% TCGA.rna.std
TCGA.rad.eigen <- eigen(S.rad)
TCGA.rna.eigen <- eigen(S.rna)

# Variance explained
(cumsum(TCGA.rad.eigen$values) / sum(TCGA.rad.eigen$values))[1:20]
(cumsum(TCGA.rna.eigen$values) / sum(TCGA.rna.eigen$values))[1:70]

# Projections onto PC
TCGA.rad.pc.proj <- TCGA.rad.std %*% TCGA.rad.eigen$vectors
TCGA.rna.pc.proj <- TCGA.rna.std %*% TCGA.rna.eigen$vectors

# Create Dataset with PC projections
colnames(TCGA.rad.pc.proj) <- paste0("pRad", 1:ncol(TCGA.rad.pc.proj))
colnames(TCGA.rna.pc.proj) <- paste0("pRNA", 1:ncol(TCGA.rna.pc.proj))
TCGA.pc.proj <- cbind(TCGA[,c("Age","Overall_survival_days","OS_event","Stage")],
	    TCGA.rad.pc.proj,TCGA.rna.pc.proj)


# Canonical Correlation Analysis (CCA) ---
# Tune hyperparameter
# Note: the following should be ran on a cluster, the experimental branch of
# the rcca package needs to be built and installed. If this is not available
# see: funcs.R for a permutation based validaiton or cross validation scheme
if (RUN_ALL) {
    source("02-cluster.R")
} else {
    load("../RData/rcca_permutation_validation.RData")
}
rcca.pv <- grid.as.matrix(rcca.pv)
rcca.pv$x <- rcca.pv$x[-(1:3)]
rcca.pv$z <- rcca.pv$z[-(1:3), ]
l1l2.min <- expand.grid(l1=rcca.pv$x, l2=rcca.pv$y)[which.min(rcca.pv$z), ]
l1 <- l1l2.min$l1
l2 <- l1l2.min$l2

# Canonical vectors and Projections
n_cca_vectors <- 3
TCGA.cca <- rcca::rCCA(TCGA.rna.std, TCGA.rad.std, l1, l2, K=n_cca_vectors)
TCGA.rna.proj <- TCGA.rna.std %*% TCGA.cca$w1
TCGA.rad.proj <- TCGA.rad.std %*% TCGA.cca$w2

# TCGA.cca.rad.names <- names(TCGA.rad.subset)[which(TCGA.cca$w2[ , 2] != 0)][2]
# TCGA.rad.std <- TCGA.rad.std[ , colnames(TCGA.rad.std) != TCGA.cca.rad.names]

# Vector sensitivity analysis
# l1l2 <- expand.grid(seq(0.5, 1.5, by=0.1), seq(1, 2, by=0.1))
# TCGA.cca.sen <- ccaSensitivity(TCGA.rna.std, TCGA.rad.std, l1l2)
# save(TCGA.cca.sen, file="../RData/rcca_sensitivity.RData")

# Datasets with projections
colnames(TCGA.rad.proj) <- paste0("rRad", 1:n_cca_vectors)
colnames(TCGA.rna.proj) <- paste0("rRNA", 1:n_cca_vectors)
TCGA.proj <- cbind(TCGA[, c("Overall_survival_days", "OS_event", 
			    "Progression_free_survival_days", "PFS_event",
			    "Stage", "Age", 
			    "RPV")], TCGA.rad.proj,TCGA.rna.proj)

# Networks ---
if (RUN_ALL) {
    rcca_and_simMat <- constructSimilarityMatrix(TCGA.rna.std, TCGA.rad.std,
						 l1, l2)
    TCGA.cca_robust <- rcca_and_simMat[2:3]
    A_bar <- rcca_and_simMat$Abar
    save(list=c("A_bar", "TCGA.cca_robust"), file="../RData/networks.RData")
} else {
    load("../RData/networks.RData")
}
mods <- SmCCNet::getMultiOmicsModules(A_bar, ncol(TCGA.rna.std), PlotTree=F)
mod.names <- colnames(cbind(TCGA.rna.std, TCGA.rad.std))

# Models, Figrures and Tables ---
source("03-models.R")
source("04-figures.R")
source("05-tables.R")

