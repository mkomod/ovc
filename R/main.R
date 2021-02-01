library(data.table)
library(survival)
library(SmCCNet)                       # Network CCA (Shi et al.)
library(rcca)                          # Rel CCA (mkomod/rcca)

RUN_ALL <- FALSE                       # Re-run the analysis, might take some time
CORES <- 32                            # Number of cluster cores to use
PERMUTATIONS <- 1000                   # Number of permutations


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
# Note: the following should be ran on a cluster. Otherwise it may take some time
if (RUN_ALL) {
    source("02-cluster.R")
} else {
    load("../RData/rcca_permutation_validation.RData")
}
rcca.pv <- grid.as.matrix(rcca.pv)
rcca.pv$x <- rcca.pv$x[-(1:3)]
rcca.pv$z <- rcca.pv$z[-(1:3), ]

# Canonical vectors and Projections
n_cca_vectors <- 3
TCGA.cca <- rcca::rCCA(TCGA.rna.std, TCGA.rad.std, 1.1, 1.6, K=n_cca_vectors)
TCGA.rna.proj <- TCGA.rna.std %*% TCGA.cca$w1
TCGA.rad.proj <- TCGA.rad.std %*% TCGA.cca$w2

# Datasets with projections
colnames(TCGA.rad.proj) <- paste0("rRad", 1:n_cca_vectors)
colnames(TCGA.rna.proj) <- paste0("rRNA", 1:n_cca_vectors)
TCGA.proj <- cbind(TCGA[, c("Age", "Overall_survival_days", "OS_event", "Stage",
			    "RPV")], TCGA.rad.proj,TCGA.rna.proj)


# Models ---
source("03-models.R")


# Networks ---
if (RUN_ALL) {
    rcca_and_simMat <- constructSimilarityMatrix(TCGA.rna.std, TCGA.rad.std,
						 1.1, 1.6)
    TCGA.cca_robust <- rcca_and_simMat[2:3]
    A_bar <- rcca_and_simMat$Abar
} else {
    load("../RData/networks.RData")
}
mods <- SmCCNet::getMultiOmicsModules(A_bar, ncol(TCGA.rna.std), PlotTree=F)


# Figures and Tables ---
source("04-figures.R")
source("05-tables.R")
