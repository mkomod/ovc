# Load data ---
TCGA.clinical <- data.table::fread("../data/TCGA_clinical.csv")
TCGA.radiomics <- data.table::fread("../data/TCGA_radiomics.csv")
TCGA.rna <- data.table::fread("../data/HT_HG-U133A")


# Construct dataset ---
# Tidy up mRNA data
TCGA.rna.IDs <- colnames(TCGA.rna)[-1]
TCGA.rna.sample.names <- unlist(TCGA.rna[ , 1])
TCGA.rna <- t(TCGA.rna[ ,-1])
TCGA.rna <- as.data.table(TCGA.rna)
colnames(TCGA.rna) <- gsub("[-|@]", "_", TCGA.rna.sample.names)
rownames(TCGA.rna) <- NULL
TCGA.rna$ID <- gsub("-", ".", TCGA.rna.IDs)

# Tidy up clinical data
colnames(TCGA.clinical) <- c(
	"ID",                                 "Age_at_diagnosis",
	"Stage",                              "Lymphatic_invasion",
	"Grade",                              "Tumor_residual_disease",
	"Surgery_outcome",                    "Venous_invasion",
	"Overall_survival_days",              "OS_event",
	"Primary_therapy_outcome_success",    "RPV",
	"eRPV",                               "Thickness",
	"Progression_free_survival_days",     "PFS_event",
	"Molecular_subtype",                  "Stromal_cell_content",
	"Tumour_cell_content"
)
# Dichotamise age 
TCGA.clinical$Age <- as.factor(TCGA.clinical$Age_at_diagnosis > 60)
TCGA.clinical$Stage <- as.factor(TCGA.clinical$Stage)

# Tidy up radiomics data
colnames(TCGA.radiomics)[1] <- "ID"

# Merge the datasets
TCGA <- merge(TCGA.clinical, TCGA.radiomics, by="ID")
TCGA <- merge(TCGA, TCGA.rna, by="ID")

# Remove outlier
TCGA.outlier <- TCGA[67, ]
TCGA <- TCGA[-67, ]
# save(TCGA, file="../RData/tcga_full.RData")

# IDs not present in radiomics dataset
# TCGA.radiomics$ID[!(TCGA.radiomics$ID %in% TCGA.rna$ID)]

# Cleanup
rm(list=c("TCGA.rna.IDs", "TCGA.rna.sample.names"))
