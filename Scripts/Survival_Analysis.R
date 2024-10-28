################################################################################
# Prostate Cancer Survival Analysis with Univariate Cox Regression
# TCGA-PRAD Dataset
# 
# Description: This script performs survival analysis on prostate cancer DEGs 
# using Cox regression to identify genes linked to poor prognosis.
#
# Note: This code requires two preprocessed inputs:
# 1. A DEGs expression matrix from "Preprocessing_and_DEA.R".
# 2. Processed clinical data with 3 columns:
#       - `barcode: only tumor tissue barcodes 
#       - `os.days`: Overall survival days, created by merging "days_to_last_follow_up" (for alive patients) and 
#                    "days_to_death" (for dead patients).
#       - `os.event`: Event status, derived from "vital_status", with 0 indicating alive and 1 indicating dead
################################################################################

#--------------------- Load Libraries ---------------------
library(dbplyr)
library(survival)
library(survminer)
library(survMisc)

setwd("D:/Project/prostate")

#--------------------- Read DEG Expression Matrix ---------------------
ex <- read.delim("Cox_regression/expression_matrix_of_DEGs.txt")

#--------------------- Cox Regression Analysis ---------------------
# Remove normal sample columns
ex <- ex[, 53:ncol(ex)]

# Transpose matrix so genes are in columns
ex <- data.frame(t(ex))

# Calculate Z-scores for each gene
exs <- data.frame(scale(ex))

# List of DEGs for Cox regression
covariates <- colnames(exs)

#--------------------- Load and Process Clinical Data ---------------------
# Read clinical data and set rownames to barcode
cli <- read.delim("Cox_regression/Processed_clinical.txt")
rownames(cli) <- cli$barcode

# Merge clinical survival data (os.days, os.event) with expression data (exs)
exs$os.days <- cli[rownames(exs), "os.days"]
exs$os.event <- cli[rownames(exs), "os.event"]

# Remove samples with missing clinical data
exs <- na.omit(exs)
write.table(exs, file = "Cox_regression/Z_score.txt",
            quote = FALSE, sep = "\t")

#--------------------- Calculate Hazard Ratios for All DEGs ---------------------
# Create Cox regression formulas for each gene
univ_formulas <- sapply(covariates, function(x) {
  as.formula(paste('Surv(time = os.days, event = os.event) ~', x))
})


# Fit Cox models for each gene
univ_models <- lapply(univ_formulas, function(x) {
  coxph(x, data = exs)
})

# Extract results from each model
univ_results <- lapply(univ_models, function(x) {
  x <- summary(x)
  p.value <- signif(x$wald["pvalue"], 4)
  wald.test <- signif(x$wald["test"], 4)
  beta <- signif(x$coef[1], 4)              # Coefficient beta
  HR <- signif(x$coef[2], 4)                # Hazard Ratio
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  c(beta = beta, `HR (95% CI for HR)` = HR, wald.test = wald.test, p.value = p.value)
            names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                       })
# Convert results to data frame and save
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
res$symbol <- rownames(res)
write.table(res, "Cox_regression/cox_regression.txt", quote = F, sep = "\t")

res <- subset(res,`HR (95% CI for HR)`>1 & p.value < 0.05)
write.table(res, "Cox_regression/Filtered_cox_regression.txt", quote = F, sep = "\t")

# After cox regression > bad prognostic genes > PPI > CytoHubba > HUBs 

################################################################################
#                                    END OF THE SCRIPT                         #
################################################################################
