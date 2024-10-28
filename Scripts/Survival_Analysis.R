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
write.table(exs, file = "D:\Project\prostate\Cox_regression/Z_score.txt",
            quote = FALSE, sep = "\t")
