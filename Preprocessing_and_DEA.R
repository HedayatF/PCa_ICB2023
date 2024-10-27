################################################################################
# Prostate Cancer vs Normal Gene Expression Analysis
# TCGA-PRAD Dataset
# 
# Description: This script performs RNA-seq analysis for differential gene expression analysis in Prostate Cancer and Normal samples using TCGA-PRAD data.
# The steps include:
#    1) Data Query and Download
#    2) Data Preparation
#    3) Filtering and Normalization
#    4) Differential Expression Analysis
#    5) DEG Visualization:
#    6) Save Processed Data
#
################################################################################

#--------------------- Load Libraries ---------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(gplots)
library(biomaRt)
library(dbplyr)
library(pheatmap)
library(survival)
library(survminer)
library(survMisc)

#--------------------- Set Working Directory ---------------------
setwd("D:/Project/prostate/")

#--------------------- Query and Download Data ---------------------
query <- GDCquery(project = "TCGA-PRAD", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  experimental.strategy = "RNA-Seq")

GDCdownload(query = query, method = "api", directory = "data/", files.per.chunk = 5)
data <- GDCprepare(query = query, directory = "data/", summarizedExperiment = TRUE)

#--------------------- Extract and Prepare Data ---------------------
# Get expression matrix and clinical data
ex <- assay(data)
write.table(ex, file = "results/expression matrix (counts).txt", quote = F, sep = "\t" )
cli <- as.data.frame(colData(data))[ , -c(25, 58, 60)] #removed treatments, primary_site, disease_type, and papers
write.table(cli, file = "results/cli.txt", quote = F, sep = "\t")

# Get gene names and filter for protein-coding genes
gene <- as.data.frame(SummarizedExperiment::rowRanges(data))
gene <- gene[, c("gene_id", "gene_name", "gene_type")]
write.table(gene, file = "results/info gene.txt", quote = F, sep = "\t", row.names = F)
