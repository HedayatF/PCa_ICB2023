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
# Extract Expression Matrix and Clinical Data
ex <- assay(data)
write.table(ex, file = "results/ExpressionMatrix(counts).txt", quote = F, sep = "\t" ) #will use this in local system to create barcode.txt
cli <- as.data.frame(colData(data))
cli <- cli[, !(names(cli) %in% c("treatments", "primary_site", "disease_type"))]
write.table(cli, file = "results/clinicalData.txt", quote = F, sep = "\t")

# Retrieve Gene Information
gene <- as.data.frame(SummarizedExperiment::rowRanges(data))
gene <- gene[, c("gene_id", "gene_name", "gene_type")]
write.table(gene, file = "results/Geneinfo.txt", quote = F, sep = "\t", row.names = F)

# Filter for Protein-Coding Genes
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                   mirror =  "asia")
gene_metadata <- getBM(attributes=c("external_gene_name",
                                    "ensembl_gene_id_version",
                                    "gene_biotype"),
              filters = c("ensembl_gene_id_version"),
              values = rownames(ex),
              mart = mart)
protein_coding_genes <- gene_metadata[gene_metadata$gene_biotype == "protein_coding",]
ex <- ex[rownames(ex) %in% protein_coding_genes$ensembl_gene_id_version,]

#Reorganize Samples Based on Barcodes 
# The barcode file:
# - The first column contains unique sample barcodes (IDs) from TCGA-PRAD.
# - The second column indicates sample type (e.g., Normal or Cancer).
bar <- read.delim("barcode/barcode.txt")
rownames(bar) <- bar$barcode
ex <- ex[, bar$barcode]

# Create Design Matrix
gr <- factor(c(rep("Normal", 52), rep("Cancer", 502)))
design <- model.matrix(~0 + gr)
colnames(design) <- levels(gr)

#--------------------- Normalize Data ---------------------
dge <- DGEList(counts = ex)
keep <- filterByExpr(dge, design = design, min.count = 10, min.pct = 0.7)
dge <- dge[keep, , keep.lib.size = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_data <- voom(dge, design = design, plot = TRUE)
exn <- voom_data$E
boxplot(exn[, 1:30]) #check if data is normalized or not 
exn <- log2(2^exn + 1)
write.table(exn, "results/ExpressionMatrix(normalized).txt", quote = F, sep = "\t")

#--------------------- Differential Expression Analysis ---------------------
fit <- lmFit(voom_data, design = design)
contrast <- makeContrasts(Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrasts = contrast)
fit2 <- eBayes(fit2)

# Save all significant genes and DEGs
tT <- topTable(fit2, sort.by = "B", number = Inf)
degs <- subset(tT, abs(logFC) > 2 & adj.P.Val < 0.01)
write.table(degs, "results/degs.txt", quote = F, sep = "\t")
degsUP <- subset(degs, logFC > 2)
degsDOWN <- subset(degs, logFC < -2)

# Annotate Differential Expression with Gene Names 
gene <- read.delim("results/Geneinfo.txt") #if not already present
rownames(gene) <- gene$gene_id
degs$symbol <- gene[as.character(rownames(degs)), "gene_name"]
degsDOWN$symbol <- gene[as.character(rownames(degsDOWN)), "gene_name"]
degsUP$symbol <- gene[as.character(rownames(degsUP)), "gene_name"]
write.table(degsDOWN, file = "degsDOWN.txt", quote = F, sep = "\t" )
write.table(degsUP, file = "degsUP.txt", quote = F, sep = "\t" )
#now you can use these gene symbols for Enrichment Analysis

#--------------------- Visualization ---------------------
# Volcano Plot
volc <- subset(tT)
volc$Significant <- "No"
volc$Significant[volc$logFC > 2 & volc$adj.P.Val < 0.01] <- "Up"
volc$Significant[volc$logFC < -2 & volc$adj.P.Val < 0.01] <- "Down"
png("plots/volcano plot.png",height = 1600, width = 2000,res=300,units = "px" )
ggplot(volc, aes(logFC, -log10(adj.P.Val),color =Significant))+ 
  geom_point(size = 1.6,shape = 19) + theme_bw()+ 
  geom_vline(xintercept=c(-2, 2), col="#EF00FF",linetype = "longdash") + 
  geom_hline(yintercept=-log10(0.01), col="#EF00FF",linetype = "longdash") +
  scale_color_manual(values=c("blue", "gray", "red")) +
  ggtitle("Volcano Plot")
dev.off()

#  Heatmap of DEGs 
exn <- read.delim("results/ExpressionMatrix(normalized).txt")
degs <- read.delim("results/degs.txt")
exn <- data.frame(exn)
degs$ID <- rownames(degs)
heatmap_data <- exn[rownames(exn) %in% degs$ID,]
rownames(heatmap_data) <- gene[as.character(rownames(degs)), "gene_name"]
degs <- degs[,-8]
rownames(degs) <- degs$symbol

#making annotation_col:
coloo <- colnames(heatmap_data) 
annotation1 <- data.frame(gr, coloo)
rownames(annotation1) <- annotation1$coloo
colnames(annotation1) = "Group"
annotation1 <- annotation1[,-2,drop=F]

annotation2 <- data.frame(degs$logFC)
rownames(annotation2) <- degs$symbol
colnames(annotation2) = "Up/Down regulated"
annotation2$`Up/Down regulated`[degs$logFC > 2 & degs$adj.P.Val < 0.01] <- "Up"
annotation2$`Up/Down regulated`[degs$logFC < -2 & degs$adj.P.Val < 0.01] <- "Down"

ann_colors = list(
  "Group" = c( Normal= "#00CC00",  Cancer= "#FF0000"),
  "Up/Down regulated" = c(Down = "#CC99FF", Up = "#ffff00"))

png("plots/heatmap_default.png",height = 3500, width = 3300,res=300,units = "px" )
pheatmap(heatmap_data,color= bluered(256), show_rownames=F ,border_color =NA
         ,annotation_col = annotation1,annotation_row = annotation2, 
         scale = "row", labels_col = gr, annotation_colors = ann_colors,
         annotation_names_row = F,annotation_names_col = F,fontsize_row = 4.5,
         fontsize_col=1,angle_col=45,
         clustering_distance_rows = "correlation",
         clustering_distance_cols= "correlation",
         clustering_method = "complete")
dev.off()

#--------------------- get expression matrix of DEGs for the next step ---------------------
exn$symbol <- gene[as.character(rownames(exn)), "gene_name"]
exn2 <- exn[exn$symbol %in% degs$symbol, ] 
rownames(exn2) <- exn2$symbol
exn2 <- exn2[, -555] #emove symbol column
write.table(exn2, "Cox_regression/expression_matrix_of_DEGs.txt", sep = "\t", quote = F)

# END OF THE SCRIPT
