library(Seurat) #load Seurat 3.1
library(tidyverse)
library(pheatmap)
library(viridis)
library(reshape2)

######################################
### TRBV11-2 CD4 T cell phenotype ####
######################################
# Generates Figures 5D,5E
# This file need the Seurat object: misc_SeuratObj_submission.rds downloaded in input folder
### within CD4 clustering #################################################################
merge <- readRDS("input/misc_SeuratObj_submission.rds")
merge.CD4 <- subset(merge, subset = coarsecelltype %in% c("CD4"))
NormalizeData(merge.CD4, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
merge.CD4 <- FindVariableFeatures(merge.CD4)
merge.CD4 <- ScaleData(merge.CD4, features = VariableFeatures(merge.CD4))
merge.CD4 <- RunPCA(merge.CD4, features = VariableFeatures(merge.CD4))
merge.CD4 <- FindNeighbors(merge.CD4, reduction = "pca", dims = 1:15)
merge.CD4 <- FindClusters(merge.CD4, resolution = 0.6, verbose = FALSE)
merge.CD4 <- RunUMAP(merge.CD4, reduction = "pca", dims = 1:15)

merge.CD4.UnSort <- subset(merge.CD4, subset = Sorted == "UnSort")
meta <- merge.CD4.UnSort@meta.data

meta$Class <- factor(meta$Class, levels = c("HC","COVID","MIS-C"))
DimPlot(merge.CD4.UnSort, group.by = "seurat_clusters", label = TRUE)

# MIS-C TRBV11-2 CD4 T cell differential marker test
merge.CD4.UnSort$TRBV[is.na(merge.CD4.UnSort$TRBV)] <- "ND"
merge.CD4.UnSort$TRBV_Class = paste(merge.CD4.UnSort$TRBV, merge.CD4.UnSort$Class, sep = "_")
merge.CD4.UnSort$TRBV112 <- merge.CD4.UnSort$TRBV == "TRBV11-2"

merge.CD4.UnSort.misc <- subset(merge.CD4.UnSort, subset = Class == "MIS-C")
Idents(merge.CD4.UnSort.misc) = "TRBV"
TRBVmarker.cite = FindMarkers(merge.CD4.UnSort.misc, ident.1 = "TRBV11-2", assay = "CITE")
TRBVmarker.rna = FindMarkers(merge.CD4.UnSort.misc, ident.1 = "TRBV11-2", assay = "RNA", logfc.threshold = 0.2, min.pct = 0.05, only.pos = TRUE)

# TRBVmarker.rna.genes <- rownames(TRBVmarker.rna[TRBVmarker.rna$p_val<0.2,])


##########################################################
### Plot average heatmap of marker genes -- Figure 5D ####
##########################################################
# MIS-C TRBV11-2 vs non-TRBV11-2
TRBVmarker.rna <- TRBVmarker.rna %>% rownames_to_column("gene")

merge.CD4.UnSort.misc$TRBV[is.na(merge.CD4.UnSort.misc$TRBV)] <- "ND"
merge.CD4.UnSort.misc$TRBV112 <- merge.CD4.UnSort.misc$TRBV == "TRBV11-2"

# plot all 4 groups, but use the DE genes from TRBV11-2_MISC vs nonTRBV11-2_MISC
Idents(merge.CD4.UnSort) <- "TRBV112_Class"
avg.CD4.UnSort <- data.frame(log1p(AverageExpression(merge.CD4.UnSort, verbose = TRUE)$RNA))
avg.CD4.UnSort$gene <- rownames(avg.CD4.UnSort)

avg.CD4.UnSort.marker <- filter(avg.CD4.UnSort, gene %in% TRBVmarker.rna[TRBVmarker.rna$p_val<0.2,]$gene)
avg.CD4.UnSort.marker <- avg.CD4.UnSort.marker[,c("FALSE_HC","FALSE_COVID","FALSE_MIS.C","TRUE_HC","TRUE_COVID","TRUE_MIS.C")]
pheatmap(avg.CD4.UnSort.marker[,1:6], cluster_cols = FALSE,
         border_color = "white", scale = "row",
         cellwidth = 25,cellheight = 8,
         # filename = file.path("output/CD4.MISC.TRBV112.marker.gene.heatmap.pdf"),
         fontsize_row = 6)
dev.off()


##########################################################
### Plot average heatmap of CITE markers -- Figure 5E ####
##########################################################
# MIS-C TRBV11-2 vs non-TRBV11-2
TRBVmarker.cite.selected <- filter(TRBVmarker.cite, p_val<0.1)
TRBVmarker.cite.sig <- filter(TRBVmarker.cite, p_val_adj<0.2)

merge.CD4.UnSort.misc$TRBV[is.na(merge.CD4.UnSort.misc$TRBV)] <- "ND"
merge.CD4.UnSort.misc$TRBV112 <- merge.CD4.UnSort.misc$TRBV == "TRBV11-2"
merge.CD4.UnSort$TRBV112 <- merge.CD4.UnSort$TRBV == "TRBV11-2"
merge.CD4.UnSort$TRBV112_Class <- paste(merge.CD4.UnSort$TRBV112, merge.CD4.UnSort$Class, sep = "_")

Idents(merge.CD4.UnSort) <- "TRBV112_Class"
avg.CD4.UnSort.cite <- data.frame(log1p(AverageExpression(merge.CD4.UnSort, verbose = FALSE)$CITE))
avg.CD4.UnSort.cite$marker <- rownames(avg.CD4.UnSort.cite)

avg.CD4.UnSort.cite.marker <- filter(avg.CD4.UnSort.cite, marker %in% features)
avg.CD4.UnSort.cite.marker <- avg.CD4.UnSort.cite.marker[,c("FALSE_HC","FALSE_COVID","FALSE_MIS.C","TRUE_HC","TRUE_COVID","TRUE_MIS.C")]

avg.CD4.UnSort.cite.marker.sig <- avg.CD4.UnSort.cite.marker[rownames(avg.CD4.UnSort.cite.marker) %in% rownames(TRBVmarker.cite.sig),]
pheatmap(avg.CD4.UnSort.cite.marker.sig[,1:6], cluster_cols = FALSE,
         border_color = "white", scale = "row",
         cellwidth = 25,cellheight = 8,
         # filename = file.path("output/CD4.TRBV112.marker.cite.heatmap.padj0.2.pdf"),
         # display_numbers = display.sig,
         fontsize_row = 6, fontsize_col = 6)
dev.off()



