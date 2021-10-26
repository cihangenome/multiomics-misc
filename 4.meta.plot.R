library("Seurat")
library("dplyr")
library("matrixStats")
library('tidyverse')
library(RColorBrewer)
library(openxlsx)

#########################################
### Plot UMAP of Final Seurat Object ####
#########################################
merge <- readRDS("input/misc_SeuratObj_submission.rds")

### Figure 4A
color_clusters<-c(brewer.pal(n = 8,name = "Set1"),
                  brewer.pal(n = 8,name = "Set2"),
                  brewer.pal(n = 12,name = "Set3"),
                  brewer.pal(n = 12,name = "Paired"),
                  brewer.pal(n = 9,name = "Pastel1"),
                  brewer.pal(n = 8,name = "Pastel2"),
                  brewer.pal(n = 8,name = "Accent"))


pdf(file = "output/UMAP_all.pdf", width = 12, height = 8)
DimPlot(object = merge,
        reduction = 'umap',
        group.by="mergedcelltype",
        label = TRUE) + scale_color_manual(values=color_clusters)
dev.off()


### within mono clustering #################################################################
DimPlot(merge, pt.size = 0.2, reduction = "umap", label = TRUE, group.by = "mergedcelltype")
DimPlot(merge, pt.size = 0.2, reduction = "umap", label = TRUE, group.by = "Sorted")

merge.mono <- subset(merge, subset = coarsecelltype %in% c("Mono"))

NormalizeData(merge.mono, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
merge.mono <- FindVariableFeatures(merge.mono)
merge.mono <- ScaleData(merge.mono, features = VariableFeatures(merge.mono))
merge.mono <- RunPCA(merge.mono, features = VariableFeatures(merge.mono))
merge.mono <- FindNeighbors(merge.mono, reduction = "pca", dims = 1:15)
merge.mono <- FindClusters(merge.mono, resolution = 0.6, verbose = FALSE)
merge.mono <- RunUMAP(merge.mono, reduction = "pca", dims = 1:15)

FeaturePlot(merge.mono, reduction = "umap",
            features = c("cite_CD163.1", "rna_S100A1", "cite_CD127"), 
            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
merge.mono.unsort <- subset(merge.mono, subset = Sorted == "UnSort")
DimPlot(merge.mono.unsort, pt.size = 0.2, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

### Figure S5C
pdf("output/UMAP.mono.Class.pdf", height = 6, width = 7.5)
DimPlot(merge.mono.unsort, pt.size = 0.2, reduction = "umap", label = FALSE, group.by = "Class")+
  scale_color_manual(values = c("HC" = "#0000AC", "COVID" = "#00AFBB", "MIS-C" = "#E7B800"))
dev.off()

merge.mono.unsort$Class <- factor(merge.mono.unsort$Class, levels = c("HC","COVID","MIS-C"))
# showing only s100 family genes with precentage expressed > 60%
s100family <- c("S100A4","S100A6","S100A8","S100A9","S100A10","S100A11","S100A12")
pdf("output/mono.scgex.group.pdf", height = 4, width = 5.5)
DotPlot(merge.mono.unsort, features = s100family, group.by = "Class") + RotatedAxis()
dev.off()
pdf("output/mono.scCD163.group.pdf", height = 4, width = 5.5)
VlnPlot(merge.mono.unsort, features = c("cite_CD163.1"), slot = "data", group.by = "Class", pt.size = 0.1)+
  scale_fill_manual(values = c("HC" = "#0000AC", "COVID" = "#00AFBB", "MIS-C" = "#E7B800"))
dev.off()










