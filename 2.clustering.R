#this is tested using R 3.6.1 on a high-performance computing node with 8 cores and at least 160 gb of ram. 
.libPaths(c("/hpcdata/sg/sg_data/users/liuc19/condaenvs/Seurat3.1/lib/R/library",
            "/sysapps/cluster/software/Anaconda2/5.3.0/lib64/R/site-library"))
library("Seurat") #load Seurat 3.1
library("dplyr")
library("matrixStats")
library('tidyverse')

#clustering based on protein (CITE) data

B4merge <- readRDS(file = "B4merge.SNG.preclust.rds")
CITEB4merge = (GetAssayData(B4merge[["CITE"]], slot = "data"))
library('parallelDist')
adt.dist <- parDist(t(CITEB4merge[rownames(B4merge[["CITE"]]),]), threads = 8)

B4merge[["adt_snn"]] <- FindNeighbors(adt.dist, nn.eps = 1)$snn
B4merge <- FindClusters(B4merge, resolution = c(2), graph.name = "adt_snn", algorithm = 1)
#Plot of different cell clusters identified based on protein expression levels
B4merge <- RunUMAP(B4merge, assay = "CITE", features = rownames(B4merge[["CITE"]]), n_neighbors=50L, min_dist=1)

pdf("init.clust.umap.pdf", width = 20, height = 18)
# DimPlot(B4merge, do.return = TRUE, pt.size = 0.5, reduction.use = "umap", group.by = "seurat_clusters", vector.friendly = TRUE, label=TRUE, no.legend=TRUE)
DimPlot(B4merge, pt.size = 0.5, reduction = "umap", group.by = "adt_snn_res.2", label=TRUE)
dev.off()
saveRDS(B4merge, file = "B4merge.SNG.postclust.rds")

table(B4merge$Batch, B4merge$adt_snn_res.2)
table(B4merge$seurat_clusters)
source("util_funs/AverageExpression_MeanOnly.r")
environment(AverageExpression_MeanOnly) <- asNamespace('Seurat') # need this to allow the function to call hidden seurat package functions
aver = AverageExpression_MeanOnly(B4merge, return.seurat=T)

quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(GetAssayData(aver[["CITE"]])), n = 101)

library("pheatmap")
library("viridis")
ggsave(filename = "averExp_res.2_CITE_infernoColor_quantileBreaks.pdf", width = 10, height = 24,
       plot = pheatmap(GetAssayData(aver[["CITE"]]), scale = "none", border_color=NA, inferno(length(mat_breaks) - 1), breaks = mat_breaks))
ggsave(filename = "pcasnn_CITE_infernoColor_quantileBreaks.pdf", width = 10, height = 24,
       plot = pheatmap(GetAssayData(aver[["CITE"]]), scale = "none", border_color=NA, inferno(length(mat_breaks) - 1), breaks = mat_breaks))
ggsave(filename = "Vln.nCountCITE.pdf", width = 10, height = 6,
       plot = VlnPlot(B4merge, features = "nCount_CITE", pt.size = 0)
)
ggsave(filename = "Vln.nCountRNA.pdf", width = 10, height = 6,
       plot = VlnPlot(B4merge, features = "nCount_RNA", pt.size = 0)
)

source("util_funs/histogram.R")

FeatureScatter(subset(B4merge, subset = seurat_clusters %in% c("2", "22", "27")), feature1="CD56", feature2="CD16")

DimPlot(B4merge, group.by = "Class", label = TRUE)
FeaturePlot(B4merge, reduction = "umap",
            features = c("cite_CD19.1", "cite_S1probe", "cite_CD20", "cite_TCRVd2"), 
            min.cutoff = "q05", max.cutoff = "q95", ncol = 2)

write.csv(file = "B4merge_CellsPerSamplebyCluster.csv", table(B4merge$adt_snn_res.2, B4merge$sample_id))

# put in cluster annotations
ClusterNames = read.csv("input/ClusterNames.B4.20201006.csv", header=TRUE)

Idents(B4merge) <- "adt_snn_res.2"
cluster.ids <- ClusterNames$cluster
B4merge$hirescelltype <- plyr::mapvalues(x = Idents(B4merge), from = cluster.ids, to = as.character(ClusterNames$hires_celltype))
B4merge$mergedcelltype <- plyr::mapvalues(x = Idents(B4merge), from = cluster.ids, to = as.character(ClusterNames$mergedcelltype))
B4merge$coarsecelltype <- plyr::mapvalues(x = Idents(B4merge), from = cluster.ids, to = as.character(ClusterNames$coarsecelltype))

B4merge_sort <- subset(B4merge, subset = Batch == "B4chip3")
B4merge_unt <- subset(B4merge, subset = Batch == "B4chip3", invert = TRUE)

Idents(B4merge_sort) <- "hirescelltype"
pdf("umap.labelled.sort.pdf", width = 10, height = 8)
DimPlot(B4merge_sort, pt.size = 0.5, reduction = "umap", label = TRUE) + NoLegend()
dev.off()
Idents(B4merge_unt) <- "hirescelltype"
pdf("umap.labelled.unt.pdf", width = 14, height = 12)
DimPlot(B4merge_unt, pt.size = 0.5, reduction = "umap", label = TRUE) + NoLegend()
dev.off()
Idents(B4merge) <- "hirescelltype"
pdf("umap.labelled.all.pdf", width = 8, height = 6)
DimPlot(B4merge, pt.size = 0.5, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

table(B4merge_unt$hirescelltype)
table(B4merge_unt$mergedcelltype)
table(B4merge_unt$coarsecelltype)


# B4merge <- DietSeurat(B4merge)
saveRDS(B4merge, file = "B4merge.Clustered.20201006.rds")
B4merge <- readRDS("B4merge.Clustered.20201006.rds")


### Cluster QC: celltypes to remove contamination cells ----------------------------------------------------------------------------------------------------------------------------
Idents(B4merge) <- "coarsecelltype"
FeatureScatter(subset(B4merge, subset = coarsecelltype == c("Mono")), feature1 = "cite_CD10", feature2 = "cite_CD24")
B4merge_QCd <- subset(B4merge, cells = colnames(B4merge[["CITE"]][,which(B4merge$coarsecelltype == c("Mono") & B4merge[["CITE"]]["CD19.1",] > 10)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("Mono") & B4merge_QCd[["CITE"]]["CD3",] > 15)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("Mono") & B4merge_QCd[["CITE"]]["CD10",] > 8)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("Mono") & B4merge_QCd[["CITE"]]["CD24",] > 6)]), invert = TRUE)

FeatureScatter(subset(B4merge_QCd, subset = coarsecelltype %in% c("Plasmablast", "B")), feature1 = "cite_CD3", feature2 = "cite_CD56")
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("B") & B4merge_QCd[["CITE"]]["CD16",] > 10)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("B") & B4merge_QCd[["CITE"]]["CD64",] > 10)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("B") & B4merge_QCd[["CITE"]]["CD3",] > 15)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("B") & B4merge_QCd[["CITE"]]["CD56",] > 8)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("Plasmablast") & B4merge_QCd[["CITE"]]["CD3",] > 15)]), invert = TRUE)

FeatureScatter(subset(B4merge_QCd, subset = coarsecelltype %in% c("CD4")), feature1 = "cite_CD8", feature2 = "cite_CD19.1")
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD4") & B4merge_QCd[["CITE"]]["CD16",] > 8)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD4") & B4merge_QCd[["CITE"]]["CD64",] > 5)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD4") & B4merge_QCd[["CITE"]]["CD8",] > 5)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD4") & B4merge_QCd[["CITE"]]["CD19.1",] > 6)]), invert = TRUE)

FeatureScatter(subset(B4merge_QCd, subset = coarsecelltype %in% c("CD8")), feature1 = "cite_CD4.1", feature2 = "cite_CD19.1")
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD8") & B4merge_QCd[["CITE"]]["CD64",] > 5)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD8") & B4merge_QCd[["CITE"]]["CD14.1",] > 6)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD8") & B4merge_QCd[["CITE"]]["CD4.1",] > 15)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("CD8") & B4merge_QCd[["CITE"]]["CD19.1",] > 6)]), invert = TRUE)

FeatureScatter(subset(B4merge_QCd, subset = coarsecelltype %in% c("NK")), feature1 = "cite_CD64", feature2 = "cite_CD14.1")
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("NK") & B4merge_QCd[["CITE"]]["CD64",] > 5)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("NK") & B4merge_QCd[["CITE"]]["CD14.1",] > 5)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("NK") & B4merge_QCd[["CITE"]]["CD3",] > 8)]), invert = TRUE)
B4merge_QCd <- subset(B4merge_QCd, cells = colnames(B4merge_QCd[["CITE"]][,which(B4merge_QCd$coarsecelltype == c("NK") & B4merge_QCd[["CITE"]]["CD19.1",] > 6)]), invert = TRUE)

B4merge <- B4merge_QCd

# filter the seurat object containing only COVID data
B4merge <- subset(B4merge, subset = Class %in% c("HC","COVID","MIS-C"))
saveRDS(B4merge, file = "B4merge.Clustered.COVID.QCd.rds")

### add more metadata
####################################################
### add ethnicity and original sample_code meta ####
####################################################
merge <- readRDS("B4merge.Clustered.COVID.QCd.rds")
# input meta file downloaded from the uploaded data
samplemeta_new <- read.xlsx("input/MIS-C_cohort_citeseq", sheet = "Sheet1")
samplemeta_new$sample_id <- paste(samplemeta_new$Subject.ID, samplemeta_new$Visit, sep = "_")
samplemeta_new$sample_code <- paste(samplemeta_new$Subject_code, samplemeta_new$Visit, sep = "_")
identical(sort(samplemeta_new$sample_id), sort(unique(merge$sample_id)))

Idents(merge) <- "sample_id"
merge$subject_code <- plyr::mapvalues(x = Idents(merge), from = samplemeta_new$sample_id, to = samplemeta_new$Subject_code)
merge$sample_code <- plyr::mapvalues(x = Idents(merge), from = samplemeta_new$sample_id, to = samplemeta_new$sample_code)
merge$Ethnicity <- plyr::mapvalues(x = Idents(merge), from = samplemeta_new$sample_id, to = samplemeta_new$Ethnicity)


# remove extra columns
merge$days_since_symptom_onset <- NULL
merge$Class <- droplevels(merge$Class)
merge$HTO_maxID <- NULL
merge$HTO_secondID <- NULL
merge$HTO_margin <- NULL
merge$HTO_classification <- NULL
merge$HTO_classification.global <- NULL
merge$TCR_both_chain_exist <- NULL
saveRDS(merge, "misc_SeuratObj_submission.rds")





