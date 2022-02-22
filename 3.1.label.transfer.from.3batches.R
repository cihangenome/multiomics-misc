library(Seurat)
library(matrixStats)
library(plyr)
library(tidyverse)
library(pheatmap)
library(viridis)
library(reshape2)

### label transfer #######################################################################
# transfer WCT labels from COVID - Brescia batches 1-3
# download the seurat object of adult data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161918
Brescia_merge <- readRDS("input/brescia_paper1_seurat.rds")
Brescia_merge <- NormalizeData(Brescia_merge)
Brescia_merge <- FindVariableFeatures(Brescia_merge, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
merge <- readRDS("misc_SeuratObj_submission.rds")
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 2000, verbose = TRUE)

# integrate 2 datasets together
COVID.list <- list("batches1_3" = Brescia_merge, "batch4" = merge)

merge.anchors <- FindTransferAnchors(reference = Brescia_merge, query = merge, 
                                     dims = 1:30)
predictions <- TransferData(anchorset = merge.anchors, refdata = Brescia_merge$WCTcoursecelltype, 
                            dims = 1:30)

merge <- AddMetaData(merge, metadata = predictions)

# plot the overlap of predicted transfer id with annotated id
## Extended Data Fig.5a
celltype <- data.frame(table(merge$mergedcelltype, merge$predicted.id)) %>%
  left_join(data.frame(table(merge$mergedcelltype)), by = "Var1") %>%
  mutate(prct = Freq.x/Freq.y*100) %>%
  reshape2::dcast(Var1 ~ Var2, value.var = "prct") %>%
  column_to_rownames("Var1")

pdf("celltype.overlapwWCTcourse.pdf", width = 5.5, height = 4.5)
pheatmap(celltype,cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

saveRDS(merge, "input/B4merge.Clustered.COVID.QCd.rds")



