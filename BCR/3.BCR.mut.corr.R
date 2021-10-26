library(Seurat) #load Seurat 3.1
library(matrixStats)
library(plyr)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(Hmisc)
####################################################################
### Correlation of mutation freq with surface protein expression ###
####################################################################
# Generates Figures S7E, S7F

### BCR and gene expression
merge <- readRDS("input/misc_SeuratObj_submission.rds")
merge.B <- subset(merge, subset = coarsecelltype %in% c("B","Plasmablast"))
meta <- merge.B@meta.data
meta$Class <- factor(meta$Class, levels = c("HC","COVID","MIS-C"))

# within BCR clustering
NormalizeData(merge.B, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
merge.B <- FindVariableFeatures(merge.B)
merge.B <- ScaleData(merge.B, features = VariableFeatures(merge.B))
merge.B <- RunPCA(merge.B, features = VariableFeatures(merge.B))
merge.B <- FindNeighbors(merge.B, reduction = "pca", dims = 1:15)
merge.B <- FindClusters(merge.B, resolution = 0.6, verbose = FALSE)
merge.B <- RunUMAP(merge.B, reduction = "pca", dims = 1:15)
DimPlot(merge.B, group.by = "seurat_clusters", label = TRUE)

### Correlation of mutation freq with surface protein expression ###
#######################
### Memory B cells ####
#######################
merge.B.mem <- subset(merge.B, subset = WCT.coarsecelltype == "B_Mem")

merge.B.cite <- data.frame(t(data.frame(merge.B.mem@assays$CITE@data))) %>%
  mutate("mu_freq_tot" = merge.B.mem@meta.data$mu_freq_tot) %>%
  drop_na()

# use Hmisc::rcorr to get correlations
merge.B.cite.cor <- Hmisc::rcorr(as.matrix(merge.B.cite), merge.B.cite$mu_freq_tot)
merge.B.cite.cor <- cbind(merge.B.cite.cor$r[,c("y")], merge.B.cite.cor$P[,c("y")])
colnames(merge.B.cite.cor) <- c("cor.r","cor.pval")
merge.B.cite.cor <- data.frame(merge.B.cite.cor) %>%
  rownames_to_column("marker") %>%
  arrange(desc(abs(cor.r))) %>%
  column_to_rownames("marker")
rownames(merge.B.cite.cor) <- gsub("\\.1", "", rownames(merge.B.cite.cor))


selected.marker <- c("IgD","CD82","CD305","CD27","CD99","IgM","CD49d","CD54","CD45RA","CD35",
                     "HLA.ABC","CD95","CD71","IgA","IgGFc","HLA.DR","S1probe")
merge.B.cite.cor.plot <- data.frame("cor.r" = merge.B.cite.cor[selected.marker,"cor.r"],
                                    "marker" = selected.marker)

p <- ggplot(merge.B.cite.cor.plot, aes(x=reorder(marker, abs(cor.r)), y=cor.r)) +
  geom_point(color=ifelse(merge.B.cite.cor.plot$cor.r>0, "#E64B35FF", "#4DBBD5FF"), size = 3) + 
  geom_segment(aes(x=marker, xend=marker, y=0, yend=cor.r), color=ifelse(merge.B.cite.cor.plot$cor.r>0, "#E64B35FF", "#4DBBD5FF"))+
  coord_flip() + ggtitle("Mem_B mutation correlation") +
  xlab("CITE_Markers")+
  theme_bw()
ggsave(p, device = "pdf", filename = "output/BCR_MemB_mut_corr_cite.pdf", width = 4, height = 4)

######################
### Plasmblasts ######
######################
merge.B.PB <- subset(merge.B, subset = WCT.coarsecelltype == "Plasmablast")

merge.B.cite <- data.frame(t(data.frame(merge.B.PB@assays$CITE@data))) %>%
  mutate("mu_freq_tot" = merge.B.PB@meta.data$mu_freq_tot) %>%
  # filter(mu_freq_tot > 0) %>%
  drop_na()

# use Hmisc::rcorr to get correlations
merge.B.cite.cor <- Hmisc::rcorr(as.matrix(merge.B.cite), merge.B.cite$mu_freq_tot)
merge.B.cite.cor <- cbind(merge.B.cite.cor$r[,c("y")], merge.B.cite.cor$P[,c("y")])
colnames(merge.B.cite.cor) <- c("cor.r","cor.pval")
merge.B.cite.cor <- data.frame(merge.B.cite.cor) %>%
  rownames_to_column("marker") %>%
  arrange(desc(abs(cor.r))) %>%
  column_to_rownames("marker")
rownames(merge.B.cite.cor) <- gsub("\\.1", "", rownames(merge.B.cite.cor))

# showing top10 markers
selected.marker <- c("CD95","CD99","CD58","CD18","CD22","CD146","CD360","HLA.DR",
                     "CD45RA","NLRP2")
merge.B.cite.cor.plot <- data.frame("cor.r" = merge.B.cite.cor[selected.marker,"cor.r"],
                                    "marker" = selected.marker)

p <- ggplot(merge.B.cite.cor.plot, aes(x=reorder(marker, abs(cor.r)), y=cor.r)) +
  geom_point(color=ifelse(merge.B.cite.cor.plot$cor.r>0, "#E64B35FF", "#4DBBD5FF"), size = 3) + 
  geom_segment(aes(x=marker, xend=marker, y=0, yend=cor.r), color=ifelse(merge.B.cite.cor.plot$cor.r>0, "#E64B35FF", "#4DBBD5FF"))+
  coord_flip() + ggtitle("Plasmablast mutation correlation") +
  xlab("CITE_Markers")+
  theme_bw()
ggsave(p, device = "pdf", filename = "output/BCR_PB_mut_corr_cite.pdf", width = 4, height = 4)




