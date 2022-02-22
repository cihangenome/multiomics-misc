# this is tested using R 3.6.1 on a high-performance computing node with 8 cores and 160 gb of ram. 
# Use seurat v3.1.0 on Conda environment
library(Seurat) #load Seurat 3.1
library(matrixStats)
library(tidyverse)
library(pheatmap)
library(viridis)
library(reshape2)
library(openxlsx)
library(fgsea)

#################################
### TRBV11-2 CD4 T cell GSEA ####
#################################
# Generates Extended Data Fig.6f
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

# Parse meta
merge.CD4.UnSort$TRBV[is.na(merge.CD4.UnSort$TRBV)] <- "ND"
merge.CD4.UnSort$TRBV_Class = paste(merge.CD4.UnSort$TRBV, merge.CD4.UnSort$Class, sep = "_")

# MIS-C only
merge.CD4.UnSort.misc <- subset(merge.CD4.UnSort, subset = Class == "MIS-C")
Idents(merge.CD4.UnSort.misc) = "TRBV"
# use no cutoff to get the gene list for GSEA
TRBVmarker.misc.rna = FindMarkers(merge.CD4.UnSort.misc, ident.1 = "TRBV11-2", assay = "RNA", logfc.threshold = 0, min.pct = 0)
saveRDS(TRBVmarker.misc.rna, "output/TRBV11-2.misc.CD4.DEallmarker.wocutoff.rds")

###################################################################################################
### enrichment test on the makers from single cell level DE analysis using Seurat FindMarkers() ###
###################################################################################################
### use fgsea rank based test
DElist.wocutoff <- readRDS("output/TRBV11-2.misc.CD4.DEallmarker.wocutoff.rds")
DElist.wocutoff$gene <- rownames(DElist.wocutoff)
geneset = readRDS("input/kegg_go_btm_reactome_foointerferon.rds")

set.seed(1)
ranks.tmp <- DElist.wocutoff$avg_logFC
names(ranks.tmp) <- DElist.wocutoff$gene
ranks.tmp <- ranks.tmp[order(ranks.tmp, decreasing = TRUE)]

fgseaRes <- fgsea(geneset, ranks.tmp, minSize=10, maxSize = 1500, nperm=10000)

fgseaRes.df <- data.frame(fgseaRes) %>% mutate(celltype = "TRBV11-2_CD4")

fgseaRes.df.filtered <- fgseaRes.df %>% filter(padj < 0.2) %>% 
  mutate(n_logp = -log10(padj)) 
fgseaRes.df$leadingEdge <- vapply(fgseaRes.df$leadingEdge, paste, collapse = ", ", character(1L))
fgseaRes.df.filtered$leadingEdge <- vapply(fgseaRes.df.filtered$leadingEdge, paste, collapse = ", ", character(1L))

write.csv(fgseaRes.df.filtered, "output/TRBV11-2.misc.CD4.fgsea.enrich.padj0.2.csv")

fgseaRes.df.apoptotic <- fgseaRes.df.filtered %>% 
  filter(grepl("apop|APOP|Apop", pathway)) %>%
  filter(pathway %in% c("reactome_Regulation of Apoptosis", "reactome_Apoptosis", 
                        "GO_POSITIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY", 
                        "GO_REGULATION_OF_PROTEIN_INSERTION_INTO_MITOCHONDRIAL_MEMBRANE_INVOLVED_IN_APOPTOTIC_SIGNALING_PATHWAY",
                        "GO_REGULATION_OF_MITOCHONDRIAL_OUTER_MEMBRANE_PERMEABILIZATION_INVOLVED_IN_APOPTOTIC_SIGNALING_PATHWAY",
                        "GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
                        "GO_APOPTOTIC_SIGNALING_PATHWAY"))

# pdf(file = "TCR/plots/TRBV11-2.misc.CD4.fgsea.enrichment.apoptotic.pdf", width = 11, height = 3)
# GSEABubblePlot(fgseaRes.df.apoptotic)
# dev.off()

# Extended Data Fig.6f
p=ggplot(fgseaRes.df.apoptotic, aes(NES,forcats::fct_reorder(pathway,NES)))+
  geom_segment(aes(xend=0,yend=pathway, color=NES), size = 1)+
  geom_point(aes(color=NES,size=n_logp))+
  # scale_size(range = c(0, 7))+
  scale_colour_gradientn(limits = c(0.5,2), colours=c("white","red3"))+
  #scale_color_viridis_c(begin = 0.5, end = 1)+
  # scale_size_continuous(range =c(2,8))+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
  theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
  ylab(NULL)
ggsave(plot = p, file = "output/plots/TRBV11-2.misc.CD4.fgsea.enrichment.apoptotic.bar.pdf", useDingbats = FALSE, width = 11, height = 3)



















