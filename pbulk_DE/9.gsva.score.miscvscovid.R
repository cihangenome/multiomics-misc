library(edgeR)
library(GSVA)
library(Biobase)
library(reshape2)
library(BiocParallel)
library(readr)
library(tidyverse)
library(ggpubr)
source("pbulk_DE/util/gsvadftoelist.R")

### get GSVA score of the filtered pseudobulk objects ##############################################
### using LE genes of MIS-C vs COVID
### input is output from filtered and normalized pseudobulk objects
### Generates Fig.4e
### UnSort ###
DGELISTS_IN_PATH_UNSORT <- c("output/pbulk_DE/all_samples/pseudobulk_dgelists_normalized/UnSort-mergedcelltype_pseudobulk_dgelist_normalized.rds")
DGELISTS_IN_PATH_UNFIL_UNSORT <- c("output/pbulk_DE/all_samples/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds")
OUT_DIR <- "output/pbulk_DE/gsva/"
dir.create(OUT_DIR, recursive = TRUE)

COMBINED.GENESETS.IN.PATH <- "input/kegg_go_btm_reactome_foointerferon.rds"
genesets <- readRDS(COMBINED.GENESETS.IN.PATH)
genesets_NFkB <- genesets["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]

pbulk_list <- readRDS(DGELISTS_IN_PATH_UNSORT)

# list of expression data by cell type
eset_list <- lapply(pbulk_list, function(dge){
  mat <- DESeq2::varianceStabilizingTransformation(dge$counts)
  meta <- dge$samples
  ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(meta))
})
eset_list <- unlist(eset_list)

# get fgsea results of LE genes from the MISC-COVID model
FGSEA_IN_PATH_MISC <- file.path("output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_covid/fgsea_tables/UnSort_misc_vs_covid_fgseares.rds")

fgsea_list_misc <- readRDS(FGSEA_IN_PATH_MISC)

fgsea_res_misc <- bind_rows(fgsea_list_misc, .id = "celltype") %>%
  mutate(leadingEdge = sapply(leadingEdge, toString))


scores_list <- lapply(1:length(eset_list), function(eset){
  cat(names(eset_list)[eset],"\n")
  celltype.gsea.res <- subset(fgsea_res_misc, celltype == names(eset_list)[eset])
  if (nrow(celltype.gsea.res) > 0) {
    celltype.leading.edge.genesets.misc <- sapply(celltype.gsea.res$leadingEdge,function(x){unique(unlist(strsplit(x,", ")))})

    names(celltype.leading.edge.genesets.misc) <- celltype.gsea.res$pathway
    module.scores <- gsva(expr = eset_list[[eset]], gset.idx.list = celltype.leading.edge.genesets.misc, method = "gsva", parallel.sz = 16, min.sz = 5)
    return(cbind(reshape2::melt(exprs(module.scores)),celltype= names(eset_list)[eset]))
  }
})



# add meta data and transform the dataframe into esetlist
unfiltered_elist <- readRDS(DGELISTS_IN_PATH_UNFIL_UNSORT)

module.scores.df <- do.call("rbind",scores_list)
colnames(module.scores.df) <- c("pathway","sample","module.score","celltype")
module.scores.df2 <- left_join(module.scores.df,fgsea_res_misc[,c("pathway","celltype","pval","padj")],by=c("pathway","celltype"))
module.scores.df2$pathway_celltype <- paste0(module.scores.df2$pathway,"--",module.scores.df2$celltype)
module.scores.df2 <- left_join(module.scores.df2, unfiltered_elist$CD4_Mem$samples[,c("Class","sample_id","days_since_admission","Timepoint2")], by=c("sample"="sample_id"))
write.csv(module.scores.df2, file.path(OUT_DIR, "miscvscovid_module_score_gsva_filtered_samples_genes_df.csv"))

# transfer to eset list with metadata
gsva_esetlist <- df_to_gsva_elist(gsva_df = module.scores.df, meta_eset = unfiltered_elist)
saveRDS(gsva_esetlist, file.path(OUT_DIR, "miscvscovid_module_score_gsva_filtered_samples_genes.rds"))

# plot pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# Fig. 4e
module.scores.df2 <- read.csv(file.path(OUT_DIR, "miscvscovid_module_score_gsva_filtered_samples_genes_df.csv"), row.names = 1)
celltype_selected <- c("CD4_Naive","CD4_Mem","CD8_Naive","CD8_Mem","NK_CD16hi","cDC","Mono_Classical","Mono_NonClassical")
module.scores.NFkB <- filter(module.scores.df2, pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                             days_since_admission<41,
                             Class %in% c("COVID","MIS-C"),
                             celltype %in% celltype_selected)
module.scores.NFkB$celltype <- factor(module.scores.NFkB$celltype, levels = celltype_selected)

p <- ggplot(module.scores.NFkB, aes(x=celltype, y=module.score, color = Class))+
  geom_boxplot(outlier.colour = NA)+
  scale_color_manual(values = c("COVID" = "#00AFBB", "MIS-C" = "#E7B800"))+
  geom_point(position = position_jitterdodge())+
  ylim(-0.8,1)+
  theme_bw()

p
# ggsave(p, filename = "output/plots/mono_NFkB_gsva.pdf", device = "pdf", width = 8.5, height = 6)















