library(limma)
library(edgeR)
library(Biobase)
library(tidyverse)
library(openxlsx)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(RColorBrewer)
source("pbulk_DE/util/plot_subject_heatmap.R")

## pseudobulk heatmaps
## Generates Fig.4e
DGE_IN_PATH <- "output/pbulk_DE/all_samples/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds"
dir.create("output/pbulk_DE/dge_eset_lists/", recursive = TRUE)
ESET_OUT_PATH <- "output/pbulk_DE/dge_eset_lists/pbulk_eset_list_normalized_UnSort_mergedcelltype_metafiltered.rds"
FGSEA_IN_PATH <- "output/pbulk_DE/sample_groups/"
UnSort.mergedcelltype <- readRDS(DGE_IN_PATH)

eset_list <- lapply(UnSort.mergedcelltype, function(dge){
  mat <- cpm(dge, log = TRUE)
  meta <- dge$samples
  ExpressionSet(assayData = mat, phenoData = AnnotatedDataFrame(meta))
})


parse_meta <- function(dge){
  dge@phenoData@data$Class <- factor(dge@phenoData@data$Class, levels = c("HC","COVID","MIS-C"))
    # filter <- colnames(dge)[!is.na(dge@phenoData@data$PC1_cat)]
  dge <- dge
}

eset_list <- lapply(eset_list, parse_meta)
saveRDS(eset_list, file = ESET_OUT_PATH)

####################
### heatmap plot ###
####################

eset_list <- readRDS(ESET_OUT_PATH)
eset_list_d40 <- lapply(eset_list, function(dge){
  # dge@phenoData@data$Class <- factor(dge@phenoData@data$Class, levels = c("HC","COVID","MIS-C"))
  filter <- colnames(dge)[dge@phenoData@data$days_since_admission<41]
  dge <- dge[,filter]
})

# read in selected_pathway_summary.xlsx
pathways <- read.xlsx("input/pathway_summary.xlsx", sheet = "Sheet1", startRow = 1)
# IFNsets <- filter(pathways, Category %in% c("Type I IFN response"))
NFkBsets <- filter(pathways, Primary_gene_sets %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"))

# read in FGSEA results tables
misc_days_onset <- readRDS(paste0(FGSEA_IN_PATH, "all_timepoints_misc_only/fgsea_tables/UnSort_days_onset_fgseares.rds"))
covid_days_onset <- readRDS(paste0(FGSEA_IN_PATH, "all_timepoints_covid_only/fgsea_tables/UnSort_days_onset_fgseares.rds"))
misc_vs_HC <- readRDS(paste0(FGSEA_IN_PATH, "tso_within_d40_misc_plus_healthy/fgsea_tables/UnSort_misc_vs_healthy_fgseares.rds"))
covid_vs_HC <- readRDS(paste0(FGSEA_IN_PATH, "tso_within_d40_covid_plus_healthy/fgsea_tables/UnSort_covid_vs_healthy_fgseares.rds"))
misc_vs_covid <- readRDS(paste0(FGSEA_IN_PATH, "tso_within_d40_misc_plus_covid/fgsea_tables/UnSort_misc_vs_covid_fgseares.rds"))

## UnSorted populations
### CD4_Mem NFkB ####
CD4_Mem <- eset_list$CD4_Mem
# showing the leading edge genes of MISC vs COVID
CD4_Mem_misc_vs_HC_NFkB <- misc_vs_HC$CD4_Mem %>% filter(pathway %in% NFkBsets$Primary_gene_sets)
CD4_Mem_covid_vs_HC_NFkB <- covid_vs_HC$CD4_Mem %>% filter(pathway %in% NFkBsets$Primary_gene_sets)
CD4_Mem_misc_vs_covid_NFkB <- misc_vs_covid$CD4_Mem %>% filter(pathway %in% NFkBsets$Primary_gene_sets)

CD4_Mem_misc_vs_HC_NFkB_LE <- unique(unlist(sapply(CD4_Mem_misc_vs_HC_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))
CD4_Mem_covid_vs_HC_NFkB_LE <- unique(unlist(sapply(CD4_Mem_covid_vs_HC_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))
CD4_Mem_misc_vs_covid_NFkB_LE <- unique(unlist(sapply(CD4_Mem_misc_vs_covid_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))
Mono_Classical_misc_vs_covid_NFkB_LE <- unique(unlist(sapply(Mono_Classical_misc_vs_covid_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))

CD4_Mem_mtx_co <- exprs(CD4_Mem)[CD4_Mem_misc_vs_covid_NFkB_LE, ]
meta <- pData(CD4_Mem) %>% dplyr::arrange(Class,days_since_admission)
CD4_Mem_mtx_co <- CD4_Mem_mtx_co[,rownames(meta)]

# pdf("output/pbulk_DE/plots/CD4_Mem_misc_vs_covid_NFkB.pdf", width = 8, height = 4)
subject.hm2(exprs_mtx = CD4_Mem_mtx_co, meta = meta, celltype = "CD4_Mem", module = "NFkB", 
            LEset1 = Mono_Classical_misc_vs_covid_NFkB_LE, LEset2 = CD4_Mem_misc_vs_covid_NFkB_LE)
dev.off()



### Mono_Classical NFkB ####
Mono_Classical <- eset_list$Mono_Classical
# showing the leading edge genes of MISC vs COVID
Mono_Classical_misc_vs_HC_NFkB <- misc_vs_HC$Mono_Classical %>% filter(pathway %in% NFkBsets$Primary_gene_sets)
Mono_Classical_covid_vs_HC_NFkB <- covid_vs_HC$Mono_Classical %>% filter(pathway %in% NFkBsets$Primary_gene_sets)
Mono_Classical_misc_vs_covid_NFkB <- misc_vs_covid$Mono_Classical %>% filter(pathway %in% NFkBsets$Primary_gene_sets)

Mono_Classical_misc_vs_HC_NFkB_LE <- unique(unlist(sapply(Mono_Classical_misc_vs_HC_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))
Mono_Classical_covid_vs_HC_NFkB_LE <- unique(unlist(sapply(Mono_Classical_covid_vs_HC_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))
Mono_Classical_misc_vs_covid_NFkB_LE <- unique(unlist(sapply(Mono_Classical_misc_vs_covid_NFkB$leadingEdge, function(x) str_split(x, pattern = ","))))

Mono_Classical_mtx_co <- exprs(Mono_Classical)[Mono_Classical_misc_vs_covid_NFkB_LE, ]
meta <- pData(Mono_Classical) %>% dplyr::arrange(Class,days_since_admission)
Mono_Classical_mtx_co <- Mono_Classical_mtx_co[,rownames(meta)]

# pdf("output/pbulk_DE/plots/Mono_Classical_misc_vs_covid_NFkB.pdf", width = 8, height = 4)
subject.hm2(exprs_mtx = Mono_Classical_mtx_co, meta = meta, celltype = "Mono_Classical", module = "NFkB", 
            LEset1 = Mono_Classical_misc_vs_covid_NFkB_LE, LEset2 = CD4_Mem_misc_vs_covid_NFkB_LE)
dev.off()









