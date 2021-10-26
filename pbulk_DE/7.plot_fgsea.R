library(tidyverse)
library(openxlsx)
source("pbulk_DE/util/plot_fgsea_functions.R")

#########################################
#### Plot fgsea #########################
#########################################
## generates Figures 4B, 4C, S5D
INPUT_PATH <- "pbulk_DE/sample_groups/"

SAMPLE_GROUP <- c("tso_within_d40_misc_plus_healthy/results/", 
                  "tso_within_d40_covid_plus_healthy/results/",
                  "tso_within_d40_misc_plus_covid/results/", 
                  "all_timepoints_misc_only/results/")


selected.pathways <- read.xlsx("input/pathway_summary.xlsx", sheet = "Sheet1")
celltypes <- c("B_Naive", "B_Mem", "Plasmablast", "CD4_Naive", "CD4_Mem", "CD4_isobinding",
               "CD8_Naive", "CD8_Mem", "MAIT", "T_Vd2", "NKT", "DNT",
               "NK_CD16hi", "NK_CD56hiCD16lo", 
               "Mono_Classical", "Mono_NonClassical", "Mono_Intermediate",
               "cDC", "pDC", "HSC", "Platelet")
cellgroups = rep(c("B","CD4","CD8","otherT","NK","Mono","other"),c(3,3,2,4,2,3,4))
cellgroup <- data.frame(celltype = celltypes, 
                        cellgroup = cellgroups)

#########################
### MIS-C vs HC d40 #####
#########################
fgseares.list.Unsort <- readRDS(paste0(INPUT_PATH, SAMPLE_GROUP[1], "fgsea_tables/UnSort_misc_vs_healthy_fgseares.rds"))
fgseares.Unsort <- RbindGseaResultList(fgseares.list.Unsort, NES_filter = 0, padj_filter = 1)  %>%
  filter(celltype %in% celltypes) %>%
  left_join(cellgroup, by = c("celltype")) %>%
  mutate(celltype = factor(.$celltype, levels = celltypes)) %>%
  mutate(cellgroup = factor(.$cellgroup, levels = unique(cellgroups))) %>%
  mutate(n_logp = replace(n_logp, padj > .2, NA))

fgseares.selected <- fgseares.Unsort %>% 
  filter(pathway %in% selected.pathways$Primary_gene_sets) %>% 
  left_join(selected.pathways, by = c("pathway" = "Primary_gene_sets"))

pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[1], "fgsea_bubble/UnSort_misc_vs_healthy_fgseares_selected.pdf"), width = 10, height = 8)
GSEABubblePlot_selected(fgseares.selected)
dev.off()


#########################
### COVID vs HC d40 #####
#########################
fgseares.list.Unsort <- readRDS(paste0(INPUT_PATH, SAMPLE_GROUP[2], "fgsea_tables/UnSort_covid_vs_healthy_fgseares.rds"))
fgseares.Unsort <- RbindGseaResultList(fgseares.list.Unsort, NES_filter = 0, padj_filter = 1)  %>%
  filter(celltype %in% celltypes) %>%
  left_join(cellgroup, by = c("celltype")) %>%
  mutate(celltype = factor(.$celltype, levels = celltypes)) %>%
  mutate(cellgroup = factor(.$cellgroup, levels = unique(cellgroups))) %>%
  mutate(n_logp = replace(n_logp, padj > .2, NA))

fgseares.selected <- fgseares.Unsort %>% 
  filter(pathway %in% selected.pathways$Primary_gene_sets) %>% 
  left_join(selected.pathways, by = c("pathway" = "Primary_gene_sets"))

pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[2], "fgsea_bubble/UnSort_covid_vs_healthy_fgseares_selected.pdf"), width = 10, height = 8)
GSEABubblePlot_selected(fgseares.selected)
dev.off()



##########################
### MIS-C vs COVID d40 ###
##########################
fgseares.list.Unsort <- readRDS(paste0(INPUT_PATH, SAMPLE_GROUP[3], "fgsea_tables/UnSort_misc_vs_covid_fgseares.rds"))
fgseares.Unsort <- RbindGseaResultList(fgseares.list.Unsort, NES_filter = 0, padj_filter = 1)  %>%
  filter(celltype %in% celltypes) %>%
  left_join(cellgroup, by = c("celltype")) %>%
  mutate(celltype = factor(.$celltype, levels = celltypes)) %>%
  mutate(cellgroup = factor(.$cellgroup, levels = unique(cellgroups))) %>%
  mutate(n_logp = replace(n_logp, padj > .2, NA))

fgseares.selected <- fgseares.Unsort %>% 
  filter(pathway %in% selected.pathways$Primary_gene_sets) %>% 
  left_join(selected.pathways, by = c("pathway" = "Primary_gene_sets"))

pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[3], "fgsea_bubble/UnSort_misc_vs_covid_fgseares_selected.pdf"), width = 10, height = 8)
GSEABubblePlot_selected(fgseares.selected)
dev.off()



#######################################
### MIS-C only all timepoints #########
#######################################
fgseares.list.Unsort <- readRDS(paste0(INPUT_PATH, SAMPLE_GROUP[4], "fgsea_tables/UnSort_days_onset_fgseares.rds"))
fgseares.Unsort <- RbindGseaResultList(fgseares.list.Unsort, NES_filter = 0, padj_filter = 1)  %>%
  filter(celltype %in% celltypes) %>%
  left_join(cellgroup, by = c("celltype")) %>%
  mutate(celltype = factor(.$celltype, levels = celltypes)) %>%
  mutate(cellgroup = factor(.$cellgroup, levels = unique(cellgroups))) %>%
  mutate(n_logp = replace(n_logp, padj > .2, NA))

fgseares.selected <- fgseares.Unsort %>% 
  filter(pathway %in% selected.pathways$Primary_gene_sets) %>% 
  left_join(selected.pathways, by = c("pathway" = "Primary_gene_sets"))

pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[4], "fgsea_bubble/UnSort_days_onset_fgseares_selected.pdf"), width = 10, height = 8)
GSEABubblePlot_selected(fgseares.selected)
dev.off()



#######################################
### COVID only all timepoints #########
#######################################
fgseares.list.Unsort <- readRDS(paste0(INPUT_PATH, SAMPLE_GROUP[5], "fgsea_tables/UnSort_days_onset_fgseares.rds"))
fgseares.Unsort <- RbindGseaResultList(fgseares.list.Unsort, NES_filter = 0, padj_filter = 1)  %>%
  filter(celltype %in% celltypes) %>%
  left_join(cellgroup, by = c("celltype")) %>%
  mutate(celltype = factor(.$celltype, levels = celltypes)) %>%
  mutate(cellgroup = factor(.$cellgroup, levels = unique(cellgroups))) %>%
  mutate(n_logp = replace(n_logp, padj > .2, NA))

fgseares.selected <- fgseares.Unsort %>% 
  filter(pathway %in% selected.pathways$Primary_gene_sets) %>% 
  left_join(selected.pathways, by = c("pathway" = "Primary_gene_sets"))

pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[5], "fgsea_bubble/UnSort_days_onset_fgseares_selected.pdf"), width = 10, height = 8)
GSEABubblePlot_selected(fgseares.selected)
dev.off()


### Figure 4C
# plot IFN time change for MISC and COVID, showing all cell populations
fgseares.selected.IFN.misctime <- dplyr::full_join(fgseares.selected.IFN.misctime, fgseares.selected.IFN.covidtime[,c(1:2,6:7)], 
                                                   by = c("celltype","pathway","cellgroup","Category"))
fgseares.selected.IFN.covidtime <- dplyr::full_join(fgseares.selected.IFN.covidtime, fgseares.selected.IFN.misctime[,c(1:2,6:7)], 
                                                   by = c("celltype","pathway","cellgroup","Category"))
pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[4], "fgsea_bubble/UnSort_days_onset_fgseares_selectedIFN.pdf"), width = 10, height = 2.3)
GSEABubblePlot_selected(fgseares.selected.IFN.misctime)
dev.off()

pdf(file = paste0(INPUT_PATH, SAMPLE_GROUP[5], "fgsea_bubble/UnSort_days_onset_fgseares_selectedIFN.pdf"), width = 10, height = 2.3)
GSEABubblePlot_selected(fgseares.selected.IFN.covidtime)
dev.off()




