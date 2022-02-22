library(tidyverse)
library(scRepertoire)
library(reshape2)
library(ggpubr)
library(ggsci)
library(ggfortify)
source("TCR/util/TCR_usage_functions.R")

############################################################
### MISC TCR analysis ######################################
############################################################
# Generates Fig.5c for CD4 TCR TRBV11-2 usage
# read in TCR data from output folder, output from 1.compile.TCR.data.R script
TCR_combined <- readRDS("output/TCR_tenx_filtered_anno_wofiltered_celltype_combined.rds")
TCR_combined_filtered <- filter(TCR_combined, !is.na(Class), coarsecelltype %in% c("CD4","CD8","MAIT","T_gd","T_Vd2"))

# read in metadata
meta_subject <- read.csv("input/B4sample.manualhto.csv", row.names = 1) %>%
  filter(Class %in% c("HC","COVID","MIS-C"), Pool != "misassigned")
meta_subject$sample_id <- paste(meta_subject$Subject, meta_subject$Timepoint, sep = "_")
meta_subject <- meta_subject[!duplicated(meta_subject$sample_id),]
meta_subject$Class <- factor(meta_subject$Class, levels = c("HC","COVID","MIS-C"))


#######################################################
### V gene usage ######################################
#######################################################

##################################
# cal ratio -- total cells #######
##################################

# TRBV gene
# days_since_admission within d40
TCR_combined_filtered_unsort <- filter(TCR_combined_filtered, Sorted == "UnSort")
TCR_combined_filtered_unsort_T40 <- filter(TCR_combined_filtered_unsort, days_since_admission<41)

TCR_combined_filtered_TRBV_subject <- cal_Vgene_subject_ratio(TCR_combined_filtered_unsort_T40, "TRBV")
p <- plot_Vgene_usage_ratio(TCR_combined_filtered_TRBV_subject, "all_T_cells", "TRBV")
# ggsave(p, device = "pdf", filename = "output/plots/TCR_allTcells_TRBVgene_T40.pdf", width = 25, height = 17.5)

# time effect
TCR_combined_filtered_TRBV_subject <- cal_Vgene_subject_ratio(TCR_combined_filtered_unsort, "TRBV")
trbv112 <- filter(TCR_combined_filtered_TRBV_subject, Var1 == "TRBV11-2", Class == "MIS-C")
p <- ggplot(trbv112, aes(x=days_since_admission, y=ratio))+
  geom_point(color = "#E7B800", size = 2) + geom_smooth(method = "lm", fill = "#E7B800", color = "#E7B800", alpha=0.3) +
  theme_bw() + stat_cor(method = "pearson") + ggtitle("total T cell, TRBV11-2") +
  xlab("Days_since_hospitalization")
ggsave(p, device = "pdf", filename = "output/plots/TCR_allTcells_TRBV112_time.pdf", width = 4, height = 3)

# plot only TRBV11-2 usage
TCR_combined_filtered_TRBV112_T40 <- filter(cal_Vgene_subject_ratio(TCR_combined_filtered_unsort_T40, "TRBV"), 
                                            Var1 == "TRBV11-2")
p <- plot_Vgene_usage_ratio(TCR_combined_filtered_TRBV112_T40, "all_T_cells", "TRBV11-2")
ggsave(p, device = "pdf", filename = "output/plots/TCR_allTcells_TRBV112_T40.pdf", width = 4, height = 3)



##################################
# cal ratio -- CD4 T cells #######
##################################
# TRBV gene
TCR_combined_filtered_CD4 <- filter(TCR_combined_filtered_unsort, coarsecelltype == "CD4")
TCR_combined_filtered_CD4_T40 <- filter(TCR_combined_filtered_CD4, days_since_admission<41)

TCR_combined_filtered_TRBV_CD4 <- cal_Vgene_subject_ratio(TCR_combined_filtered_CD4_T40, "TRBV")
p <- plot_Vgene_usage_ratio(TCR_combined_filtered_TRBV_CD4, "CD4", "TRBV")
ggsave(p, device = "pdf", filename = "output/plots/TCR_CD4_TRBVgene_T40.pdf", width = 25, height = 17.5)

# time effect Fig.5c
TCR_combined_filtered_TRBV_CD4 <- cal_Vgene_subject_ratio(TCR_combined_filtered_CD4, "TRBV")
cd4_trbv112 <- filter(TCR_combined_filtered_TRBV_CD4, Var1 == "TRBV11-2", Class == "MIS-C")
p <- ggplot(cd4_trbv112, aes(x=days_since_admission, y=ratio))+
  geom_point(color = "#e6b906", size = 2) + geom_smooth(method = "lm", fill = "#e6b906", color = "#e6b906", alpha=0.3) +
  theme_bw() + stat_cor(method = "pearson") + ggtitle("CD4 T cell, TRBV11-2") +
  xlab("Days_since_hospitalization")
ggsave(p, device = "pdf", filename = "output/plots/TCR_CD4Tcells_TRBV112_time.pdf", width = 4, height = 3)

# plot only TRBV11-2 usage for figures
TCR_combined_filtered_TRBV112_CD4_T40 <- filter(cal_Vgene_subject_ratio(TCR_combined_filtered_CD4_T40, "TRBV"), 
                                            Var1 == "TRBV11-2")
p <- plot_Vgene_usage_ratio(TCR_combined_filtered_TRBV112_CD4_T40, "CD4_T_cells", "TRBV11-2")
ggsave(p, device = "pdf", filename = "output/TCR_CD4Tcells_TRBV112_T40.pdf", width = 4, height = 3)


















