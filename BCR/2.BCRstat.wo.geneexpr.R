library(Seurat)
library(tidyverse)
library(reshape2)
library(ggsci)
library(ggfortify)
library(ggpubr)
source("BCR/util/BCR_usage_functions.R")

#######################################################
### MISC BCR analysis #################################
#######################################################
# Generates Fig.5d, Extended Data Fig.5b, 5d
# This file need the BCR data: immcant_bcr_ighl_filtered.csv downloaded in input folder 
# and Seurat object downloaded from Zenodo
BCR_imm_combined <- read.csv("input/immcant_bcr_ighl_filtered.csv", header = TRUE, row.names = 1)
BCR_imm_combined_filtered <- BCR_imm_combined %>% filter(COARSECELLTYPE %in% c("B","Plasmablast")) %>% 
  mutate(CLASS = factor(CLASS, levels = c("HC","COVID","MIS-C")))

# read in metadata
meta_subject <- read.csv("input/B4sample.manualhto.csv", row.names = 1) %>%
  filter(Class %in% c("HC","COVID","MIS-C"), Pool != "misassigned")
meta_subject$sample_id <- paste(meta_subject$Subject, meta_subject$Timepoint, sep = "_")
meta_subject <- meta_subject[!duplicated(meta_subject$sample_id),]
meta_subject$Class <- factor(meta_subject$Class, levels = c("HC","COVID","MIS-C"))
merge <- readRDS("input/misc_SeuratObj_submission.rds")
meta <- merge@meta.data
meta$Class <- factor(meta$Class, levels = c("HC","COVID","MIS-C"))
meta_short <- dplyr::select(meta, NewBarcode, Sorted)

BCR_imm_combined_filtered_unsort <- BCR_imm_combined_filtered %>% left_join(meta_short, by = c("BARCODES" = "NewBarcode")) %>%
  filter(Sorted == "UnSort")

#######################################################
### V gene usage ######################################
#######################################################
# Extended Data Fig.5b
# heavy chain
BCR_imm_combined_subject_hc <- data.frame(table(BCR_imm_combined_filtered_unsort$V_CALL_10X, BCR_imm_combined_filtered_unsort$SAMPLE_ID)) %>% 
  dplyr::group_by(Var2) %>%
  dplyr::mutate(ratio = Freq/sum(Freq)) %>%
  drop_na() %>%
  # dcast(Var2 ~ Var1, value.var = "ratio") %>%
  left_join(meta_subject, by = c("Var2" = "sample_id")) %>%
  filter(days_since_admission < 41)

# plot only IGHV4-34 of days_since_admission within d40
BCR_imm_combined_subject_hc_selected <- filter(BCR_imm_combined_subject_hc, Var1 %in% c("IGHV4-34"))
p <- plot_Vgene_usage_ratio(BCR_imm_combined_subject_hc_selected, "total unsorted B", "IGHV4-34")
p
# ggsave(p, device = "pdf", filename = "output/plots/BCR_allBcells_VHgene_selected_T40.pdf", width = 3.5, height = 3.2)

#######################################################
### Mutation Frequency ################################
#######################################################
# showing all IgH sequenced
db_filtered_IgH <- readRDS("input/imgt_mut_load_igh_cellfiltered.rds")
db_filtered_IgH$CLASS <- factor(db_filtered_IgH$CLASS, levels = c("HC","COVID","MIS-C"))
db_filtered_PB <- filter(db_filtered_IgH, WCT.MERGEDCELLTYPE %in% c("Plasmablast"))
my_comparisons <- list( c("HC", "COVID"), c("COVID", "MIS-C"), c("HC", "MIS-C"))

# Fig.5d
p <- ggplot(db_filtered_PB, aes(x=CLASS, y=mu_freq_tot, color=CLASS)) +
  theme_bw() + ggtitle("Mutation quantification_Plasmablast") +
  scale_color_manual(values = c("HC" = "#0800ac", "COVID" = "#1cb0bb", "MIS-C" = "#e6ba07"))+
  ylab("Mutation frequency") +
  geom_boxplot()+
  geom_point(position = "jitter", size = 1)+
  stat_compare_means(comparisons = my_comparisons, mmethod = "wilcox.test")
p
# ggsave("output/plots/igh_mut_PB.pdf",p, device = "pdf", width = 4.5, height = 3.7)

# Extended Data Fig.5d
db_filtered_Mem <- filter(db_filtered_IgH, MERGEDCELLTYPE %in% c("B_Mem"))
p <- ggplot(db_filtered_Mem, aes(x=CLASS, y=mu_freq_tot, color=CLASS)) +
  theme_bw() + ggtitle("Mutation quantification_Mem B cells") +
  scale_color_manual(values = c("HC" = "#0800ac", "COVID" = "#1cb0bb", "MIS-C" = "#e6ba07"))+
  ylab("Mutation frequency") +
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, mmethod = "wilcox.test")
p
# ggsave("output/plots/igh_mut_MemB.pdf",p, device = "pdf", width = 4.5, height = 3.7)

         



