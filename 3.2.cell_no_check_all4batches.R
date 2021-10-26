library(Seurat)
library(matrixStats)
library(plyr)
library(tidyverse)
library(openxlsx)
library(reshape2)
library(ggpubr)
library(dplyr)
source("util_funs/cell_freq_check_util_funs.R")

### cell # check of transferred labels #####################################################
### transferred labels = WCT.coursecelltype from batches 1-3 ###
### and combine with adult data ############################################################
### Download Seurat object into input/ folder
merge <- readRDS("input/misc_SeuratObj_submission.rds")
merge_unt <- subset(merge, subset = Batch == "B4chip3", invert = TRUE)
merge_sort <- subset(merge, subset = Batch == "B4chip3")

### read in metadata
samplemetadata = read.csv("input/B4sample.manualhto.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
samplemetadata <- samplemetadata %>%
  filter(Timepoint != "na", Class %in% c("HC","MIS-C","COVID")) %>%
  mutate(sample_id = paste(Subject, Timepoint, sep = "_")) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)
# write.csv(samplemetadata, "B4sample.meta.csv")

### cell ratios unsort
merge_cell_UnSort_WCTcourse <- data.frame(table(merge_unt$sample_id, merge_unt$predicted.id)) %>% 
  dplyr::group_by(Var1) %>% dplyr::mutate(ratio = Freq/sum(Freq))

dir.create("./label_transfer", recursive = TRUE)
saveRDS(merge_cell_UnSort_WCTcourse, "label_transfer/unt_transfWCTcourse_sample_freq.rds")
merge_cell_UnSort_WCTcourse <- readRDS("label_transfer/unt_transfWCTcourse_sample_freq.rds")

merge_cell_UnSort_WCTcourse <- addcellmeta(merge_cell_UnSort_WCTcourse)
merge_cell_UnSort_WCTcourse$Class <- replace(as.character(merge_cell_UnSort_WCTcourse$Class), 
                                             merge_cell_UnSort_WCTcourse$Class == "HC", "ped.HC")
merge_cell_UnSort_WCTcourse$Class <- replace(as.character(merge_cell_UnSort_WCTcourse$Class), 
                                             merge_cell_UnSort_WCTcourse$Class == "COVID", "ped.COVID")
merge_cell_UnSort_WCTcourse$Timepoint <- merge_cell_UnSort_WCTcourse$Timepoint2
merge_cell_UnSort_WCTcourse <- merge_cell_UnSort_WCTcourse %>%
  dplyr::select(Var1, Var2, ratio, Subject, Class, Timepoint, days_since_symptom_onset) %>%
  mutate(group = "Pediatric")

# read in first 3 batches adult data, downloaded from https://github.com/niaid/covid19-time-resolved/tree/main/A1/input
batch1_3_freq_UnSort_mtx <- read.csv("input/ACOVID_final_full_cell_freq_UnSort_mtx.20200818.csv",
                                       header = TRUE, row.names = 1, check.names = FALSE)
colnames(batch1_3_freq_UnSort_mtx)
batch1_3_freq_UnSort_mtx_select <- dplyr::select(batch1_3_freq_UnSort_mtx, c(1:33)) %>%
  mutate(B_Naive = batch1_3_freq_UnSort_mtx$B_Naive_to_parent*batch1_3_freq_UnSort_mtx$CD19) %>%
  mutate(B_Mem = batch1_3_freq_UnSort_mtx$B_Mem_to_parent*batch1_3_freq_UnSort_mtx$CD19) %>%
  mutate(PB_Plasmablasts = batch1_3_freq_UnSort_mtx$PB_Plasmablasts_to_parent*batch1_3_freq_UnSort_mtx$CD19) %>%
  mutate(CD8_Naive = batch1_3_freq_UnSort_mtx$CD8_Naive_to_parent*batch1_3_freq_UnSort_mtx$CD8) %>%
  mutate(CD8_Mem = batch1_3_freq_UnSort_mtx$CD8_Mem_to_parent*batch1_3_freq_UnSort_mtx$CD8) %>%
  mutate(CD4_Naive = batch1_3_freq_UnSort_mtx$CD4_Naive_to_parent*batch1_3_freq_UnSort_mtx$CD4) %>%
  mutate(CD4_Mem = batch1_3_freq_UnSort_mtx$CD4_Mem_to_parent*batch1_3_freq_UnSort_mtx$CD4)

batch1_3_freq_UnSort_mtx_select <- reshape2::melt(batch1_3_freq_UnSort_mtx_select, 
                                        id.vars = c("Var1", "Timepoint", "sample_id",
                                                    "Batch", "severity", "days_since_symptoms_onset",
                                                    "PC1", "PC1class", "Subject", "severity_outcome"),
                                        variable.name = "Var2",
                                        value.name = "ratio") %>%
  dplyr::select(Var1, Var2, ratio, Subject, PC1class, Timepoint, days_since_symptoms_onset) %>%
  filter(!is.na(PC1class)) %>%
  mutate(group = "Adult")

colnames(batch1_3_freq_UnSort_mtx_select) <- colnames(merge_cell_UnSort_WCTcourse)

batch1_4_freq <- rbind(batch1_3_freq_UnSort_mtx_select, merge_cell_UnSort_WCTcourse)
colnames(batch1_4_freq) <- c("sample_id","variable","value","Subject","Class","Timepoint","days_since_symptom_onset","group")
batch1_4_freq$Class2 <- batch1_4_freq$Class
# batch1_4_freq$Class <- replace(batch1_4_freq$Class, batch1_4_freq$Class %in% c("PC1_high","PC1_low"), "adult.COVID")
# batch1_4_freq$Class <- factor(batch1_4_freq$Class, levels = c("HC","ped.HC","adult.COVID","ped.COVID","MIS-C"))
batch1_4_freq$Class <- factor(batch1_4_freq$Class, levels = c("ped.HC", "ped.COVID","MIS-C","HC","PC1_low","PC1_high"))

batch1_4_freq$days_since_symptom_onset <- as.numeric(as.character(batch1_4_freq$days_since_symptom_onset))


### use only the first timepoint for each subject (T0)
### Figure S5B
df_T0 <- filter(batch1_4_freq, Timepoint %in% c("T0","HC")) %>%
  filter(variable %in% c("Mono_Classical","Mono_NonClassical","NK_CD16hi","gammadeltaT",
                         "MAIT","pDC","cDC","TissuResMemT", "Platelets", "B_Naive",
                         "B_Mem", "PB_Plasmablast","CD8_Naive","CD8_Mem", "CD4_Naive","CD4_Mem"))
my_comparisons <- list(c("ped.HC", "ped.COVID"), c("ped.COVID", "MIS-C"), c("PC1_low","PC1_high"), 
                       c("ped.HC","MIS-C"), c("PC1_low","ped.COVID"))
p1 <- ggplot(df_T0, aes(x = Class, y = value, color = Class))+
  geom_boxplot(outlier.shape=NA)+geom_point(size = 2)+
  # ggsci::scale_fill_npg(name="Class",alpha = 0.5) +
  scale_color_manual(values = c("HC" = "#6664c0", "PC1_low" = "#73bde0", "PC1_high" = "#e1bd76",
                                "ped.HC" = "#0000AC", "ped.COVID" = "#00AFBB", "MIS-C" = "#E7B800"))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", )+
  facet_wrap(~variable, scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30,hjust=1), text = element_text(size=18))
p1$layers[[2]]$aes_params$textsize <- 9
ggsave("output/batch1_4_WCTcourse.class.tototal.T0.pdf", p1, device = "pdf", width = 18, height = 18)




