library(Seurat)
library(tidyverse)
library(scRepertoire)
library(reshape2)

############################################################
### MISC compile scTCR data with Seurat scMeta data ########
############################################################
# This file need the 
# Seurat object: misc_SeuratObj_submission.rds
# TCR raw data: TCR_filtered_contig_annotations.csv 
# downloaded in input folder
# But the TCR data is already included in the Seurat object, this is just an illustration of the process
SEURAT_META_IN_PATH <- "input/misc_SeuratObj_submission.rds"
B4_TENX_IN_PATH <- "input/TCR_filtered_contig_annotations.csv"
OUT_PATH <- "output/TCR_tenx_filtered_anno_wofiltered_celltype.rds"

merge <- readRDS(SEURAT_META_IN_PATH)
b4_tenx <- read.csv(B4_TENX_IN_PATH, stringsAsFactors = FALSE)

barcode_pre <- substr(b4_tenx$barcode, 1, 16)
barcode_suffix <- sapply(strsplit(b4_tenx$barcode, "-"), `[[`, 2)
for(i in seq_along(barcode_suffix)){
  if(nchar(barcode_suffix[[i]]) == 1){
    barcode_suffix[[i]] <- paste0("0", barcode_suffix[[i]])
  }
}
barcodes <- paste(barcode_pre, barcode_suffix, sep = "")

b4_tenx$tenx_barcode <- barcodes

meta <- merge@meta.data
meta$Class <- factor(meta$Class, levels = c("HC","COVID","MIS-C"))
out_dat <- left_join(b4_tenx, meta, by = c("tenx_barcode" = "NewBarcode"))
saveRDS(out_dat, OUT_PATH)

####################################################################################
### format TCR filtered contig-file to have combined a and b chains for one cell ###
####################################################################################

b4_tenx <- read.csv("input/TCR_filtered_contig_annotations.csv", header = TRUE)
b4_tenx_filtered <- filter(b4_tenx, 
                           is_cell == "True", 
                           high_confidence == "True", 
                           productive == "True")

tcr_combined_filtered <- scRepertoire::combineTCR(b4_tenx_filtered, 
                                                  samples = "batch4", 
                                                  cells ="T-AB",
                                                  ID = "batch4",
                                                  filterMulti = TRUE)
tcr_combined_filtered <- tcr_combined_filtered$batch4_batch4

barcode_pre <- substr(tcr_combined_filtered$barcode, 15, 30)
barcode_suffix <- sapply(strsplit(tcr_combined_filtered$barcode, "-"), `[[`, 2)
for(i in seq_along(barcode_suffix)){
  if(nchar(barcode_suffix[[i]]) == 1){
    barcode_suffix[[i]] <- paste0("0", barcode_suffix[[i]])
  }
}
barcodes <- paste(barcode_pre, barcode_suffix, sep = "")

tcr_combined_filtered$tenx_barcode <- barcodes
# split TCRa and b genes
tcr_combined_filtered <- tcr_combined_filtered %>% mutate(TCRa = tcr_combined_filtered$TCR1) %>% 
  separate(TCRa, c("TRAV", "TRAJ", NA), sep = "\\.")
tcr_combined_filtered <- tcr_combined_filtered %>% mutate(TCRb = tcr_combined_filtered$TCR2) %>% 
  separate(TCRb, c("TRBV", "TRBJ", "TRBC"), sep = "\\.")
# rename columns to differentiation BCR/TCR
colnames(tcr_combined_filtered) <- c("barcode", "sample", "ID", "TCR1", "cdr3_aa1_tcr", "cdr3_nt1_tcr",
                                     "TCR2", "cdr3_aa2_tcr", "cdr3_nt2_tcr", "CTgene_tcr", "CTnt_tcr",
                                     "CTaa_tcr", "CTstrict_tcr", "cellType_tcr", "tenx_barcode", 
                                     "TRAV", "TRAJ", "TRBV", "TRBJ", "TRBC")

meta <- merge@meta.data[,c(1:70, 84:87)] # didn't select the TCR columns in the scMeta
meta$Class <- factor(meta$Class, levels = c("HC","COVID","MIS-C"))
out_dat_combined <- left_join(tcr_combined_filtered, meta, by = c("tenx_barcode" = "NewBarcode"))

out_dat_combined$cdr3_length1 <- nchar(out_dat_combined$cdr3_aa1)
out_dat_combined$cdr3_length2 <- nchar(out_dat_combined$cdr3_aa2)

saveRDS(out_dat_combined, "output/TCR_tenx_filtered_anno_wofiltered_celltype_combined.rds")






