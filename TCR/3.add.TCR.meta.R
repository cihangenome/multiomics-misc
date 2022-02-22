library(Seurat) #load Seurat 3.1
library(matrixStats)
library(plyr)
library(tidyverse)

######################################
### add TCR meta to Seurat object ####
######################################

### Seurat object from output of BCR/1.add.BCR.meta.R
# the TCR data is already included in the Seurat object, this is just an illustration of the process
merge <- readRDS("input/B4merge.Clustered.COVID.QCd.wBCRmeta.rds")

# read in TCR data
# used tenx output TCR data from script 1.compile.TCR.data.R
TCR_combined <- readRDS("output/TCR_tenx_filtered_anno_wofiltered_celltype_combined.rds")
TCR_combined_filtered <- filter(TCR_combined, !is.na(Class), coarsecelltype %in% c("CD4","CD8","MAIT","T_gd","T_Vd2"))
TCR_combined_filtered_abboth <- filter(TCR_combined_filtered, !is.na(TCR1), !is.na(TCR2))
TCR_combined_filtered$TCR_both_chain_exist <- TCR_combined_filtered$tenx_barcode %in% TCR_combined_filtered_abboth$tenx_barcode
TCR_freq <- sort(table(TCR_combined_filtered_abboth$CTgene_tcr), decreasing = TRUE)
# set clone id from frequency high to low
clone_id <- seq(1, length(TCR_freq))
names(clone_id) <- names(TCR_freq)
TCR_combined_filtered$TCR_clone_id <- as.numeric(plyr::mapvalues(x = TCR_combined_filtered$CTgene_tcr, from = names(clone_id), to = clone_id))

TCR_combined_filtered_short <- select(TCR_combined_filtered, 
                                      tenx_barcode, TCR1, cdr3_aa1_tcr, cdr3_nt1_tcr, 
                                      TCR2, cdr3_aa2_tcr, cdr3_nt2_tcr, CTgene_tcr, 
                                      TRAV, TRAJ, TRBV, TRBJ, TRBC,
                                      TCR_both_chain_exist, TCR_clone_id)
meta <- merge@meta.data %>% left_join(TCR_combined_filtered_short, by = c("NewBarcode"="tenx_barcode"))
rownames(meta) <- meta$NewBarcode
merge <- AddMetaData(merge, meta[,c(rownames(TCR_combined_filtered_short))])
merge@meta.data <- merge@meta.data[,-which(str_detect(colnames(merge@meta.data),"prediction.score"))[-31]]

# this is the seurat object with both TCR and BCR meta -> misc_SeuratObj_submission.rds
saveRDS(merge, "input/B4merge.Clustered.COVID.QCd.wBCRTCRmeta.rds")



