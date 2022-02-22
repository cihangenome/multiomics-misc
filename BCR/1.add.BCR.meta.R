library(Seurat) #load Seurat 3.1
library(matrixStats)
library(plyr)
library(tidyverse)

### add BCR meta to Seurat object
### input Seurat Object is the output from 3.1.label.transfer.from.3batches.R
# the BCR data is already included in the Seurat object, this is just an illustration of the process
merge <- readRDS("input/misc_SeuratObj_submission.rds")
Idents(merge) <- "sample"
samplemetadata = read.csv("input/B4sample.manualhto.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(Class %in% c("HC","COVID","MIS-C"), Pool != "misassigned")
identical(sort(samplemetadata$x), sort(unique(merge$sample)))
merge$Timepoint2 <- plyr::mapvalues(x = Idents(merge), from = samplemetadata$x, to = samplemetadata$Timepoint2)


bcr_db <- read.csv("input/immcant_bcr_ighl_filtered.csv", header = TRUE, row.names = 1)
bcr_db_short <- bcr_db %>% dplyr::select(BARCODES,
                                  V_CALL, D_CALL, J_CALL, V_CALL_10X, D_CALL_10X, J_CALL_10X,
                                  V_CALL_10X_light, J_CALL_10X_light,
                                  GERMLINE_V_CALL, GERMLINE_D_CALL, GERMLINE_J_CALL,
                                  LOCUS_light, clone_id_hierac_thrsld_0.35, clonefamily,
                                  mu_freq_tot, mu_freq_tot_light, mu_freq_HL)
meta <- merge@meta.data %>% left_join(bcr_db_short, by = c("NewBarcode"="BARCODES"))
rownames(meta) <- meta$NewBarcode
merge <- AddMetaData(merge, meta[,c(colnames(bcr_db_short))])

saveRDS(merge, "input/B4merge.Clustered.COVID.QCd.wBCRmeta.rds")




