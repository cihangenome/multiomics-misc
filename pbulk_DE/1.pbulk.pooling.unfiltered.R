library(Seurat)
library(matrixStats)
library(limma)
library(reshape2)
library(edgeR)

#############################################
### pseudobulk pooling -- mergedcelltype ####
#############################################

source("pbulk_DE/util/pbulk_pooling_functions.R")
SEURAT_IN_PATH <- "input/"
DGELISTS_OUT_PATH <- "output/pbulk_DE/all_samples/pseudobulk_dgelists_unfiltered/"
dir.create(DGELISTS_OUT_PATH, recursive = TRUE)

# Sample filtering parameters
KEEP_SORTED <- c("UnSort", "Sort")

# Parameters for pseudobulk pooling
LIBSIZE_FILTER <- 30000
MIN_CELLS_PER_POOL <- 4
MIN_SAMPLES_PER_CELLTYPE <- 5

# Cell annotation column
CELL_ANNOTATION_COLUMN <- "mergedcelltype"

# read in data
merged <- readRDS(SEURAT_IN_PATH)

# filter to just covid samples
keep_cells_UnSort <- WhichCells(merged, expression = Sorted == KEEP_SORTED[1])
keep_cells_Sort <- WhichCells(merged, expression = Sorted == KEEP_SORTED[2])

seurat_obj <- subset(merged, cells = keep_cells_UnSort)

meta <- seurat_obj@meta.data

rna <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
# rm(seurat_obj)

sample_cols <- c("Age", "Gender", "Status", 
                 "Timepoint", "Subject", "Class", "days_since_admission",
                 "sample_id", "Timepoint2")

sample_cols <- c(sample_cols, CELL_ANNOTATION_COLUMN)

print("sample_cols that aren't in metadata")
sample_cols[!sample_cols %in% colnames(meta)]

samples <- paste(meta$Subject, meta$Timepoint, sep = "_")

pseudobulk_list <- getPseudobulkList(mat = rna, celltypes = meta[[CELL_ANNOTATION_COLUMN]], 
                                     meta = meta,
                                     samples = samples, 
                                     min_cells_per_pool = MIN_CELLS_PER_POOL, 
                                     min_samples_per_celltype = MIN_SAMPLES_PER_CELLTYPE,
                                     barcode_col_name = "NewBarcode", 
                                     sample_level_meta_cols = sample_cols,
                                     #cell_level_meta_cols = cell_cols, 
                                     pooling_function = "sum",
                                     output_type = "DGEList")


pseudobulk_list_libfiltered <- lapply(pseudobulk_list, function(dge){
  dge[, dge$samples$lib.size > LIBSIZE_FILTER]
  
})


n_samples_per_celltype <- sapply(pseudobulk_list, ncol)
print("n_samples_per_celltype before final filtering")
print(n_samples_per_celltype)

pseudobulk_list_libfiltered <- pseudobulk_list_libfiltered[n_samples_per_celltype > MIN_SAMPLES_PER_CELLTYPE]

n_samples_per_celltype <- sapply(pseudobulk_list_libfiltered, ncol)
print("n_samples_per_celltype after final filtering")
print(n_samples_per_celltype)

saveRDS(pseudobulk_list_libfiltered, paste0(DGELISTS_OUT_PATH, "UnSort-mergedcelltype.rds"))


##############################################
### subset pseudobulk Timepoint ##############
##############################################
PBULKLIST_IN_PATH = "output/pbulk_DE/all_samples/pseudobulk_dgelists_unfiltered/"
PBULKLIST_OUT_PATH = "output/pbulk_DE/sample_groups/"
dir.create("output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_healthy/pseudobulk_dgelists_unfiltered/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/tso_within_d40_covid_plus_healthy/pseudobulk_dgelists_unfiltered/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_covid/pseudobulk_dgelists_unfiltered/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/all_timepoints_misc_only/pseudobulk_dgelists_unfiltered/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/all_timepoints_covid_only/pseudobulk_dgelists_unfiltered/", recursive = TRUE)

# ---1
pbulklist <- readRDS(paste0(PBULKLIST_IN_PATH, "UnSort-mergedcelltype.rds"))

pbulklist_d40_misc_plus_healthy <- lapply(pbulklist, function(dge){
  dge[, dge$samples$days_since_admission < 41 & dge$samples$Class %in% c("MIS-C","HC")]
})

pbulklist_d40_covid_plus_healthy <- lapply(pbulklist, function(dge){
  dge[, dge$samples$days_since_admission < 41 & dge$samples$Class %in% c("COVID","HC")]
})

pbulklist_d40_misc_plus_covid <- lapply(pbulklist, function(dge){
  dge[, dge$samples$days_since_admission < 41 & dge$samples$Class %in% c("MIS-C","COVID")]
})

pbulklist_all_misc_only <- lapply(pbulklist, function(dge){
  dge[, dge$samples$Class %in% c("MIS-C")]
})

pbulklist_all_covid_only <- lapply(pbulklist, function(dge){
  dge[, dge$samples$Class %in% c("COVID")]
})


# ---1
saveRDS(pbulklist_d40_misc_plus_healthy, paste0(PBULKLIST_OUT_PATH, "tso_within_d40_misc_plus_healthy/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_d40_covid_plus_healthy, paste0(PBULKLIST_OUT_PATH, "tso_within_d40_covid_plus_healthy/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_d40_misc_plus_covid, paste0(PBULKLIST_OUT_PATH, "tso_within_d40_misc_plus_covid/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_all_misc_only, paste0(PBULKLIST_OUT_PATH, "all_timepoints_misc_only/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_all_covid_only, paste0(PBULKLIST_OUT_PATH, "all_timepoints_covid_only/pseudobulk_dgelists_unfiltered/UnSort-mergedcelltype.rds"))

