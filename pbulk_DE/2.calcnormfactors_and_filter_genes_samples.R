library(limma)
library(edgeR)

###############################################################################
### caculate norm factors, normalize pseudobulk libraries and filter genes ####
###############################################################################
source("pbulk_DE/util/limma_DE_functions.R")

DGELIST_IN_PATH <- "output/pbulk_DE/all_samples/pseudobulk_dgelists_unfiltered/"

dir.create("output/pbulk_DE/all_samples/pseudobulk_dgelists_normalized/", recursive = TRUE)
DGELIST_OUT_PATH <- "output/pbulk_DE/all_samples/pseudobulk_dgelists_normalized/"

dge_list <- readRDS(paste0(DGELIST_IN_PATH, "UnSort-mergedcelltype.rds"))

dge_list <- lapply(dge_list, function(dge){NormalizePseudobulk(dge)})

saveRDS(dge_list, paste0(DGELIST_OUT_PATH, "UnSort-mergedcelltype_pseudobulk_dgelist_normalized.rds"))


##############################################
### subset pseudobulk Timepoint ##############
##############################################
PBULKLIST_IN_PATH = "output/pbulk_DE/all_samples/pseudobulk_dgelists_normalized/"
PBULKLIST_OUT_PATH = "output/pbulk_DE/sample_groups/"
dir.create("output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_healthy/pseudobulk_dgelists_normalized/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/tso_within_d40_covid_plus_healthy/pseudobulk_dgelists_normalized/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_covid/pseudobulk_dgelists_normalized/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/all_timepoints_misc_only/pseudobulk_dgelists_normalized/", recursive = TRUE)
dir.create("output/pbulk_DE/sample_groups/all_timepoints_covid_only/pseudobulk_dgelists_normalized/", recursive = TRUE)

# ---1
pbulklist <- readRDS(paste0(PBULKLIST_IN_PATH, "UnSort-mergedcelltype_pseudobulk_dgelist_normalized.rds"))

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
saveRDS(pbulklist_d40_misc_plus_healthy, paste0(PBULKLIST_OUT_PATH, "tso_within_d40_misc_plus_healthy/pseudobulk_dgelists_normalized/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_d40_covid_plus_healthy, paste0(PBULKLIST_OUT_PATH, "tso_within_d40_covid_plus_healthy/pseudobulk_dgelists_normalized/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_d40_misc_plus_covid, paste0(PBULKLIST_OUT_PATH, "tso_within_d40_misc_plus_covid/pseudobulk_dgelists_normalized/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_all_misc_only, paste0(PBULKLIST_OUT_PATH, "all_timepoints_misc_only/pseudobulk_dgelists_normalized/UnSort-mergedcelltype.rds"))
saveRDS(pbulklist_all_covid_only, paste0(PBULKLIST_OUT_PATH, "all_timepoints_covid_only/pseudobulk_dgelists_normalized/UnSort-mergedcelltype.rds"))






