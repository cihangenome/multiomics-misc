library(limma)
library(edgeR)

# run limma differential expression model
#############################
### run DE ##################
#############################
source("pbulk_DE/util/limma_DE_functions.R")

DGELIST_IN_PATH <- "output/pbulk_DE/sample_groups/"
SAMPLE_GROUP <- c("tso_within_d40_misc_plus_healthy/", 
                  "tso_within_d40_covid_plus_healthy/", 
                  "tso_within_d40_misc_plus_covid/",
                  "all_timepoints_misc_only/", 
                  "all_timepoints_covid_only/")
PATIENT_ID_COL <- "Subject"

#########################
### MIS-C vs HC d40 #####
#########################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[1], "pseudobulk_dgelists_normalized/", "UnSort-mergedcelltype.rds"))
design_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_design.rds"))

dge <- dge_Unsort
design <- design_Unsort

celltypes = names(dge)
fit.list = list()

###
for (i in celltypes) {
  
#model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the DGEList
dge[[i]] <- dge[[i]][, rownames(design[[i]])]

fit <- RunVoomLimma(dge[[i]], design_matrix = design[[i]], add_random_effect = TRUE, 
                    block_this_variable = PATIENT_ID_COL,
                    do_contrast_fit = FALSE, my_contrast_matrix = NULL)

fit.list[[i]] <- fit

}
###

saveRDS(fit.list, "output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_healthy/UnSort_misc_vs_healthy_fit.rds")



#########################
### COVID vs HC d40 #####
#########################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[2], "pseudobulk_dgelists_normalized/", "UnSort-mergedcelltype.rds"))
design_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_design.rds"))

dge <- dge_Unsort
design <- design_Unsort

celltypes = names(dge)
fit.list = list()

###
for (i in celltypes) {
  
  #model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the DGEList
  dge[[i]] <- dge[[i]][, rownames(design[[i]])]
  
  fit <- RunVoomLimma(dge[[i]], design_matrix = design[[i]], add_random_effect = TRUE, 
                      block_this_variable = PATIENT_ID_COL,
                      do_contrast_fit = FALSE, my_contrast_matrix = NULL)
  
  fit.list[[i]] <- fit
  
}
###

saveRDS(fit.list, "output/pbulk_DE/sample_groups/tso_within_d40_covid_plus_healthy/UnSort_covid_vs_healthy_fit.rds")



############################
### MIS-C vs COVID d40 #####
############################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[3], "pseudobulk_dgelists_normalized/", "UnSort-mergedcelltype.rds"))
design_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_design.rds"))

dge <- dge_Unsort
design <- design_Unsort

celltypes = names(dge)
fit.list = list()

###
for (i in celltypes) {
  
  #model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the DGEList
  dge[[i]] <- dge[[i]][, rownames(design[[i]])]
  
  fit <- RunVoomLimma(dge[[i]], design_matrix = design[[i]], add_random_effect = TRUE, 
                      block_this_variable = PATIENT_ID_COL,
                      do_contrast_fit = FALSE, my_contrast_matrix = NULL)
  
  fit.list[[i]] <- fit
  
}
###

saveRDS(fit.list, "output/pbulk_DE/sample_groups/tso_within_d40_misc_plus_covid/UnSort_misc_vs_covid_fit.rds")




#######################################
### MIS-C only all timepoints #########
#######################################

dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[4], "pseudobulk_dgelists_normalized/", "UnSort-mergedcelltype.rds"))
design_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_design.rds"))

dge <- dge_Unsort
design <- design_Unsort

# filter out celltypes with <=3 subjects
# remove B_Naive -- error correlation is 1 or -1, so the model is degenerate
celltypes = names(design)[sapply(design, nrow) > 3 & names(design) != "B_Naive"]
fit.list = list()

###
for (i in celltypes) {
  
  #model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the DGEList
  dge[[i]] <- dge[[i]][, rownames(design[[i]])]
  
  fit <- RunVoomLimma(dge[[i]], design_matrix = design[[i]], add_random_effect = TRUE, 
                      block_this_variable = PATIENT_ID_COL,
                      do_contrast_fit = FALSE, my_contrast_matrix = NULL)
  
  fit.list[[i]] <- fit
  
}
###

saveRDS(fit.list, "output/pbulk_DE/sample_groups/all_timepoints_misc_only/UnSort_days_onset_fit.rds")


#######################################
### COVID only all timepoints #########
#######################################

dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[5], "pseudobulk_dgelists_normalized/", "UnSort-mergedcelltype.rds"))
design_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_design.rds"))

dge <- dge_Unsort
design <- design_Unsort

# filter out celltypes with <3 subjects
celltypes = names(design)[sapply(design, nrow) > 3]

fit.list = list()

###
for (i in celltypes) {
  
  #model.matrix will remove rows with NA's for any of the columns, so I subset to get rid of those in the DGEList
  dge[[i]] <- dge[[i]][, rownames(design[[i]])]
  
  fit <- RunVoomLimma(dge[[i]], design_matrix = design[[i]], add_random_effect = TRUE, 
                      block_this_variable = PATIENT_ID_COL,
                      do_contrast_fit = FALSE, my_contrast_matrix = NULL)
  
  fit.list[[i]] <- fit
  
}
###

saveRDS(fit.list, "output/pbulk_DE/sample_groups/all_timepoints_covid_only/UnSort_days_onset_fit.rds")

