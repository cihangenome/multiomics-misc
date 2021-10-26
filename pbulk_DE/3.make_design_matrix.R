library(limma)
library(edgeR)

# Make design matrix for each celltype 
#############################
### make design matrices ####
#############################
source("pbulk_DE/util/limma_DE_functions.R")
source("pbulk_DE/util/fix_single_level_contrasts.R")

DGELIST_IN_PATH <- "pbulk_DE/sample_groups/"
SAMPLE_GROUP <- c("tso_within_d40_misc_plus_healthy/pseudobulk_dgelists_normalized/", 
                  "tso_within_d40_covid_plus_healthy/pseudobulk_dgelists_normalized/",
                  "tso_within_d40_misc_plus_covid/pseudobulk_dgelists_normalized/", 
                  "all_timepoints_misc_only/pseudobulk_dgelists_normalized/", 
                  "all_timepoints_covid_only/pseudobulk_dgelists_normalized/")

# Make design matrix ----------------------------------------------
#########################
### MIS-C vs HC d40 #####
#########################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[1], "UnSort-mergedcelltype.rds"))

FORMULA <- "~0 + Class + Age"
FORMULA <- as.formula(FORMULA)

dge <- dge_Unsort
dge <- dge_All

celltypes = names(dge)
design.matrix.list = list()

###
for (i in celltypes) {

meta <- dge[[i]]$samples
# print(colnames(meta))
# formulate the age to be numeric
meta$Age <- as.numeric(as.character(meta$Age))
meta$Class <- as.factor(replace(as.character(meta$Class), meta$Class == "MIS-C", "MISC"))

# drop unused factor levels
for(nm in colnames(meta)){
  if(is.factor(meta[[nm]])){
    meta[[nm]] <- droplevels(meta[[nm]])
  }
}

FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)

design <- model.matrix(FORMULA, meta)

colnames(design) <- gsub(":", "X", colnames(design))

if (identical(rownames(design), colnames(dge[[i]]$counts)) == FALSE){
  stop(paste0(celltype, ": subjects(rows of meta) don't match with count columns"))
}

design.matrix.list[[i]] = design

}
###

saveRDS(design.matrix.list, "pbulk_DE/sample_groups/tso_within_d40_misc_plus_healthy/results/UnSort_misc_vs_healthy_design.rds")



#########################
### COVID vs HC d40 #####
#########################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[2], "UnSort-mergedcelltype.rds"))

FORMULA <- "~0 + Class + Age"
FORMULA <- as.formula(FORMULA)

dge <- dge_Unsort
dge <- dge_All

celltypes = names(dge)
design.matrix.list = list()

###
for (i in celltypes) {
  
  meta <- dge[[i]]$samples
  # print(colnames(meta))
  # formulate the age to be numeric
  meta$Age <- as.numeric(as.character(meta$Age))
  
  # drop unused factor levels
  for(nm in colnames(meta)){
    if(is.factor(meta[[nm]])){
      meta[[nm]] <- droplevels(meta[[nm]])
    }
  }
  
  FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)
  
  design <- model.matrix(FORMULA, meta)
  
  colnames(design) <- gsub(":", "X", colnames(design))
  
  if (identical(rownames(design), colnames(dge[[i]]$counts)) == FALSE){
    stop(paste0(celltype, ": subjects(rows of meta) don't match with count columns"))
  }
  
  design.matrix.list[[i]] = design
  
}
###

saveRDS(design.matrix.list, "pbulk_DE/sample_groups/tso_within_d40_covid_plus_healthy/results/UnSort_covid_vs_healthy_design.rds")




############################
### MIS-C vs COVID d40 #####
############################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[3], "UnSort-mergedcelltype.rds"))

# taking TSO into account
FORMULA <- "~0 + Class + days_since_admission + Age"
FORMULA <- as.formula(FORMULA)

dge <- dge_Unsort
dge <- dge_All

celltypes = names(dge)
design.matrix.list = list()

###
for (i in celltypes) {
  
  meta <- dge[[i]]$samples
  # print(colnames(meta))
  # formulate the age to be numeric
  meta$Age <- as.numeric(as.character(meta$Age))
  meta$Class <- as.factor(replace(as.character(meta$Class), meta$Class == "MIS-C", "MISC"))
  
  # drop unused factor levels
  for(nm in colnames(meta)){
    if(is.factor(meta[[nm]])){
      meta[[nm]] <- droplevels(meta[[nm]])
    }
  }
  
  FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)
  
  design <- model.matrix(FORMULA, meta)
  
  colnames(design) <- gsub(":", "X", colnames(design))
  
  if (identical(rownames(design), colnames(dge[[i]]$counts)) == FALSE){
    stop(paste0(celltype, ": subjects(rows of meta) don't match with count columns"))
  }
  
  design.matrix.list[[i]] = design
  
}
###

saveRDS(design.matrix.list, "pbulk_DE/sample_groups/tso_within_d40_misc_plus_covid/results/UnSort_misc_vs_covid_design.rds")




#######################################
### MIS-C only all timepoints #########
#######################################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[4], "UnSort-mergedcelltype.rds"))

FORMULA <- "~days_since_admission + Age"
FORMULA <- as.formula(FORMULA)

dge <- dge_Unsort
dge <- dge_All

celltypes = names(dge)
design.matrix.list = list()

###
for (i in celltypes) {
  
  meta <- dge[[i]]$samples
  # print(colnames(meta))
  # formulate the age to be numeric
  meta$Age <- as.numeric(as.character(meta$Age))
  meta$Class <- as.factor(replace(as.character(meta$Class), meta$Class == "MIS-C", "MISC"))
  
  # drop unused factor levels
  for(nm in colnames(meta)){
    if(is.factor(meta[[nm]])){
      meta[[nm]] <- droplevels(meta[[nm]])
    }
  }
  
  FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)
  
  design <- model.matrix(FORMULA, meta)
  
  colnames(design) <- gsub(":", "X", colnames(design))
  
  if (identical(rownames(design), colnames(dge[[i]]$counts)) == FALSE){
    stop(paste0(celltype, ": subjects(rows of meta) don't match with count columns"))
  }
  
  design.matrix.list[[i]] = design
  
}
###

saveRDS(design.matrix.list, "pbulk_DE/sample_groups/all_timepoints_misc_only/results/UnSort_days_onset_design.rds")


#######################################
### COVID only all timepoints #########
#######################################
dge_Unsort <- readRDS(paste0(DGELIST_IN_PATH, SAMPLE_GROUP[5], "UnSort-mergedcelltype.rds"))

FORMULA <- "~days_since_admission + Age"
FORMULA <- as.formula(FORMULA)

dge <- dge_Unsort
dge <- dge_All

# remove HSC since only 1 or 2 samples available
celltypes = names(dge)[-which(names(dge) == "HSC")]
design.matrix.list = list()

###
for (i in celltypes) {
  
  meta <- dge[[i]]$samples
  # print(colnames(meta))
  # formulate the age to be numeric
  meta$Age <- as.numeric(as.character(meta$Age))
  
  # drop unused factor levels
  for(nm in colnames(meta)){
    if(is.factor(meta[[nm]])){
      meta[[nm]] <- droplevels(meta[[nm]])
    }
  }
  
  FORMULA <- update_formula_remove_single_level_contrasts(FORMULA, meta)
  
  design <- model.matrix(FORMULA, meta)
  
  colnames(design) <- gsub(":", "X", colnames(design))
  
  if (identical(rownames(design), colnames(dge[[i]]$counts)) == FALSE){
    stop(paste0(celltype, ": subjects(rows of meta) don't match with count columns"))
  }
  
  design.matrix.list[[i]] = design
  
}
###

saveRDS(design.matrix.list, "pbulk_DE/sample_groups/all_timepoints_covid_only/results/UnSort_days_onset_design.rds")



