library(limma)
library(edgeR)
library(tidyverse)

# Make contrast matrix for each celltype 
# contrast fit
# get top table of DE genes
##########################################################################
### make contrast matrix -- contrast fit -- get top table of DE genes ####
##########################################################################

DESIGN_MATRIX_PATH <- "output/pbulk_DE/sample_groups/"
SAMPLE_GROUP <- c("tso_within_d40_misc_plus_healthy/", 
                  "tso_within_d40_covid_plus_healthy/",
                  "tso_within_d40_misc_plus_covid/",
                  "all_timepoints_misc_only/", 
                  "all_timepoints_covid_only/")


#########################
### MIS-C vs HC d40 #####
#########################
design_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_design.rds"))
fit_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_fit.rds"))

#---1
design <- design_Unsort
fit <- fit_Unsort

CONTRAST <- "ClassMISC - ClassHC"
CONTRAST_NAME <- "MISC-HC"

celltypes = names(design)
contrast.matrix.list = list()
cfit.list = list()
toptab.list = list()

for (i in celltypes) {

contrast_mat <- makeContrasts(contrasts = CONTRAST, levels = colnames(design[[i]])) 
colnames(contrast_mat) <- CONTRAST_NAME

contrast.matrix.list[[i]] = contrast_mat

# contrast_fit
cfit <- contrasts.fit(fit[[i]], contrasts = contrast_mat)
cfit <- eBayes(cfit)

cfit.list[[i]] <- cfit

toptab <- topTable(cfit, number = nrow(cfit), coef = CONTRAST_NAME) %>%
  rownames_to_column(var = "gene")

toptab.list[[i]] <- toptab
}

#---1
saveRDS(contrast.matrix.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_contrast_mat.rds"))
saveRDS(cfit.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_cfit.rds"))
saveRDS(toptab.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_toptab.rds"))



#########################
### COVID vs HC d40 #####
#########################
design_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_design.rds"))
fit_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_fit.rds"))

#---1
design <- design_Unsort
fit <- fit_Unsort

CONTRAST <- "ClassCOVID - ClassHC"
CONTRAST_NAME <- "COVID-HC"

celltypes = names(design)
contrast.matrix.list = list()
cfit.list = list()
toptab.list = list()

for (i in celltypes) {
  
  contrast_mat <- makeContrasts(contrasts = CONTRAST, levels = colnames(design[[i]])) 
  colnames(contrast_mat) <- CONTRAST_NAME
  
  contrast.matrix.list[[i]] = contrast_mat
  
  # contrast_fit
  cfit <- contrasts.fit(fit[[i]], contrasts = contrast_mat)
  cfit <- eBayes(cfit)
  
  cfit.list[[i]] <- cfit
  
  toptab <- topTable(cfit, number = nrow(cfit), coef = CONTRAST_NAME) %>%
    rownames_to_column(var = "gene")
  
  toptab.list[[i]] <- toptab
}

#---1
saveRDS(contrast.matrix.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_contrast_mat.rds"))
saveRDS(cfit.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_cfit.rds"))
saveRDS(toptab.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_toptab.rds"))



##########################
### MIS-C vs COVID d40 ###
##########################
design_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_design.rds"))
fit_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_fit.rds"))

#---1
design <- design_Unsort
fit <- fit_Unsort

CONTRAST <- "ClassMISC - ClassCOVID"
CONTRAST_NAME <- "MISC-COVID"

celltypes = names(design)
contrast.matrix.list = list()
cfit.list = list()
toptab.list = list()

for (i in celltypes) {
  
  contrast_mat <- makeContrasts(contrasts = CONTRAST, levels = colnames(design[[i]])) 
  colnames(contrast_mat) <- CONTRAST_NAME
  
  contrast.matrix.list[[i]] = contrast_mat
  
  # contrast_fit
  cfit <- contrasts.fit(fit[[i]], contrasts = contrast_mat)
  cfit <- eBayes(cfit)
  
  cfit.list[[i]] <- cfit
  
  toptab <- topTable(cfit, number = nrow(cfit), coef = CONTRAST_NAME) %>%
    rownames_to_column(var = "gene")
  
  toptab.list[[i]] <- toptab
}

#---1
saveRDS(contrast.matrix.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_contrast_mat.rds"))
saveRDS(cfit.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_cfit.rds"))
saveRDS(toptab.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_toptab.rds"))



#######################################
### MIS-C only all timepoints #########
#######################################
design_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_design.rds"))
fit_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_fit.rds"))

#---1
design <- design_Unsort
fit <- fit_Unsort

CONTRAST <- "days_since_admission"
CONTRAST_NAME <- "days_since_admission"

celltypes = names(fit)
contrast.matrix.list = list()
cfit.list = list()
toptab.list = list()

for (i in celltypes) {
  
  contrast_mat <- makeContrasts(contrasts = CONTRAST, levels = colnames(design[[i]])) 
  colnames(contrast_mat) <- CONTRAST_NAME
  
  contrast.matrix.list[[i]] = contrast_mat
  
  # contrast_fit
  cfit <- contrasts.fit(fit[[i]], contrasts = contrast_mat)
  cfit <- eBayes(cfit)
  
  cfit.list[[i]] <- cfit
  
  toptab <- topTable(cfit, number = nrow(cfit), coef = CONTRAST_NAME) %>%
    rownames_to_column(var = "gene")
  
  toptab.list[[i]] <- toptab
}

#---1
saveRDS(contrast.matrix.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_contrast_mat.rds"))
saveRDS(cfit.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_contrast_cfit.rds"))
saveRDS(toptab.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_toptab.rds"))


#######################################
### COVID only all timepoints #########
#######################################
design_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_design.rds"))
fit_Unsort <- readRDS(paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_fit.rds"))

#---1
design <- design_Unsort
fit <- fit_Unsort

CONTRAST <- "days_since_admission"
CONTRAST_NAME <- "days_since_admission"

celltypes = names(fit)
contrast.matrix.list = list()
cfit.list = list()
toptab.list = list()

for (i in celltypes) {
  
  contrast_mat <- makeContrasts(contrasts = CONTRAST, levels = colnames(design[[i]])) 
  colnames(contrast_mat) <- CONTRAST_NAME
  
  contrast.matrix.list[[i]] = contrast_mat
  
  # contrast_fit
  cfit <- contrasts.fit(fit[[i]], contrasts = contrast_mat)
  cfit <- eBayes(cfit)
  
  cfit.list[[i]] <- cfit
  
  toptab <- topTable(cfit, number = nrow(cfit), coef = CONTRAST_NAME) %>%
    rownames_to_column(var = "gene")
  
  toptab.list[[i]] <- toptab
}

#---1
saveRDS(contrast.matrix.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_contrast_mat.rds"))
saveRDS(cfit.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_contrast_cfit.rds"))
saveRDS(toptab.list, paste0(DESIGN_MATRIX_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_toptab.rds"))





