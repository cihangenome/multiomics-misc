
library(limma)
library(variancePartition)
library(fgsea)
library(data.table)

# run fgsea for the pseudobulk DE results
#Set paths ---------------------------------------------------------------
TOPTAB_IN_PATH <- "pbulk_DE/sample_groups/"

SAMPLE_GROUP <- c("tso_within_d40_misc_plus_healthy/results/", 
                  "tso_within_d40_covid_plus_healthy/results/", 
                  "tso_within_d40_misc_plus_covid/results/", 
                  "all_timepoints_misc_only/results/", 
                  "all_timepoints_covid_only/results/")

COMBINED.GENESETS.IN.PATH <- "input/kegg_go_btm_reactome_foointerferon.rds"


#Load data ---------------------------------------------------------------
geneset.list <- readRDS(COMBINED.GENESETS.IN.PATH)

#########################
### MIS-C vs HC d40 #####
#########################

toptab_Unsort <- readRDS(paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[1], "UnSort_misc_vs_healthy_toptab.rds"))

toptab <- toptab_Unsort

celltypes = names(toptab)
fgseares.list = list()
FORMULA <- "~0 + Class + Age"
CONTRAST_NAME <- "MISC-HC"

for (i in celltypes) {
  
  t.stat <- toptab[[i]]$t
  names(t.stat) <- toptab[[i]]$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = t.stat,
                    minSize=10,
                    maxSize=500,
                    nperm=50000)
  topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  TSV.OUT.PATH <- paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[1], "fgsea_tables/UnSort/", i, "--", CONTRAST_NAME, "--fgsea.tsv")

  fwrite(fgseaRes, file=TSV.OUT.PATH, sep="\t", sep2=c("", " ", ""), nThread = 1)
  fgseares.list[[i]] <- fgseaRes
}

saveRDS(fgseares.list, paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[1], "fgsea_tables/UnSort_misc_vs_healthy_fgseares.rds"))


#########################
### COVID vs HC d40 #####
#########################

toptab_Unsort <- readRDS(paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[2], "UnSort_covid_vs_healthy_toptab.rds"))

toptab <- toptab_Unsort

celltypes = names(toptab)
fgseares.list = list()
FORMULA <- "~0 + Class + Age"
CONTRAST_NAME <- "COVID-HC"

for (i in celltypes) {
  
  t.stat <- toptab[[i]]$t
  names(t.stat) <- toptab[[i]]$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = t.stat,
                    minSize=10,
                    maxSize=500,
                    nperm=50000)
  topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  TSV.OUT.PATH <- paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[2], "fgsea_tables/UnSort/", i, "--", CONTRAST_NAME, "--fgsea.tsv")

  fwrite(fgseaRes, file=TSV.OUT.PATH, sep="\t", sep2=c("", " ", ""), nThread = 1)
  fgseares.list[[i]] <- fgseaRes
}

saveRDS(fgseares.list, paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[2], "fgsea_tables/UnSort_covid_vs_healthy_fgseares.rds"))


##########################
### MIS-C vs COVID d40 ###
##########################

toptab_Unsort <- readRDS(paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[3], "UnSort_misc_vs_covid_toptab.rds"))

toptab <- toptab_Unsort

celltypes = names(toptab)
fgseares.list = list()
FORMULA <- "~0 + Class + days_since_admission + Age"
CONTRAST_NAME <- "MISC-COVID"

for (i in celltypes) {
  
  t.stat <- toptab[[i]]$t
  names(t.stat) <- toptab[[i]]$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = t.stat,
                    minSize=10,
                    maxSize=500,
                    nperm=50000)
  topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  TSV.OUT.PATH <- paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[3], "fgsea_tables/UnSort/", i, "--", CONTRAST_NAME, "--fgsea.tsv")

  fwrite(fgseaRes, file=TSV.OUT.PATH, sep="\t", sep2=c("", " ", ""), nThread = 1)
  fgseares.list[[i]] <- fgseaRes
}

saveRDS(fgseares.list, paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[3], "fgsea_tables/UnSort_misc_vs_covid_fgseares.rds"))


#######################################
### MIS-C only all timepoints #########
#######################################

toptab_Unsort <- readRDS(paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[4], "UnSort_days_onset_toptab.rds"))

toptab <- toptab_Unsort

celltypes = names(toptab)
fgseares.list = list()
FORMULA <- "~days_since_admission + Age"
CONTRAST_NAME <- "days_since_admission"

for (i in celltypes) {
  
  t.stat <- toptab[[i]]$t
  names(t.stat) <- toptab[[i]]$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = t.stat,
                    minSize=10,
                    maxSize=500,
                    nperm=50000)
  topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  TSV.OUT.PATH <- paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[4], "fgsea_tables/UnSort/", i, "--", CONTRAST_NAME, "--fgsea.tsv")

  fwrite(fgseaRes, file=TSV.OUT.PATH, sep="\t", sep2=c("", " ", ""), nThread = 1)
  fgseares.list[[i]] <- fgseaRes
}

saveRDS(fgseares.list, paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[4], "fgsea_tables/UnSort_days_onset_fgseares.rds"))


#######################################
### COVID only all timepoints #########
#######################################

toptab_Unsort <- readRDS(paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[5], "UnSort_days_onset_toptab.rds"))

toptab <- toptab_Unsort
toptab <- toptab_All

celltypes = names(toptab)
fgseares.list = list()
FORMULA <- "~days_since_admission + Age"
CONTRAST_NAME <- "days_since_admission"

for (i in celltypes) {
  
  t.stat <- toptab[[i]]$t
  names(t.stat) <- toptab[[i]]$gene
  
  set.seed(1)
  fgseaRes <- fgsea(pathways = geneset.list, 
                    stats = t.stat,
                    minSize=10,
                    maxSize=500,
                    nperm=50000)
  topPathwaysUp <- fgseaRes[ES > 0][head(rev(order(NES)), n=20), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(NES), n=20), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  TSV.OUT.PATH <- paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[5], "fgsea_tables/UnSort/", i, "--", CONTRAST_NAME, "--fgsea.tsv")

  fwrite(fgseaRes, file=TSV.OUT.PATH, sep="\t", sep2=c("", " ", ""), nThread = 1)
  fgseares.list[[i]] <- fgseaRes
}

saveRDS(fgseares.list, paste0(TOPTAB_IN_PATH, SAMPLE_GROUP[5], "fgsea_tables/UnSort_days_onset_fgseares.rds"))












