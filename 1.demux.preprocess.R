#this is tested using R 3.6.1 on a high-performance computing node with 8 cores and at least 160 gb of ram. 
library("Seurat") #load Seurat 3.1
library("dplyr")
library("matrixStats")
library('tidyverse')

# 10x lanes for import
# B4Lanes = c("01","02","03","04","05","06","07","08")
B4Lanes = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15",
            "17","18","19","20","21","22","23","24")

B4data = list()
B4CSPdata = list()
B4SeuratObj = list()

### Import the data
# the code below assumes the cellranger output data has been generated under a directory "Counts"
for(i in 1:length(B4Lanes)){
  B4data[[i]] = Read10X(paste("Counts/COVID_CITE_RNA5P",B4Lanes[[i]],"/outs/filtered_feature_bc_matrix", sep=""))
  B4SeuratObj[[i]] <- CreateSeuratObject(counts = B4data[[i]]$'Gene Expression', assay = "RNA", min.feature = 5)
  B4SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B4data[[i]]$`Antibody Capture`[11:203, ])
  B4SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B4data[[i]]$`Antibody Capture`[c(1:2), ])
  B4SeuratObj[[i]] <- RenameCells(B4SeuratObj[[i]], new.names = paste(substr(colnames(B4SeuratObj[[i]]), start = 1, stop = 16), B4Lanes[[i]], sep = ""))
}

names(B4data) = names(B4SeuratObj) = B4Lanes

#merge objects
B4merge1 = merge(B4SeuratObj[[1]], B4SeuratObj[2:8])
B4merge2 = merge(B4SeuratObj[[9]], B4SeuratObj[10:15])
B4merge3 = merge(B4SeuratObj[[16]], B4SeuratObj[18:23])

B4merge1$Batch <- rep("B4chip1", length(colnames(B4merge1)))
B4merge1$Sorted <- rep("UnSort", length(colnames(B4merge1)))
B4merge2$Batch <- rep("B4chip2", length(colnames(B4merge2)))
B4merge2$Sorted <- rep("UnSort", length(colnames(B4merge2)))
B4merge3$Batch <- rep("B4chip3", length(colnames(B4merge3)))
B4merge3$Sorted <- rep("Sort", length(colnames(B4merge3)))
B4merge <- merge(merge(B4merge1, B4merge2),B4merge3)

# QC and filter
quantile(B4merge$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))

quantile(B4merge$nCount_HTO, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))

B4merge <- subset(B4merge, subset = nCount_HTO < 50000)
B4merge[["percent.mt"]] <- PercentageFeatureSet(B4merge, pattern = "^MT-")
head(B4merge@meta.data, 5)
VlnPlot(B4merge, features = c("nFeature_RNA", "percent.mt", "nCount_HTO"), ncol = 3, pt.size =0)

B4merge <- subset(B4merge, subset = nFeature_RNA < 4000 & percent.mt < 20)

# Do HTODemux to get the negatives list for DSB normalization
B4merge[["HTO"]] <- SetAssayData(B4merge[["HTO"]], slot = "data", new.data = log(as.matrix(GetAssayData(B4merge[["HTO"]], slot = "counts"))+10))
B4merge <- ScaleData(B4merge, assay = "HTO", model.use = "negbinom")
B4merge = HTODemux(B4merge, positive.quantile = 0.95)
table(B4merge$HTO_classification.global)
Idents(B4merge) <- "hash.ID"

RidgePlot(B4merge, assay = "HTO", features = levels(B4merge$HTO_classification), ncol = 2)
FeatureScatter(B4merge, feature1="hto_HTO1", feature2="hto_HTO2")

# bring in SNP demuxlet calls from input folder
B4demuxbestList = list()
for(i in 1:length(B4Lanes)){
  B4demuxbestList[[i]] = read.table(paste("input/demuxlet/COVID_CITE_RNA5P",B4Lanes[[i]],"/demuxed.best", sep = ""), sep = "\t", header = TRUE)
}

for(i in 1:length(B4Lanes)){
  B4demuxbestList[[i]]$NewBarcode = paste(substr(B4demuxbestList[[i]]$BARCODE, start = 1, stop = 16),B4Lanes[[i]], sep = "")
}

B4demuxbestdf <- plyr::ldply(B4demuxbestList, data.frame) %>%
  mutate(SNG.BEST.GUESS = sapply(str_split(.$SNG.BEST.GUESS, pattern = "/"), function(x)x[12])) %>%
  mutate(SNG.BEST.GUESS = sapply(str_split(.$SNG.BEST.GUESS, pattern = "_"), function(x)x[1]))

length(which(colnames(B4merge) %in% B4demuxbestdf$NewBarcode))
setdiff(colnames(B4merge), B4demuxbestdf$NewBarcode)
B4merge <- B4merge[,-which(colnames(B4merge) %in% setdiff(colnames(B4merge), B4demuxbestdf$NewBarcode))]
setdiff(colnames(B4merge), B4demuxbestdf$NewBarcode)

rownames(B4demuxbestdf) <- B4demuxbestdf$NewBarcode

B4merge <- AddMetaData(B4merge, metadata = B4demuxbestdf[colnames(B4merge),])

B4merge.AMB <- subset(B4merge, subset = DROPLET.TYPE == "AMB")
B4merge.SNG <- subset(B4merge, subset = DROPLET.TYPE == "SNG")

# Do HTODemux to get the negatives list for DSB normalization
B4merge.SNG[["HTO"]] <- SetAssayData(B4merge.SNG[["HTO"]], slot = "data", new.data = log(as.matrix(GetAssayData(B4merge.SNG[["HTO"]], slot = "counts"))+10))
B4merge.SNG <- ScaleData(B4merge.SNG, assay = "HTO", model.use = "negbinom")
B4merge.SNG = HTODemux(B4merge.SNG, positive.quantile = 0.99)
Idents(B4merge.SNG) <- "hash.ID"
RidgePlot(B4merge.SNG, assay = "HTO", features = levels(B4merge.SNG$HTO_classification), ncol = 2)

table(B4merge.SNG$hash.ID, B4merge.SNG$DROPLET.TYPE)
B4merge.SNG$BEST.GUESS <- droplevels(B4merge.SNG$BEST.GUESS)
B4merge.SNG$SNG.BEST.GUESS <- droplevels(B4merge.SNG$SNG.BEST.GUESS)
table(B4merge.SNG$hash.ID, B4merge.SNG$SNG.BEST.GUESS)

# get the negatives list for DSB normalization
B4NegativesObj = subset(B4merge.SNG, subset = nFeature_RNA < 200 & nCount_CITE >0 & hash.ID == "Negative")
quantile(B4merge.SNG$nFeature_RNA, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
quantile(B4merge.SNG$nCount_CITE, probs = c(0, .1, .25, .5, .75, .95, .99, .999, .9999, 1))
B4merge.SNG = subset(B4merge.SNG, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 20 & nCount_CITE < 200000 & DROPLET.TYPE == "SNG")
VlnPlot(B4merge.SNG, features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mt"), pt.size = 0)


# DSB Normalization
isotype.control.name.vec = c("IgG1Kiso", "IgG2aKiso", "IgG2bKiso", "ratIgG2bKiso")
source("util_funs/dsb_normalization_functions.R")

neg_adt_matrix = as.matrix(GetAssayData(B4NegativesObj[["CITE"]], slot = "counts"))
positive_adt_matrix = as.matrix(GetAssayData(B4merge.SNG[["CITE"]], slot = "counts"))
limma.norm.ADT.list.B4merge = DSBNormalizeProtein(cell.columns.protein.matrix = positive_adt_matrix,
                                                  control.protein.matrix = neg_adt_matrix,
                                                  define.pseudocount  = TRUE, pseudocount.use = 10, 
                                                  isotype.control.name.vec =  isotype.control.name.vec)
B4merge.SNG[["CITE"]] <- SetAssayData(B4merge.SNG[["CITE"]], slot = "data", new.data = limma.norm.ADT.list.B4merge$denoised_adt)
B4merge.SNG <- AddMetaData(object = B4merge.SNG, metadata = limma.norm.ADT.list.B4merge$cellwise_background_mean, col.name = "ADTmclust1mean")
B4merge.SNG = subset(B4merge.SNG, subset = ADTmclust1mean < as.numeric(quantile(B4merge.SNG$ADTmclust1mean, 0.99)))
RidgePlot(B4merge.SNG, assay = "CITE", features = c("CD4","CD19","CD8","S1probe"), ncol = 2)

# ReDo HTODemux on DSB norm HTO
FeatureScatter(B4merge.SNG, feature1="hto_HTO1", feature2="hto_HTO2")
table(B4merge.SNG$hash.ID,B4merge.SNG$DROPLET.TYPE)
hist(log1p(GetAssayData(B4merge.SNG[["HTO"]], slot = "counts")["HTO1",]+9), breaks = 100)
hist(log1p(GetAssayData(subset(B4merge.SNG, subset = hash.ID == "HTO1")[["HTO"]], slot = "counts")["HTO1",]+9), breaks = 100)
hist(log1p(GetAssayData(B4merge.SNG[["HTO"]], slot = "counts")["HTO2",]+9), breaks = 100)
hist(log1p(GetAssayData(subset(B4merge.SNG, subset = hash.ID == "HTO2")[["HTO"]], slot = "counts")["HTO2",]+9), breaks = 100)
VlnPlot(B4merge.SNG, features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mt"), pt.size = 0)

##### manual thresholding hash calls
HTONames = rownames(B4merge.SNG[["HTO"]])
B4merge.SNG <- ScaleData(B4merge.SNG, assay = "HTO")

pdf("plots/HTOdist.check.B4.pdf")
for (i in 1:2) {
  hist(as.matrix(GetAssayData(B4merge.SNG[["HTO"]], slot = "data"))[i,], breaks=100)
}
dev.off()


# manual threshold HTO calls
hithres = c(5.4, 5.9)
lothres = c(5.35, 5.85)
automaticHTOsHiThres = character(length=length(colnames(B4merge.SNG)))
for(i in 1:length(colnames(B4merge.SNG))){
  automaticHTOsHiThres[i] = paste(HTONames[GetAssayData(B4merge.SNG[["HTO"]], slot = "data")[,i] > hithres], collapse="+")
}

# manual lower threshold HTO calls, for doublet removal
automaticHTOsloThres = character(length=length(colnames(B4merge.SNG)))
for(i in 1:length(colnames(B4merge.SNG))){
  automaticHTOsloThres[i] = paste(HTONames[GetAssayData(B4merge.SNG[["HTO"]], slot = "data")[,i] > lothres], collapse="+")
}

singlets = automaticHTOsHiThres %in% HTONames & automaticHTOsloThres %in% HTONames
autoHashcalls = ifelse(singlets, automaticHTOsHiThres, "NonSinglet")
B4merge.SNG <- AddMetaData(object = B4merge.SNG, metadata = autoHashcalls, col.name = "autoHashcalls")
Idents(B4merge.SNG) <- "hash.ID"
FeatureScatter(subset(B4merge.SNG, idents = c("HTO1","HTO2")), feature1="hto_HTO1", feature2="hto_HTO2")
RidgePlot(B4merge.SNG, assay = "HTO", features = c("HTO1","HTO2"), ncol = 2)

Idents(B4merge.SNG) <- "autoHashcalls"
FeatureScatter(subset(B4merge.SNG, idents = c("HTO1","HTO2")), feature1="hto_HTO1", feature2="hto_HTO2")
table(B4merge.SNG$hash.ID, B4merge.SNG$autoHashcalls)
table(B4merge.SNG$autoHashcalls)
RidgePlot(B4merge.SNG, assay = "HTO", features = c("HTO1","HTO2"), ncol = 2)

B4merge.SNG.NonSinglet <- subset(B4merge.SNG, idents = "NonSinglet")
B4merge.SNG.singlet <- subset(B4merge.SNG, idents = "NonSinglet", invert = TRUE)
Idents(B4merge.SNG.singlet) <- "hash.ID"
RidgePlot(B4merge.SNG.singlet, assay = "HTO", features = c("HTO1", "HTO2"), ncol = 2)
Idents(B4merge.SNG.singlet) <- "autoHashcalls"
RidgePlot(B4merge.SNG.singlet, assay = "HTO", features = c("HTO1", "HTO2"), ncol = 2)

### add in metadata ------------------------------------------------------------------------
B4merge.SNG$sample <- paste(B4merge.SNG$SNG.BEST.GUESS, B4merge.SNG$autoHashcalls, sep = "_")
Idents(B4merge.SNG) <- "sample"
table(Idents(B4merge.SNG))
write.csv(unique(B4merge.SNG$sample), "B4sample.csv")

samplenames = unique(Idents(B4merge.SNG))
samplemetadata = read.csv("input/B4sample.manualhto.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
identical(samplemetadata$x, sort(unique(B4merge.SNG$sample)))
B4merge.SNG$Age <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Age)
B4merge.SNG$Gender <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Gender)
B4merge.SNG$Status <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Status)
B4merge.SNG$Pool <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Pool)
B4merge.SNG$Timepoint <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Timepoint)
B4merge.SNG$Subject <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Subject)
B4merge.SNG$Class <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$Class)
B4merge.SNG$days_since_admission <- plyr::mapvalues(x = Idents(B4merge.SNG), from = samplemetadata$x, to = samplemetadata$days_since_admission)
B4merge.SNG$days_since_admission <- as.numeric(as.character(B4merge.SNG$days_since_admission))
B4merge.SNG$sample_id <- paste(B4merge.SNG$Subject, B4merge.SNG$Timepoint, sep = "_")

B4merge.SNG.withmissassigned <- B4merge.SNG
B4merge.SNG <- subset(B4merge.SNG, subset = Pool == "misassigned", invert = TRUE)
table(B4merge.SNG$sample_id, B4merge.SNG$Pool) 

saveRDS(object = B4merge.SNG, file = "B4merge.SNG.preclust.rds")
saveRDS(object = B4merge.SNG.withmissassigned, file = "B4merge.wmisassigned.preclust.rds")


