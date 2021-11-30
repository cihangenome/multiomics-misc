

packages <-  c("GENIE3","Biobase","pheatmap","writexl","readxl","ggpubr", "dplyr","psych", "corrplot", "ComplexHeatmap","grDevices","reshape2","immunarch","tidyverse","readr","ggplot2","stringr","ggtext","gtools")
invisible(lapply(packages, library, character.only = TRUE))

df_v112_and_biomarkers_timeline_cortico_ivig<-read_excel("biomarker_data_input.xlsx",sheet="treatment_timeline")

df_biomarkers<-read_excel("biomarker_data_input.xlsx",sheet="biomarker_levels")

df_v112<-read_excel("biomarker_data_input.xlsx",sheet ="TRBV11_2_usage")

####all samples are MIS-C, we can remove the column with the diagnosis
df_v112_and_biomarkers_timeline_cortico_ivig<-df_v112_and_biomarkers_timeline_cortico_ivig[,-c(2)]

####indices of the samples with no corticosteroid treatment timeline data
ind_out<-grep("never",df_v112_and_biomarkers_timeline_cortico_ivig$Glucocorticoid_interval)

####samples with corticosteroid treatment timeline data
df_samples_with_cortico_dates<-df_v112_and_biomarkers_timeline_cortico_ivig[-ind_out,c(1,3)]


####indices of the samples with missing corticosteroid or ivig treatment timeline data
ind_out1<-grep("never",df_v112_and_biomarkers_timeline_cortico_ivig$Glucocorticoid_interval)
ind_out2<-grep("never",df_v112_and_biomarkers_timeline_cortico_ivig$IVIG_interval)
####samples with both corticosteroid and ivig treatment timeline dates
df_samples_with_cortico_ivig_dates<-df_v112_and_biomarkers_timeline_cortico_ivig[-union(ind_out1,ind_out2),c(1,3,5)]

###to find samples with corticosteroid but no ivig dates
samples_with_cortico_no_ivig<-setdiff(df_samples_with_cortico_dates$Sample_name,df_samples_with_cortico_ivig_dates$Sample_name)
ind_keep<- which( df_v112_and_biomarkers_timeline_cortico_ivig$Sample_name   %in% samples_with_cortico_no_ivig )
####samples with cortico but no ivig dates
df_samples_with_cortico_no_ivig_dates<-df_v112_and_biomarkers_timeline_cortico_ivig[ind_keep,c(1,3)]



df_v112_and_biomarkers<-merge(df_biomarkers,df_v112[,c(1:2)],by="Sample_name",all=TRUE)

df_v112_and_biomarkers<-merge(df_v112_and_biomarkers,df_v112_and_biomarkers_timeline_cortico_ivig[,c(1,6)],by="Sample_name")

###Expand the names of two columns
colnames(df_v112_and_biomarkers)[grep("Usage",colnames(df_v112_and_biomarkers))]<-"TRBV11-2 usage"
colnames(df_v112_and_biomarkers)[match("interval",colnames(df_v112_and_biomarkers))]<-"time interval"

df_v112_and_biomarkers<-df_v112_and_biomarkers[,c(1,53,2,54,3:52)]

df_usage_biomarker_cols<-data.frame(cols=colnames(df_v112_and_biomarkers))

##columns to log10 transform
more_cols<-setdiff(c(1:nrow(df_usage_biomarker_cols)),c(1:4))

###for samples with no gene usage data, set the NAs to zero so that the data can be transformed
df_v112_and_biomarkers[is.na(df_v112_and_biomarkers)] <- 0

#log10 transformation the biomarker levels
df_v112_and_biomarkers[,more_cols]<-log10(df_v112_and_biomarkers[,more_cols])

df_v112_and_biomarkers_selec<-df_v112_and_biomarkers[,c(1,2,4,more_cols)]

df_v112_and_biomarkers_selec[is.na(df_v112_and_biomarkers_selec)] <- 0

####combine the biomarkers with timeline data
df_v112_and_biomarkers_selec<-merge(df_v112_and_biomarkers_selec,df_v112_and_biomarkers_timeline_cortico_ivig,by="Sample_name")

colnames(df_v112_and_biomarkers_selec)[grep("_interval",colnames(df_v112_and_biomarkers_selec))]<-gsub("_interval"," interval",colnames(df_v112_and_biomarkers_selec)[grep("_interval",colnames(df_v112_and_biomarkers_selec))])


#this is where we start processing df_v112_and_biomarkers_selec with the transpose step
#remove the sample name column

sample_names_to_record<-df_v112_and_biomarkers_selec$Sample_name

df_v112_and_biomarkers_selec<-df_v112_and_biomarkers_selec[,-c(1)]
df_v112_and_biomarkers_selec<-t(df_v112_and_biomarkers_selec)
colnames(df_v112_and_biomarkers_selec)<-sample_names_to_record

####this step is very important, set the Inf entries to 0
df_v112_and_biomarkers_selec[grep("Inf",df_v112_and_biomarkers_selec)]<-0


#####Start subsetting df_v112_and_biomarkers_selec here for separate analyses

####which samples have cortico dates?
ind_col_cortico_dates<-which(colnames(df_v112_and_biomarkers_selec) %in% df_samples_with_cortico_dates$Sample_name)

###these are the samples for which we have TRBV gene usage data
ind_nonzero_usage_samples<-which(as.numeric(df_v112_and_biomarkers_selec[1,])>0)

dum<-df_v112_and_biomarkers_selec[c(1:52,54),intersect(ind_col_cortico_dates,ind_nonzero_usage_samples)]
gene_names<-rownames(dum)

dum<-as.data.frame(dum)
dum<-as.matrix(sapply(dum, as.numeric))

df_v112_and_biomarkers_selec_cortico_scaled<-scale(dum, center = FALSE, scale=TRUE)
rownames(df_v112_and_biomarkers_selec_cortico_scaled)<-gene_names

#The GENIE3 algorithm outputs a matrix containing the weights of the putative regulatory links, with higher #weights corresponding to more likely regulatory links. weightMat[i,j] is the weight of the link directed #from the i-th gene to j-th gene

set.seed(123)
#Glucocorticoid interval,time interval,TRBV11-2 usage are the non-biomarker variables
network_cortico_scaled<- GENIE3(df_v112_and_biomarkers_selec_cortico_scaled,nTrees = 1000)
links_network_cortico_scaled<-getLinkList(network_cortico_scaled)

correl_cortico_scaled<-corr.test(t(df_v112_and_biomarkers_selec_cortico_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)
correl_cortico_scaled<-as.data.frame(correl_cortico_scaled[["r"]])

func_sign <- function(a) {
  ifelse(a > 0, 1, -1)  #'1' and '-1'
}
correl_cortico_scaled_sign<-apply(correl_cortico_scaled,c(1,2), func_sign)

links_network_cortico_scaled_target_usage<-links_network_cortico_scaled[which(links_network_cortico_scaled$targetGene=="TRBV11-2 usage"),]
thres_0p5<-quantile(links_network_cortico_scaled_target_usage$weight, probs = c(0.50))
weight_thres<-as.data.frame(thres_0p5)$thres_0p5
links_network_cortico_scaled_target_usage_selec<-links_network_cortico_scaled_target_usage[which(links_network_cortico_scaled_target_usage$weight>=weight_thres),]
selected_interactors_cortico<-c(as.character(links_network_cortico_scaled_target_usage_selec$regulatoryGene),"TRBV11-2 usage")

network_cortico_scaled_to_plot<-network_cortico_scaled[selected_interactors_cortico,selected_interactors_cortico]
network_cortico_scaled_to_plot[which(network_cortico_scaled_to_plot=="-Inf")]<-0

# Define color palette
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")


p_network_cortico_scaled<-ggballoonplot(network_cortico_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: 92 samples with steroid treatment timeline data",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)   #c("blue", "white", "red")
 #\n, , inputs are scaled

Fig_Ext_5G<-ggpar(p_network_cortico_scaled,font.xtickslab= c(10, "plain", "black"),font.ytickslab= c(10, "plain", "black"), orientation = "horiz",font.main = c(10,"bold"),font.legend = c(8))
Fig_Ext_5G<-Fig_Ext_5G[["data"]]
colnames(Fig_Ext_5G)[c(1,2)]<-c("predictor","target")

pdf(file = "Figure_Ext_5G.pdf", width=9, height=6)
ggpar(p_network_cortico_scaled,font.xtickslab= c(10, "plain", "black"),font.ytickslab= c(10, "plain", "black"), orientation = "horiz",font.main = c(10,"bold"),font.legend = c(8))
dev.off()

col2 = colorRampPalette(c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF"))

pdf(file = "Figure_Ext_5F.pdf", width=10, height=8)
corrplot(as.matrix(correl_cortico_scaled[selected_interactors_cortico,selected_interactors_cortico]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficient values (92 samples with steroid treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)
dev.off()


Fig_Ext_5f_data<-corrplot(as.matrix(correl_cortico_scaled[selected_interactors_cortico,selected_interactors_cortico]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficient values (92 samples with steroid treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)
Fig_Ext_5f_data<-as.data.frame(Fig_Ext_5f_data[["corrPos"]])



Fig_Ext_5g_data<-ggpar(p_network_cortico_scaled,font.xtickslab= c(10, "plain", "black"),font.ytickslab= c(10, "plain", "black"), orientation = "horiz",font.main = c(10,"bold"),font.legend = c(8))
Fig_Ext_5g_data<-Fig_Ext_5g_data[["data"]]
colnames(Fig_Ext_5g_data)[c(1,2)]<-c("predictor","target")


#################
#################
ind_col_cortico_ivig_dates<-which(colnames(df_v112_and_biomarkers_selec) %in% df_samples_with_cortico_ivig_dates$Sample_name)


dum<-df_v112_and_biomarkers_selec[c(1:52,54,56),intersect(ind_col_cortico_ivig_dates,ind_nonzero_usage_samples)]
gene_names<-rownames(dum)

dum<-as.data.frame(dum)
dum<-as.matrix(sapply(dum, as.numeric))


df_v112_and_biomarkers_selec_cortico_ivig_scaled<-scale(dum, center = FALSE, scale=TRUE)
rownames(df_v112_and_biomarkers_selec_cortico_ivig_scaled)<-gene_names

set.seed(123)
#Glucocorticoid interval,time interval,IVIG interval,TRBV11-2 usage are non-biomarker variables
network_cortico_ivig_scaled<- GENIE3(df_v112_and_biomarkers_selec_cortico_ivig_scaled,nTrees = 1000)
links_network_cortico_ivig_scaled<-getLinkList(network_cortico_ivig_scaled)

correl_cortico_ivig_scaled<-corr.test(t(df_v112_and_biomarkers_selec_cortico_ivig_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)
correl_cortico_ivig_scaled<-as.data.frame(correl_cortico_ivig_scaled[["r"]])


correl_cortico_ivig_scaled_sign<-apply(correl_cortico_ivig_scaled,c(1,2), func_sign)

links_network_cortico_ivig_scaled_target_usage<-links_network_cortico_ivig_scaled[which(links_network_cortico_ivig_scaled$targetGene=="TRBV11-2 usage"),]

thres_0p5<-quantile(links_network_cortico_ivig_scaled_target_usage$weight, probs = c(0.50))

weight_thres<-as.data.frame(thres_0p5)$thres_0p5

links_network_cortico_ivig_scaled_target_usage_selec<-links_network_cortico_ivig_scaled_target_usage[which(links_network_cortico_ivig_scaled_target_usage$weight>=weight_thres),]

selected_interactors_cortico_ivig<-c(as.character(links_network_cortico_ivig_scaled_target_usage_selec$regulatoryGene),"TRBV11-2 usage")

network_cortico_ivig_scaled_to_plot<-network_cortico_ivig_scaled[selected_interactors_cortico_ivig,selected_interactors_cortico_ivig]
network_cortico_ivig_scaled_to_plot[which(network_cortico_ivig_scaled_to_plot=="-Inf")]<-0


p_network_cortico_ivig_scaled<-ggballoonplot(network_cortico_ivig_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: 52 samples with steroid+IVIG treatment timeline data \n (Top 50th percentile predictors of TRBV11-2 usage)",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)

################################
######CORTICO+NO IVIG#######
################################

ind_col_cortico_no_ivig_dates<-which(colnames(df_v112_and_biomarkers_selec) %in% df_samples_with_cortico_no_ivig_dates$Sample_name)


dum<-df_v112_and_biomarkers_selec[c(1:52,54),intersect(ind_col_cortico_no_ivig_dates,ind_nonzero_usage_samples)]
gene_names<-rownames(dum)
dum<-as.data.frame(dum)


dum<-as.matrix(sapply(dum, as.numeric))

df_v112_and_biomarkers_selec_cortico_no_ivig_scaled<-scale(dum, center = FALSE, scale=TRUE)
rownames(df_v112_and_biomarkers_selec_cortico_no_ivig_scaled)<-gene_names

set.seed(123)

#Glucocorticoid interval,time interval,TRBV11-2 usage are the non-biomarker variables
network_cortico_no_ivig_scaled<- GENIE3(df_v112_and_biomarkers_selec_cortico_no_ivig_scaled,nTrees = 1000)
links_network_cortico_no_ivig_scaled<-getLinkList(network_cortico_no_ivig_scaled)


correl_cortico_no_ivig_scaled<-corr.test(t(df_v112_and_biomarkers_selec_cortico_no_ivig_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)
correl_cortico_no_ivig_scaled<-as.data.frame(correl_cortico_no_ivig_scaled[["r"]])


correl_cortico_no_ivig_scaled_sign<-apply(correl_cortico_no_ivig_scaled,c(1,2), func_sign)


links_network_cortico_no_ivig_scaled_target_usage<-links_network_cortico_no_ivig_scaled[which(links_network_cortico_no_ivig_scaled$targetGene=="TRBV11-2 usage"),]

thres_0p5<-quantile(links_network_cortico_no_ivig_scaled_target_usage$weight, probs = c(0.50))

weight_thres<-as.data.frame(thres_0p5)$thres_0p5

links_network_cortico_no_ivig_scaled_target_usage_selec<-links_network_cortico_no_ivig_scaled_target_usage[which(links_network_cortico_no_ivig_scaled_target_usage$weight>=weight_thres),]

selected_interactors_cortico_no_ivig<-c(as.character(links_network_cortico_no_ivig_scaled_target_usage_selec$regulatoryGene),"TRBV11-2 usage")

network_cortico_no_ivig_scaled_to_plot<-network_cortico_no_ivig_scaled[selected_interactors_cortico_no_ivig,selected_interactors_cortico_no_ivig]
network_cortico_no_ivig_scaled_to_plot[which(network_cortico_no_ivig_scaled_to_plot=="-Inf")]<-0


p_network_cortico_no_ivig_scaled<-ggballoonplot(network_cortico_no_ivig_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: 40 samples with steroid (but no IVIG) treatment timeline data \n (Top 50th percentile predictors of TRBV11-2 usage)",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)

####################################
####################################

correl_cortico_scaled_new<-corr.test(t(df_v112_and_biomarkers_selec_cortico_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)

correl_cortico_ivig_scaled_new<-corr.test(t(df_v112_and_biomarkers_selec_cortico_ivig_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)

correl_cortico_no_ivig_scaled_new<-corr.test(t(df_v112_and_biomarkers_selec_cortico_no_ivig_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)

links_network_cortico_scaled_regulatory_usage<-links_network_cortico_scaled[which(links_network_cortico_scaled$regulatoryGene=="TRBV11-2 usage"),]
thres_0p5<-quantile(links_network_cortico_scaled_regulatory_usage$weight, probs = c(0.50))
weight_thres<-as.data.frame(thres_0p5)$thres_0p5
links_network_cortico_scaled_regulatory_usage_selec<-links_network_cortico_scaled_regulatory_usage[which(links_network_cortico_scaled_regulatory_usage$weight>=weight_thres),]
selected_interactors_cortico_new<-c(as.character(links_network_cortico_scaled_regulatory_usage_selec$targetGene),"TRBV11-2 usage")



links_network_cortico_ivig_scaled_regulatory_usage<-links_network_cortico_ivig_scaled[which(links_network_cortico_ivig_scaled$regulatoryGene=="TRBV11-2 usage"),]
thres_0p5<-quantile(links_network_cortico_ivig_scaled_regulatory_usage$weight, probs = c(0.50))
weight_thres<-as.data.frame(thres_0p5)$thres_0p5
links_network_cortico_ivig_scaled_regulatory_usage_selec<-links_network_cortico_ivig_scaled_regulatory_usage[which(links_network_cortico_ivig_scaled_regulatory_usage$weight>=weight_thres),]
selected_interactors_cortico_ivig_new<-c(as.character(links_network_cortico_ivig_scaled_regulatory_usage_selec$targetGene),"TRBV11-2 usage")
length(intersect(selected_interactors_cortico_ivig_new,selected_interactors_cortico))/length(selected_interactors_cortico)


links_network_cortico_no_ivig_scaled_regulatory_usage<-links_network_cortico_no_ivig_scaled[which(links_network_cortico_no_ivig_scaled$regulatoryGene=="TRBV11-2 usage"),]
thres_0p5<-quantile(links_network_cortico_no_ivig_scaled_regulatory_usage$weight, probs = c(0.50))
weight_thres<-as.data.frame(thres_0p5)$thres_0p5
links_network_cortico_no_ivig_scaled_regulatory_usage_selec<-links_network_cortico_no_ivig_scaled_regulatory_usage[which(links_network_cortico_no_ivig_scaled_regulatory_usage$weight>=weight_thres),]
selected_interactors_cortico_no_ivig_new<-c(as.character(links_network_cortico_no_ivig_scaled_regulatory_usage_selec$targetGene),"TRBV11-2 usage")



dum<-df_v112_and_biomarkers_selec[c(2:52,54),ind_col_cortico_dates]
gene_names<-rownames(dum)

dum<-as.data.frame(dum)
dum<-as.matrix(sapply(dum, as.numeric))

df_v112_and_biomarkers_selec_cortico_nousage_scaled<-scale(dum, center = FALSE, scale=TRUE)
rownames(df_v112_and_biomarkers_selec_cortico_nousage_scaled)<-gene_names


set.seed(123)
#Glucocorticoid interval and time interval
network_cortico_nousage_scaled<- GENIE3(df_v112_and_biomarkers_selec_cortico_nousage_scaled,nTrees = 1000)
links_network_cortico_nousage_scaled<-getLinkList(network_cortico_nousage_scaled)

correl_cortico_nousage_scaled<-corr.test(t(df_v112_and_biomarkers_selec_cortico_nousage_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)
correl_cortico_nousage_scaled<-as.data.frame(correl_cortico_nousage_scaled[["r"]])


correl_cortico_nousage_scaled_sign<-apply(correl_cortico_nousage_scaled,c(1,2), func_sign)


links_network_cortico_nousage_scaled_target_usage<-links_network_cortico_nousage_scaled[which(links_network_cortico_nousage_scaled$targetGene=="sST2"),]
thres_0p5<-quantile(links_network_cortico_nousage_scaled_target_usage$weight, probs = c(0.00))
weight_thres<-as.data.frame(thres_0p5)$thres_0p5
links_network_cortico_nousage_scaled_target_usage_selec<-links_network_cortico_nousage_scaled_target_usage[which(links_network_cortico_nousage_scaled_target_usage$weight>=weight_thres),]
selected_interactors_cortico<-c(as.character(links_network_cortico_nousage_scaled_target_usage_selec$regulatoryGene),"sST2")

network_cortico_nousage_scaled_to_plot<-network_cortico_nousage_scaled[selected_interactors_cortico,selected_interactors_cortico]
network_cortico_nousage_scaled_to_plot[which(network_cortico_nousage_scaled_to_plot=="-Inf")]<-0


p_network_cortico_nousage_scaled<-ggballoonplot(network_cortico_nousage_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: 161 samples with steroid treatment timeline data",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)


 col2 = colorRampPalette(c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF"))

 #################
 #################
 #################

 dum<-df_v112_and_biomarkers_selec[c(2:52,54,56),ind_col_cortico_ivig_dates]
 gene_names<-rownames(dum)

 dum<-as.data.frame(dum)
 dum<-as.matrix(sapply(dum, as.numeric))

 df_v112_and_biomarkers_selec_cortico_ivig_nousage_scaled<-scale(dum, center = FALSE, scale=TRUE)
 rownames(df_v112_and_biomarkers_selec_cortico_ivig_nousage_scaled)<-gene_names

 set.seed(123)
#time interval,Glucocorticoid interval,IVIG interval
network_cortico_ivig_nousage_scaled<- GENIE3(df_v112_and_biomarkers_selec_cortico_ivig_nousage_scaled,nTrees = 1000)
links_network_cortico_ivig_nousage_scaled<-getLinkList(network_cortico_ivig_nousage_scaled)

correl_cortico_ivig_nousage_scaled<-corr.test(t(df_v112_and_biomarkers_selec_cortico_ivig_nousage_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)
correl_cortico_ivig_nousage_scaled<-as.data.frame(correl_cortico_ivig_nousage_scaled[["r"]])

correl_cortico_ivig_nousage_scaled_sign<-apply(correl_cortico_ivig_nousage_scaled,c(1,2), func_sign)

links_network_cortico_ivig_nousage_scaled_target_usage<-links_network_cortico_ivig_nousage_scaled[which(links_network_cortico_ivig_nousage_scaled$targetGene=="sTNFRI"),]

thres_0p5<-quantile(links_network_cortico_ivig_nousage_scaled_target_usage$weight, probs = c(0.00))

weight_thres<-as.data.frame(thres_0p5)$thres_0p5

links_network_cortico_ivig_nousage_scaled_target_usage_selec<-links_network_cortico_ivig_nousage_scaled_target_usage[which(links_network_cortico_ivig_nousage_scaled_target_usage$weight>=weight_thres),]

selected_interactors_cortico_ivig_nousage<-c(as.character(links_network_cortico_ivig_nousage_scaled_target_usage_selec$regulatoryGene),"sTNFRI")

network_cortico_ivig_nousage_scaled_to_plot<-network_cortico_ivig_nousage_scaled[selected_interactors_cortico_ivig_nousage,selected_interactors_cortico_ivig_nousage]
network_cortico_ivig_nousage_scaled_to_plot[which(network_cortico_ivig_nousage_scaled_to_plot=="-Inf")]<-0

pdf(file = "Figure_S4b.pdf", width=11, height=12)
corrplot(as.matrix(correl_cortico_ivig_nousage_scaled[selected_interactors_cortico_ivig_nousage,selected_interactors_cortico_ivig_nousage]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficients (101 samples with steroid and IVIG treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)
dev.off()

FigS4b<-corrplot(as.matrix(correl_cortico_ivig_nousage_scaled[selected_interactors_cortico_ivig_nousage,selected_interactors_cortico_ivig_nousage]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficients (103 samples with steroid and IVIG treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)
FigS4b<-as.data.frame(FigS4b[["corrPos"]])

################################
######CORTICO+NO IVIG#######
################################
dum<-df_v112_and_biomarkers_selec[c(2:52,54),ind_col_cortico_no_ivig_dates]
gene_names<-rownames(dum)
dum<-as.data.frame(dum)

dum<-as.matrix(sapply(dum, as.numeric))


df_v112_and_biomarkers_selec_cortico_no_ivig_nousage_scaled<-scale(dum, center = FALSE, scale=TRUE)
rownames(df_v112_and_biomarkers_selec_cortico_no_ivig_nousage_scaled)<-gene_names


set.seed(123)
#time interval,Glucocorticoid interval
network_cortico_no_ivig_nousage_scaled<- GENIE3(df_v112_and_biomarkers_selec_cortico_no_ivig_nousage_scaled,nTrees = 1000)
links_network_cortico_no_ivig_nousage_scaled<-getLinkList(network_cortico_no_ivig_nousage_scaled)


correl_cortico_no_ivig_nousage_scaled<-corr.test(t(df_v112_and_biomarkers_selec_cortico_no_ivig_nousage_scaled),  use = "pairwise",method="pearson",adjust="none", alpha=.05,ci=TRUE)
correl_cortico_no_ivig_nousage_scaled<-as.data.frame(correl_cortico_no_ivig_nousage_scaled[["r"]])


correl_cortico_no_ivig_nousage_scaled_sign<-apply(correl_cortico_no_ivig_nousage_scaled,c(1,2), func_sign)

links_network_cortico_no_ivig_nousage_scaled_target_usage<-links_network_cortico_no_ivig_nousage_scaled[which(links_network_cortico_no_ivig_nousage_scaled$targetGene=="sST2"),]

thres_0p5<-quantile(links_network_cortico_no_ivig_nousage_scaled_target_usage$weight, probs = c(0.0))

weight_thres<-as.data.frame(thres_0p5)$thres_0p5

links_network_cortico_no_ivig_nousage_scaled_target_usage_selec<-links_network_cortico_no_ivig_nousage_scaled_target_usage[which(links_network_cortico_no_ivig_nousage_scaled_target_usage$weight>=weight_thres),]

selected_interactors_cortico_no_ivig_nousage<-c(as.character(links_network_cortico_no_ivig_nousage_scaled_target_usage_selec$regulatoryGene),"sST2")

network_cortico_no_ivig_nousage_scaled_to_plot<-network_cortico_no_ivig_nousage_scaled[selected_interactors_cortico_no_ivig_nousage,selected_interactors_cortico_no_ivig_nousage]
network_cortico_no_ivig_nousage_scaled_to_plot[which(network_cortico_no_ivig_nousage_scaled_to_plot=="-Inf")]<-0


pdf(file = "Figure_S4a.pdf", width=11, height=12)
corrplot(as.matrix(correl_cortico_no_ivig_nousage_scaled[selected_interactors_cortico_no_ivig_nousage,selected_interactors_cortico_no_ivig_nousage]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficients (57 samples with steroid (but no IVIG) treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)
dev.off()


FigS4a<-corrplot(as.matrix(correl_cortico_no_ivig_nousage_scaled[selected_interactors_cortico_no_ivig_nousage,selected_interactors_cortico_no_ivig_nousage]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficients (59 samples with steroid (but no IVIG) treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)
FigS4a<-as.data.frame(FigS4a[["corrPos"]])



p_network_cortico_ivig_nousage_scaled<-ggballoonplot(network_cortico_ivig_nousage_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: Samples with steroid+IVIG treatment timeline data",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)

data_network_cortico_ivig_nousage_scaled<-ggpar(p_network_cortico_ivig_nousage_scaled,font.xtickslab= c(10, "plain", "black"),font.ytickslab= c(10, "plain", "black"), orientation = "horiz",font.main = c(10,"bold"),font.legend = c(8))
data_network_cortico_ivig_nousage_scaled<-data_network_cortico_ivig_nousage_scaled[["data"]]
colnames(data_network_cortico_ivig_nousage_scaled)[c(1,2)]<-c("predictor","target")


data_network_cortico_ivig_nousage_scaled<-data_network_cortico_ivig_nousage_scaled[order(-data_network_cortico_ivig_nousage_scaled$value),]


p_network_cortico_no_ivig_nousage_scaled<-ggballoonplot(network_cortico_no_ivig_nousage_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: Samples with steroid (but no IVIG) treatment timeline data",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)

data_network_cortico_no_ivig_nousage_scaled<-ggpar(p_network_cortico_no_ivig_nousage_scaled,font.xtickslab= c(10, "plain", "black"),font.ytickslab= c(10, "plain", "black"), orientation = "horiz",font.main = c(10,"bold"),font.legend = c(8))
data_network_cortico_no_ivig_nousage_scaled<-data_network_cortico_no_ivig_nousage_scaled[["data"]]
colnames(data_network_cortico_no_ivig_nousage_scaled)[c(1,2)]<-c("predictor","target")

data_network_cortico_no_ivig_nousage_scaled<-data_network_cortico_no_ivig_nousage_scaled[order(-data_network_cortico_no_ivig_nousage_scaled$value),]


p_network_cortico_scaled<-ggballoonplot(network_cortico_scaled_to_plot, fill = "value",xlab="Targets",ylab="Predictors",main = "Pairwise interaction strengths derived from random forest regression (x-axis:predictors & y-axis:targets) \nInput data: 92 samples with steroid treatment timeline data",size.range = c(1,5))+scale_fill_gradientn(colors = my_cols)

data_network_cortico_scaled_to_plot<-ggpar(p_network_cortico_scaled,font.xtickslab= c(10, "plain", "black"),font.ytickslab= c(10, "plain", "black"), orientation = "horiz",font.main = c(10,"bold"),font.legend = c(8))
data_network_cortico_scaled_to_plot<-data_network_cortico_scaled_to_plot[["data"]]
colnames(data_network_cortico_scaled_to_plot)[c(1,2)]<-c("predictor","target")


data_network_cortico_scaled_to_plot<-data_network_cortico_scaled_to_plot[order(-data_network_cortico_scaled_to_plot$value),]

colnames(data_network_cortico_ivig_nousage_scaled)[3]<-"interaction_strength"
ranked_predictors_network_cortico_ivig_nousage_scaled<-data_network_cortico_ivig_nousage_scaled %>% group_by(predictor)  %>% mutate(sum_predictor_strength=sum(interaction_strength))
ranked_predictors_network_cortico_ivig_nousage_scaled<-unique(ranked_predictors_network_cortico_ivig_nousage_scaled[,c(1,4)])
ranked_predictors_network_cortico_ivig_nousage_scaled<-ranked_predictors_network_cortico_ivig_nousage_scaled[order(-ranked_predictors_network_cortico_ivig_nousage_scaled$sum_predictor_strength),]


colnames(data_network_cortico_no_ivig_nousage_scaled)[3]<-"interaction_strength"
ranked_predictors_network_cortico_no_ivig_nousage_scaled<-data_network_cortico_no_ivig_nousage_scaled %>% group_by(predictor)  %>% mutate(sum_predictor_strength=sum(interaction_strength))
ranked_predictors_network_cortico_no_ivig_nousage_scaled<-unique(ranked_predictors_network_cortico_no_ivig_nousage_scaled[,c(1,4)])
ranked_predictors_network_cortico_no_ivig_nousage_scaled<-ranked_predictors_network_cortico_no_ivig_nousage_scaled[order(-ranked_predictors_network_cortico_no_ivig_nousage_scaled$sum_predictor_strength),]

colnames(data_network_cortico_scaled_to_plot)[3]<-"interaction_strength"
ranked_predictors_network_cortico_scaled_to_plot<-data_network_cortico_scaled_to_plot %>% group_by(predictor)  %>% mutate(sum_predictor_strength=sum(interaction_strength)) ranked_predictors_network_cortico_scaled_to_plot<-unique(ranked_predictors_network_cortico_scaled_to_plot[,c(1,4)])
ranked_predictors_network_cortico_scaled_to_plot<-ranked_predictors_network_cortico_scaled_to_plot[order(-ranked_predictors_network_cortico_scaled_to_plot$sum_predictor_strength),]

####this is the data with the usage included
#links_network_cortico_ivig_scaled
#links_network_cortico_no_ivig_scaled

colnames(links_network_cortico_no_ivig_scaled)[c(1,3)]<-c("predictor","interaction_strength")

ranked_predictors_network_cortico_no_ivig_with_usage_scaled<-links_network_cortico_no_ivig_scaled %>% group_by(predictor)  %>% mutate(sum_predictor_strength=sum(interaction_strength))

ranked_predictors_network_cortico_no_ivig_with_usage_scaled<-unique(ranked_predictors_network_cortico_no_ivig_with_usage_scaled[,c(1,4)])

ranked_predictors_network_cortico_no_ivig_with_usage_scaled<-ranked_predictors_network_cortico_no_ivig_with_usage_scaled[order(-ranked_predictors_network_cortico_no_ivig_with_usage_scaled$sum_predictor_strength),]


colnames(links_network_cortico_ivig_scaled)[c(1,3)]<-c("predictor","interaction_strength")

ranked_predictors_network_cortico_ivig_with_usage_scaled<-links_network_cortico_ivig_scaled %>% group_by(predictor)  %>% mutate(sum_predictor_strength=sum(interaction_strength))

ranked_predictors_network_cortico_ivig_with_usage_scaled<-unique(ranked_predictors_network_cortico_ivig_with_usage_scaled[,c(1,4)])

ranked_predictors_network_cortico_ivig_with_usage_scaled<-ranked_predictors_network_cortico_ivig_with_usage_scaled[order(-ranked_predictors_network_cortico_ivig_with_usage_scaled$sum_predictor_strength),]

df1<-ranked_predictors_network_cortico_no_ivig_nousage_scaled #52
df2<-ranked_predictors_network_cortico_no_ivig_with_usage_scaled #53
df3<-ranked_predictors_network_cortico_ivig_nousage_scaled #53
df4<-ranked_predictors_network_cortico_ivig_with_usage_scaled #54

setdiff(df2$predictor,df1$predictor) #"TRBV11-2 usage"
setdiff(df4$predictor,df3$predictor) #"TRBV11-2 usage"
setdiff(df4$predictor,df2$predictor) #"IVIG interval"
setdiff(df4$predictor,df1$predictor) #"IVIG interval"  "TRBV11-2 usage"

df1_dum<-data.frame(predictor=c("TRBV11-2 usage","IVIG interval"),sum_predictor_strength=c(0,0))
df1<-rbind(df1,df1_dum)

df2_dum<-data.frame(predictor=c("IVIG interval"),sum_predictor_strength=c(0))
df2<-rbind(df2,df2_dum)

df3_dum<-data.frame(predictor=c("TRBV11-2 usage"),sum_predictor_strength=c(0))
df3<-rbind(df3,df3_dum)

colnames(df1)[2]<-"Steroid_no_IVIG"
colnames(df2)[2]<-"Steroid_no_IVIG_TCR"
colnames(df3)[2]<-"Steroid_IVIG"
colnames(df4)[2]<-"Steroid_IVIG_TCR"

####combining all data frames by their "predictor" column
#######################
#######################
df_summed_predictor_strength<-merge(df1,df2,by="predictor")
df_summed_predictor_strength<-merge(df_summed_predictor_strength,df3,by="predictor")
df_summed_predictor_strength<-merge(df_summed_predictor_strength,df4,by="predictor")
#######################
#######################

row_variable_names<-df_summed_predictor_strength$predictor
df_summed_predictor_strength<-df_summed_predictor_strength[,-c(1)]
df_summed_predictor_strength<-as.matrix(df_summed_predictor_strength)
rownames(df_summed_predictor_strength)<-row_variable_names

colnames(df_summed_predictor_strength)<-c("Steroid only","Steroid only (TCR)","Steroid+IVIG","Steroid+IVIG (TCR)")

#####data frame with usage importance > 0, derived from samples with TCR data
df_summed_predictor_strength_usage<-df_summed_predictor_strength[,c(2,4)]

######only selecting the two larger groups (among four) regardless of the TCR data
df_summed_predictor_strength<-df_summed_predictor_strength[,c(1,3)]

df_summed_predictor_strength_annot<-as.data.frame(matrix(nrow=length(colnames(df_summed_predictor_strength)),ncol=1))
colnames(df_summed_predictor_strength_annot)<-c("Condition")

rownames(df_summed_predictor_strength_annot)<-colnames(df_summed_predictor_strength)

df_summed_predictor_strength_annot$Condition<-c("Steroid only","Steroid+IVIG")

ann_colors = list(Condition=c(`Steroid only`="burlywood",`Steroid+IVIG`="gold"))

table_summed_predictor_strength<-as.data.frame(df_summed_predictor_strength)
table_summed_predictor_strength$predictor<-rownames(df_summed_predictor_strength)

###########################
###########################
###########################

ranked_median_predictors_network_cortico_ivig_nousage_scaled<-data_network_cortico_ivig_nousage_scaled %>% group_by(predictor)  %>% mutate(median_predictor_strength=median(interaction_strength))
ranked_median_predictors_network_cortico_ivig_nousage_scaled<-unique(ranked_median_predictors_network_cortico_ivig_nousage_scaled[,c(1,4)])
ranked_median_predictors_network_cortico_ivig_nousage_scaled<-ranked_median_predictors_network_cortico_ivig_nousage_scaled[order(-ranked_median_predictors_network_cortico_ivig_nousage_scaled$median_predictor_strength),]


ranked_median_predictors_network_cortico_no_ivig_nousage_scaled<-data_network_cortico_no_ivig_nousage_scaled %>% group_by(predictor)  %>% mutate(median_predictor_strength=median(interaction_strength))
ranked_median_predictors_network_cortico_no_ivig_nousage_scaled<-unique(ranked_median_predictors_network_cortico_no_ivig_nousage_scaled[,c(1,4)])
ranked_median_predictors_network_cortico_no_ivig_nousage_scaled<-ranked_median_predictors_network_cortico_no_ivig_nousage_scaled[order(-ranked_median_predictors_network_cortico_no_ivig_nousage_scaled$median_predictor_strength),]


ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled<-links_network_cortico_no_ivig_scaled %>% group_by(predictor)  %>% mutate(median_predictor_strength=median(interaction_strength))
ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled<-unique(ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled[,c(1,4)])
ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled<-ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled[order(-ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled$median_predictor_strength),]

ranked_median_predictors_network_cortico_ivig_with_usage_scaled<-links_network_cortico_ivig_scaled %>% group_by(predictor)  %>% mutate(median_predictor_strength=median(interaction_strength))
ranked_median_predictors_network_cortico_ivig_with_usage_scaled<-unique(ranked_median_predictors_network_cortico_ivig_with_usage_scaled[,c(1,4)])
ranked_median_predictors_network_cortico_ivig_with_usage_scaled<-ranked_median_predictors_network_cortico_ivig_with_usage_scaled[order(-ranked_median_predictors_network_cortico_ivig_with_usage_scaled$median_predictor_strength),]


df1_median<-ranked_median_predictors_network_cortico_no_ivig_nousage_scaled #52
df2_median<-ranked_median_predictors_network_cortico_no_ivig_with_usage_scaled #53
df3_median<-ranked_median_predictors_network_cortico_ivig_nousage_scaled #53
df4_median<-ranked_median_predictors_network_cortico_ivig_with_usage_scaled #54


df1_median_dum<-data.frame(predictor=c("TRBV11-2 usage","IVIG interval"),median_predictor_strength=c(0,0))
df1_median<-rbind(df1_median,df1_median_dum)

df2_median_dum<-data.frame(predictor=c("IVIG interval"),median_predictor_strength=c(0))
df2_median<-rbind(df2_median,df2_median_dum)

df3_median_dum<-data.frame(predictor=c("TRBV11-2 usage"),median_predictor_strength=c(0))
df3_median<-rbind(df3_median,df3_median_dum)

colnames(df1_median)[2]<-"Steroid_no_IVIG"
colnames(df2_median)[2]<-"Steroid_no_IVIG_TCR"
colnames(df3_median)[2]<-"Steroid_IVIG"
colnames(df4_median)[2]<-"Steroid_IVIG_TCR"

####combining all data frames by their "predictor" column
#######################
#######################
df_median_predictor_strength<-merge(df1_median,df2_median,by="predictor")
df_median_predictor_strength<-merge(df_median_predictor_strength,df3_median,by="predictor")
df_median_predictor_strength<-merge(df_median_predictor_strength,df4_median,by="predictor")
#######################
#######################

row_variable_names<-df_median_predictor_strength$predictor
df_median_predictor_strength<-df_median_predictor_strength[,-c(1)]
df_median_predictor_strength<-as.matrix(df_median_predictor_strength)
rownames(df_median_predictor_strength)<-row_variable_names

colnames(df_median_predictor_strength)<-c("Steroid only","Steroid only (TCR)","Steroid+IVIG","Steroid+IVIG (TCR)")

######only selecting the two larger groups (among four) regardless of the TCR data
df_median_predictor_strength<-df_median_predictor_strength[,c(1,3)]

df_median_predictor_strength_annot<-as.data.frame(matrix(nrow=length(colnames(df_median_predictor_strength)),ncol=1))
colnames(df_median_predictor_strength_annot)<-c("Condition")

rownames(df_median_predictor_strength_annot)<-colnames(df_median_predictor_strength)

df_median_predictor_strength_annot$Condition<-c("Steroid only","Steroid+IVIG")

ann_colors = list(Condition=c(`Steroid only`="burlywood",`Steroid+IVIG`="gold"))

ind_remov<-match("TRBV11-2 usage",rownames(df_median_predictor_strength))
df_median_predictor_strength<-df_median_predictor_strength[-ind_remov,]

######this is the updated version with cluster_cols = FALSE
p_df_median_predictor_strength_updated<-pheatmap(mat=df_median_predictor_strength, show_colnames = F,cutree_rows = 4,main="Median predictive importance values \n derived from RF regression ",fontsize_row=8, annotation_col = df_median_predictor_strength_annot,annotation_colors = ann_colors, fontsize = 8,cluster_cols = FALSE, clustering_distance_rows = "euclidean",clustering_method = "average",scale="none")

pdf(file="Figure_3g.pdf",width=5,height=8)
print(p_df_median_predictor_strength_updated) # this is ComplexHeatmap::pheatmap
dev.off()

Fig_3g<-as.data.frame(df_median_predictor_strength)
Fig_3g$predictor<-rownames(Fig_3g)
Fig_3g<-Fig_3g[,c(3,1,2)]
rownames(Fig_3g)<-c(1:nrow(Fig_3g))
####predictor	Steroid only	Steroid+IVIG

ind_remov<-match("TRBV11-2 usage",rownames(df_median_predictor_strength))
df_median_predictor_strength<-df_median_predictor_strength[-ind_remov,]


Supp_4a<-corrplot(as.matrix(correl_cortico_no_ivig_nousage_scaled[selected_interactors_cortico_no_ivig_nousage,selected_interactors_cortico_no_ivig_nousage]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficients (57 samples with steroid (but no IVIG) treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)

Supp_4a<-as.data.frame(Supp_4a[["corrPos"]])

Supp_4b<-corrplot(as.matrix(correl_cortico_ivig_nousage_scaled[selected_interactors_cortico_ivig_nousage,selected_interactors_cortico_ivig_nousage]),col = col2(200),diag=FALSE,method="circle",tl.cex=0.9,tl.col="black",cl.ratio=0.2,cl.cex=0.9,number.cex=0.9, title = "Pearson correlation coefficients (101 samples with steroid and IVIG treatment timeline data)", mar=c(0,0,2,0), order = 'hclust', addrect = 2)

Supp_4b<-as.data.frame(Supp_4b[["corrPos"]])
