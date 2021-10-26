# multiomics-misc
This repository contains the R scripts that were used in the immune repertoire and gene expression analysis presented in "Multiomics approach identifies novel age-, time- and treatment-related immunopathological signatures in MIS-C and pediatric COVID-19". The single cell expression and repertoire analysis scripts have been contributed by Can Liu.


# Script descriptions

Input files are in folder -- input folder
Output files are default to -- output folder
Utility functions used in the codes are in folder -- util_fun folder
pesudobulk scripts -- pbulk_DE folder
BCR scripts -- BCR folder
TCR scripts -- TCR folder

## 1.demux.preprocess.R

	This script takes the processed data from cellranger and demuxlet files. Creats a Seurat object, performs QC filtering of the cells, hashtag demultiplexing and DSB normalization of the CITEseq data.
	
## 2.clustering.R

	This script takes the Seurat object output by 1.demux.preprocess.R, and performs clustering using the cell surface antibody data (ADT/CITE) and adds the cell annotations.
	
## 3.1.label.transfer.from.3batches.R
	
	Generates Figure S5A
	This script uses the data from adult COVID data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161918) and integrates and predicts celltype labels of this dataset from annotated adult COVID dataset.

## 3.2.cell_no_check_all4batches.R

	Generates Figure S5B
	Cell frequency summary taking taking the predicted labels and adult COVID dataset to compare interested cell populations in both adult and pediatric patients

## 4.meta.plot.R

  Generates Figure 4A and Figure S5C
  This script takes the Seurat object (uploaded) and plots the UMAP visulization of total cells and Monocytes. Also generates single cell CD163 expression of monocytes and S100A family genes expression plots.



# pesudobulk scripts are in pbulk_DE folder
	
	Main workflow for the pseudobulk differential expresssion analyses are as follows:


## 1.pbulk.pooling.unfiltered.R
	
	Pool reads from single cells together per celltype per sample

## 2.calcnormfactors_and_filter_genes_samples.R

	Caculate normalization factors, normalize pseudobulk libraries and filter genes and samples

## 3.make_design_matrix.R

	make design matrix for differential expression analyses

## 4.run_de_limma.R
	
	Run limma differential expression model

## 5.make_contrast_fit.R

	Make contrast matrix for each celltype and get top table of differentially expressed genes

## 6.run_fgsea.R
	
	Run fgsea for the pseudobulk differential expresssion results using the t-statistic ranks

## 7.plot_fgsea.R

	Generates Figures 4B, 4C, S5D
	Fgsea results visulization of selected pathways

## 8.DEplot.subject.hm.R

	Generates Figure 4D
	Heatmaps of HALLMARK_TNFa_Signaling_via_NFkB gene set leading edge genes. Showing CD4 Memory and Classical Monocytes as examples.

## 9.gsva.score.miscvscovid.R
	
	Generates Figure 4E
	Calculate GSVA scores using leading edge genes of MIS-C vs p.COVID-19. Take input data from normalized pseudobulk objects.
	Use GSVA scores to genereate HALLMARK_TNFa_Signaling_via_NFkB gene set signature score boxplot.



# BCR scripts are in BCR folder

## 1.add.BCR.meta

	Add BCR metadata to Seurat object.

## 2.BCRstat.wo.geneexpr.R

	Generates Figures 5G, S7B, S7D
	BCR IgH usage and B cell mutation frequency summary.

## 3.BCR.mut.corr.R

	Generates Figures S7E, S7F
	Memory B cell and Plasmablast mutation frequencies correlation with cell surface protein.



# TCR scripts are in TCR folder

## 1.compile.TCR.data.R
	
	Compile raw TCR data and add metadata to the cells.

## 2.TCRstat.wo.geneexpr.R

	Generates Figure 5C
	CD4 TCR TRBV11-2 usage summary.

## 3.add.TCR.meta.R

	Add TCR metadata to Seurat object.

## 4.TRBV112CD4.gex.R
	
	Generates Figures 5D,5E
	TRBV11-2 MIS-C CD4 T cell phenotype, compared to other TRBV11-2 negative CD4 T cells within MIS-C patients. Heatmaps of genes and surface proteins markers of MIS-C TRBV11-2 CD4 T cell.

## 5.TRBV112CD4.gsea.R

	Generates Figures 5F
	Fgsea test on the genes from single cell level differential expression analysis of TRBV11-2 CD4 T cells within MIS-C patients. Plot of apoptosis related gene sets enrichment.
