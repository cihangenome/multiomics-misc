####################
### utility funs ###
####################
# plot time effect LE genes
# use:
# subject.hm(exprs_mtx = Mono_Classical_mtx, meta = meta, celltype = "Mono_Classical", module = "IFN")
subject.hm <- function(exprs_mtx, meta, celltype, module){
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#1F9F88FF", "#000000", "#FFFF00")) # green to yellow
  timeonset.color <- colorRamp2(c(0, 40), c("white", "#6a2c70"))
  Class.color <- c("HC" = "#0000AC", "COVID" = "#00AFBB", "MIS-C" = "#E7B800")

  # set how many genes to label, due to space limit
  if(nrow(exprs_mtx) < 20){
    nlabel = 10
  }else if (nrow(exprs_mtx) >= 20 & nrow(exprs_mtx) < 60){
    nlabel = 15
  }else {nlabel = round(nrow(exprs_mtx)/2.5)}
  label <- c(1:nlabel, which(rownames(exprs_mtx) %in% c("IL15RA","IL1B","ICOSL")))
  label_genes <- rownames(exprs_mtx)[label]
  # set heatmap annotation
  ha = HeatmapAnnotation(
    Class = meta$Class,
    Time = meta$days_since_admission,
    col = list(Class = Class.color,
               Time = timeonset.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"),
               show_column_names = FALSE,
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$Class,
               col = col_fun,
               column_gap = unit(2, "mm")
  ) +
    rowAnnotation(foo = anno_mark(at = label,
                                  labels = label_genes,
                                  labels_gp = gpar(fontsize = 12)))
  return(p)
}


### function for showing all genes and also highlight intersect
subject.hm2 <- function(exprs_mtx, meta, celltype, module, LEset1, LEset2){
  # set marked genesets
  coLE <- rownames(exprs_mtx) %in% intersect(LEset1, LEset2)
  names(coLE) <- rownames(exprs_mtx)
  
  # set colors
  col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#1F9F88FF", "#000000", "#FFFF00")) # green to yellow
  timeonset.color <- colorRamp2(c(0, 40), c("white", "#6a2c70"))
  Class.color <- c("HC" = "#0000AC", "COVID" = "#00AFBB", "MIS-C" = "#E7B800")

  # set how many genes to label, due to space limit
  label <- c(1:2, which(rownames(exprs_mtx) %in% intersect(LEset1, LEset2)))
  label_genes <- rownames(exprs_mtx)[label]
  # set heatmap annotation
  ha = HeatmapAnnotation(
    Class = meta$Class,
    Time = meta$days_since_admission,
    col = list(Class = Class.color,
               Time = timeonset.color)
  )
  
  # plot
  p <- Heatmap(pheatmap:::scale_rows(exprs_mtx), name = paste(celltype, module, sep = "_"),
               show_column_names = FALSE,
               top_annotation = ha, cluster_columns = FALSE,
               column_split = meta$Class,
               col = col_fun,
               # row_names_gp = gpar(fontsize = 8),
               column_gap = unit(2, "mm")
  )+
    Heatmap(coLE + 0, name = "common_LE", col = c("0" = "white", "1" = "#C6D57E"), 
            show_heatmap_legend = FALSE, width = unit(4, "mm")) +
    rowAnnotation(foo = anno_mark(at = label,
                                  labels = label_genes,
                                  labels_gp = gpar(fontsize = 11)))
  
  return(p)
}

