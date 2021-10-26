### functions for parsing and ploting fgsea results
### utility funs #############################################################################
# merge GSEA list into dataframe 
RbindGseaResultList = function(gsea_result_list, NES_filter = 0, padj_filter = 0.05){
  
  score = lapply(gsea_result_list, function(x){ 
    x = x %>% 
      dplyr::select(pathway, padj, NES) %>% 
      dplyr::filter( abs(NES) > NES_filter ) 
  })
  score = dplyr::bind_rows(score, .id = "celltype") %>% filter(padj < padj_filter) %>% mutate(n_logp = -log10(padj))
  return(score)
}


GSEABubblePlot <- function(fgseaRes.df){
  p = ggplot(fgseaRes.df, aes(y = pathway, x = celltype, fill = NES, size = n_logp)) + 
    geom_point(shape = 21) + 
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",midpoint = 0) +
    theme_bw() +
    scale_x_discrete(position = "top") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    theme(axis.title.y = element_blank()) +
    # theme(legend.position = "bottom") + 
    labs(fill = 'Normalized Enrichment Score', size = '-log10(adjusted P value)') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    guides(shape = guide_legend(override.aes = list(size = 5))) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))
  return(p)
}

GSEABubblePlot_selected <- function(fgseaRes.df){
  p = ggplot(fgseaRes.df, aes(y = pathway, x = celltype, fill = NES, size = n_logp)) + 
    geom_point(shape = 21) + 
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red3",midpoint = 0) +
    theme_bw() +
    scale_x_discrete(position = "top") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 0)) + 
    theme(axis.title.y = element_blank()) +
    # theme(legend.position = "bottom") + 
    labs(fill = 'Normalized Enrichment Score', size = '-log10(adjusted P value)') +
    theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
    theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
    theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8), panel.spacing = unit(0, "lines"))+
    guides(shape = guide_legend(override.aes = list(size = 5))) + 
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    facet_grid(Category~ cellgroup, scales = "free", space = "free") +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))+
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank()
    )
  
  return(p)
}













