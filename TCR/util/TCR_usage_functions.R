
##############################
### utility functions TCR ####
##############################
# calculate V-gene usage ratio of each subject
cal_Vgene_subject_ratio <- function(df, chain){
  #chain = c("TRAV","TRAJ","TRBV","TRBJ","TRBC","IGHV","IGHJ","IGHC","IGLV","IGLJ","IGLCons")
  ratio_df <- data.frame(table(df[,chain], df[,"sample_id"])) %>% 
    dplyr::group_by(Var2) %>%
    dplyr::mutate(ratio = Freq/sum(Freq)) %>%
    drop_na() %>%
    # dcast(Var2 ~ Var1, value.var = "ratio") %>%
    left_join(meta_subject, by = c("Var2" = "sample_id"))
  return(ratio_df)
}


# plot V-gene usage ratio of different class
plot_Vgene_usage_ratio <- function(df, celltype, chain){
  my_comparisons <- list( c("HC", "COVID"), c("HC", "MIS-C"), c("COVID", "MIS-C") )
  p <- ggplot(df, aes(x = Class, y = ratio, color = Class)) +
    geom_boxplot()+
    geom_point(position = position_jitterdodge())+
    facet_wrap(~Var1, scales = "free")+
    stat_compare_means(comparisons = my_comparisons, mmethod = "wilcox.test", size=2)+
    # ggsci::scale_color_npg(name="Class",guide=guide_legend(title.position = "top"))+
    scale_color_manual(values = c("HC" = "#0000AC", "COVID" = "#00AFBB", "MIS-C" = "#E7B800"))+
    theme_bw()+
    ggtitle(paste(celltype, chain, sep = "-"))
  return(p)
}

