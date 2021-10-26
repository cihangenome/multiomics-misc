### untility funs ############################################################################
### untility funs for the cell frequency check ###############################################

### adding basic meta data
addcellmeta <- function(df){
  df$Subject <- sapply(str_split(df$Var1, pattern = "_"), function(x)x[1])
  df$Age <- plyr::mapvalues(x = df$Var1, 
                            from = samplemetadata$sample_id, to = samplemetadata$Age)
  df$Gender <- plyr::mapvalues(x = df$Var1, 
                               from = samplemetadata$sample_id, to = samplemetadata$Gender)
  df$Class <- plyr::mapvalues(x = df$Var1, 
                              from = samplemetadata$sample_id, to = samplemetadata$Class)
  df$Timepoint <- plyr::mapvalues(x = df$Var1, 
                                  from = samplemetadata$sample_id, to = samplemetadata$Timepoint)
  df$Timepoint2 <- plyr::mapvalues(x = df$Var1, 
                                   from = samplemetadata$sample_id, to = samplemetadata$Timepoint2)
  df$days_since_symptom_onset <- plyr::mapvalues(x = df$Var1, 
                                                 from = samplemetadata$sample_id, to = samplemetadata$days_since_symptom_onset)
  # df$B_Mem.shannon <- plyr::mapvalues(x = df$Var1, 
  #                                     from = samplemetadata$sample_id, to = samplemetadata$B_Mem.shannon)
  # df$B_Naive.shannon <- plyr::mapvalues(x = df$Var1, 
  #                                       from = samplemetadata$sample_id, to = samplemetadata$B_Naive.shannon)
  df$Class <- factor(df$Class, levels = c("HC","COVID","MIS-C"))
  # df$B_Mem.shannon <- as.numeric(as.character(df$B_Mem.shannon))
  # df$B_Naive.shannon <- as.numeric(as.character(df$B_Naive.shannon))
  return(df)
}



plotcelltype <- function(df, name, ratioto, width = 16, height = 10){
  
  # Class
  df_T0 <- filter(df, Timepoint %in% c("T0","HC"))
  p1 <- ggplot(df_T0, aes(x = Class, y = value, fill = Class))+
    geom_boxplot(outlier.shape=NA)+geom_point()+
    ggsci::scale_fill_npg(name="Class") +
    # theme(axis.text.x = element_text(angle = 25),
    # text = element_text(size=20))+
    facet_wrap(~variable, scales = "free")+
    theme_bw()
  ggsave(filename = paste(name, "class", ratioto, "pdf", sep = "."), 
         plot = p1, device = "pdf", width = width+6, height = height+6)
  
  # timecourse
  tmp <- filter(df, Timepoint != "HC")
  tmp.HC <- filter(df, Timepoint == "HC")
  means <- plyr::ddply(tmp.HC, .(variable), summarise, median = median(value, na.rm = TRUE), 
                       quantile.25 = quantile(value, .25, na.rm = TRUE), 
                       quantile.75 = quantile(value, .75, na.rm = TRUE))
  p2 <- ggplot(tmp)+
    # geom_rect(aes(xmin = 17, xmax = 23, ymin = -Inf, ymax = Inf), fill = "#EADCFA", alpha = 0.1)+
    geom_point(alpha=.4, shape=21,
               aes(x = days_since_symptom_onset, y = value, fill=Class), size = 1.2)+
    ggsci::scale_fill_npg(name="Class") +
    geom_line(aes(x = days_since_symptom_onset, y = value, group = Subject), alpha = 0.08)+
    stat_summary(aes(x = days_since_symptom_onset, y = value, group = Class), fun.y = median, alpha=0)+
    # scale_shape_manual(values = c(15:18))+
    # stat_smooth(aes(x = days_since_symptoms_onset, y = value, group = PC1class, color=PC1class), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, size = 0.5)+# group = 1 overwrite group id
    stat_smooth(aes(x = days_since_symptom_onset, y = value, group = Class, color=Class), se = FALSE, size = 0.7)+# group = 1 overwrite group id
    ggsci::scale_color_npg(name="Class")+
    facet_wrap(~variable, scales = "free")+
    # geom_rect(data = means, aes(xmin = -Inf, xmax = Inf, ymin = quantile.25, ymax = quantile.75), fill = "#00BA38", alpha = 0.1)+
    geom_hline(aes(yintercept = median), data = means, alpha = 0.7, linetype = "dashed", color = "#00BA38")+
    geom_vline(xintercept = 20, color = "grey60", linetype='dashed')+
    theme(axis.text.x = element_text(angle = 90))+
    theme_bw()+
    ggtitle("green line:HC median; shaded: HC 25-75% interval")
  ggsave(filename = paste(name, "samplevsonset", ratioto, "pdf", sep = "."),
         plot = p2, device = "pdf", width = width, height = height)
  
}











