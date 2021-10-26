
library(tidyverse)
library(magrittr)
# get metadata
md = B1merge@meta.data %>% select(adt_snn_res.2)


# get data
adt = GetAssayData(B1merge[["CITE"]]) %>%
  t %>%
  as.data.frame %>%
  rownames_to_column("cell")
md %<>% rownames_to_column("cell")
# idont know if you need to do this.
md = md[match(x = md$cell, table = adt$cell),  ]
adt %<>%
  select(-c(cell)) %>%
  mutate(cell_type = md$adt_snn_res.2) %>%
  select(cell_type, everything())


# get hclust order of proteins and cell types
    x = pheatmap::pheatmap(GetAssayData(aver[["CITE"]]),
             cluster_rows = T, cluster_cols = T,
             fontsize_col = 10, fontsize_row = 8, border_color = NA)

    prot_order = rownames(GetAssayData(aver[["CITE"]])[x$tree_row$order, ])
    celltype_order = colnames(GetAssayData(aver[["CITE"]])[,x$tree_col$order ])

library('genefilter')
#filter out proteins that do not have at least 1 cluster with an average over 2
    f1 <- kOverA(1, 2)
    ffun <- filterfun(f1)
    filtaver.2 <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, ], ffun)

    prot_plot = prot_order[which(filtaver.2)]


adt.l = adt %>%
  gather(key = prot, value = normalized_count, CD80.1:DR3) %>%
  filter(prot %in% prot_plot)

library('DescTools')
  adt.l$prot = factor(adt.l$prot)
  adt.l$prot =  reorder.factor(adt.l$prot, new.order = prot_plot)
  adt.l$cell_type = factor(adt.l$cell_type)
  adt.l$cell_type =  reorder.factor(adt.l$cell_type, new.order = celltype_order)


suppressMessages(library(ggridges))
library(viridis)
library(scico)
p =  ggplot(adt.l, aes(x=normalized_count, y=prot, fill = ..x..)) +
   geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
     scale_fill_viridis(name = "normalized count", option = "B") +
    scale_x_continuous(limits = c(-1,25)) +
    theme_ridges(font_size = 8, line_size = 0.01) +
  facet_grid(~cell_type)

pdf("res.2.protHistogram.pdf", 24, 16)
print(p)
dev.off()


filtaver.2.CD4s <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, which(GetAssayData(aver[["CITE"]])["CD3",] > 3 & GetAssayData(aver[["CITE"]])["CD4.1",] > 3) ], ffun)

adt.l.CD4s = adt.l[which(adt.l$cell_type %in% names(which(GetAssayData(aver[["CITE"]])["CD3",] > 3 & GetAssayData(aver[["CITE"]])["CD4.1",] > 3))),] %>% filter(prot %in% prot_order[which(filtaver.2.CD4s)])

p.CD4s =  ggplot(adt.l.CD4s, aes(x=normalized_count, y=prot, fill = ..x..)) +
   geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
     scale_fill_viridis(name = "normalized count", option = "B") +
    scale_x_continuous(limits = c(-1,15)) +
    theme_ridges(font_size = 8, line_size = 0.1) +
  facet_grid(~cell_type)

  pdf("res.2.CD4s.protHistogram.pdf", 10, 7)
  print(p.CD4s)
  dev.off()

#CD8 hist
filtaver.2.CD8s <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, which(GetAssayData(aver[["CITE"]])["CD3",] > 3 & GetAssayData(aver[["CITE"]])["CD8",] > 3) ], ffun)

  adt.l.CD8s = adt.l[which(adt.l$cell_type %in% names(which(GetAssayData(aver[["CITE"]])["CD3",] > 3 & GetAssayData(aver[["CITE"]])["CD8",] > 3))),] %>% filter(prot %in% prot_order[which(filtaver.2.CD8s)])

  p.CD8s =  ggplot(adt.l.CD8s, aes(x=normalized_count, y=prot, fill = ..x..)) +
     geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
       scale_fill_viridis(name = "normalized count", option = "B") +
      scale_x_continuous(limits = c(-1,15)) +
      theme_ridges(font_size = 8, line_size = 0.1) +
    facet_grid(~cell_type)

    pdf("res.2.CD8s.protHistogram.pdf", 10, 7)
    print(p.CD8s)
    dev.off()

#Mono Hist
#FeatureScatter(aver, feature1="cite_CD14", feature2="cite_CD64") #CD64 better for mono marker, to get the non-classicals
    filtaver.2.Monos <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, which(GetAssayData(aver[["CITE"]])["CD64",] > 3) ], ffun)

    adt.l.Monos = adt.l[which(adt.l$cell_type %in% names(which(GetAssayData(aver[["CITE"]])["CD64",] > 3))),] %>% filter(prot %in% prot_order[which(filtaver.2.Monos)])

    p.Monos =  ggplot(adt.l.Monos, aes(x=normalized_count, y=prot, fill = ..x..)) +
       geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
         scale_fill_viridis(name = "normalized count", option = "B") +
        scale_x_continuous(limits = c(-1,15)) +
        theme_ridges(font_size = 8, line_size = 0.1) +
      facet_grid(~cell_type)

      pdf("res.2.Monos.protHistogram.pdf", 10, 7)
      print(p.Monos)
      dev.off()

#NK Hist
#FeatureScatter(merged, feature1="cite_CD56", feature2="cite_CD16")
          filtaver.2.NK <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, which(GetAssayData(aver[["CITE"]])["CD56",] > 3) ], ffun)

          adt.l.NK = adt.l[which(adt.l$cell_type %in% names(which(GetAssayData(aver[["CITE"]])["CD56",] > 3))),] %>% filter(prot %in% prot_order[which(filtaver.2.NK)])

          p.NK =  ggplot(adt.l.NK, aes(x=normalized_count, y=prot, fill = ..x..)) +
             geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
               scale_fill_viridis(name = "normalized count", option = "B") +
              scale_x_continuous(limits = c(-1,15)) +
              theme_ridges(font_size = 8, line_size = 0.1) +
            facet_grid(~cell_type)

            pdf("res.2.NKs.protHistogram.pdf", 10, 5)
            print(p.NK)
            dev.off()

#B Hist
#FeatureScatter(aver, feature1="cite_CD20", feature2="cite_CD19")
                      filtaver.2.B <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, which(GetAssayData(aver[["CITE"]])["CD19.1",] > 3) ], ffun)

                      adt.l.B = adt.l[which(adt.l$cell_type %in% names(which(GetAssayData(aver[["CITE"]])["CD19.1",] > 3))),] %>% filter(prot %in% prot_order[which(filtaver.2.B)])

                      p.B =  ggplot(adt.l.B, aes(x=normalized_count, y=prot, fill = ..x..)) +
                         geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
                           scale_fill_viridis(name = "normalized count", option = "B") +
                          scale_x_continuous(limits = c(-1,15)) +
                          theme_ridges(font_size = 8, line_size = 0.1) +
                        facet_grid(~cell_type)

                        pdf("res.2.Bcells.protHistogram.pdf", 10, 7)
                        print(p.B)
                        dev.off()


#UnConvT Hist
          filtaver.2.UnConvT <- genefilter(GetAssayData(aver[["CITE"]])[x$tree_row$order, which(GetAssayData(aver[["CITE"]])["CD3",] > 3 & GetAssayData(aver[["CITE"]])["CD4.1",] < 3 & GetAssayData(aver[["CITE"]])["CD8",] < 3) ], ffun)

          adt.l.UnConvT = adt.l[which(adt.l$cell_type %in% names(which(GetAssayData(aver[["CITE"]])["CD3",] > 3 & GetAssayData(aver[["CITE"]])["CD4.1",] < 3 & GetAssayData(aver[["CITE"]])["CD8",] < 3))),] %>% filter(prot %in% prot_order[which(filtaver.2.UnConvT)])

          p.UnConvT =  ggplot(adt.l.UnConvT, aes(x=normalized_count, y=prot, fill = ..x..)) +
                  geom_density_ridges_gradient(alpha = 0.7, scale = 3, rel_min_height = 0.02) +
                    scale_fill_viridis(name = "normalized count", option = "B") +
                      scale_x_continuous(limits = c(-1,15)) +
                        theme_ridges(font_size = 8, line_size = 0.1) +
                          facet_grid(~cell_type)

          pdf("res.2.DNT.protHistogram.pdf", 10, 5)
          print(p.UnConvT)
          dev.off()
