
library(monocle3)

##B1
cds<-readRDS("5.monocle2/8.TAM_B1/TAM_B1_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p1<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 1)')


cluster_p1<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =RColorBrewer::brewer.pal(3,"Set1") )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 1)')

#挑选top10画图展示
Track_genes_sig <- c("MARCO")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"MARCO (Branch 1)"
plot_genes1<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))



Track_genes_sig <- c("FN1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"FN1 (Branch 1)"
plot_genes2<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


Track_genes_sig <- c("VSIG4")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"VSIG4 (Branch 1)"
plot_genes3<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))






##B2
cds<-readRDS("5.monocle2/9.TAM_B2//TAM_B2_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 2)')


cluster_p2<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =c(RColorBrewer::brewer.pal(3,"Set1")[1:2],RColorBrewer::brewer.pal(3,"Set2")[1]) )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 2)')


#挑选top10画图展示
Track_genes_sig <- c("MARCO")
cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GC","GCPM","PE")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"MARCO (Branch 2)"
plot_genes4<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))



Track_genes_sig <- c("FN1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"FN1 (Branch 2)"
plot_genes5<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


Track_genes_sig <- c("VSIG4")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"VSIG4 (Branch 2)"
plot_genes6<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


genes_for_branchus1<-plot_genes1+plot_genes4+plot_genes2+plot_genes5+plot_genes3+plot_genes6+plot_layout(ncol=2)
genes_for_branchus1







##B1
cds<-readRDS("5.monocle2/8.TAM_B1/TAM_B1_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p1<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 1)')


cluster_p1<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =RColorBrewer::brewer.pal(3,"Set1") )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 1)')

#挑选top10画图展示
Track_genes_sig <- c("CD9")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"CD9 (Branch 1)"
plot_genes1<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))



Track_genes_sig <- c("TREM2")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"TREM2 (Branch 1)"
plot_genes2<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


Track_genes_sig <- c("APOC1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"APOC1 (Branch 1)"
plot_genes3<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))






##B2
cds<-readRDS("5.monocle2/9.TAM_B2//TAM_B2_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 2)')


cluster_p2<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =c(RColorBrewer::brewer.pal(3,"Set1")[1:2],RColorBrewer::brewer.pal(3,"Set2")[1]) )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 2)')


#挑选top10画图展示
Track_genes_sig <- c("CD9")
cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GC","GCPM","PE")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"CD9 (Branch 2)"
plot_genes4<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))



Track_genes_sig <- c("TREM2")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"TREM2 (Branch 2)"
plot_genes5<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


Track_genes_sig <- c("APOC1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"APOC1 (Branch 2)"
plot_genes6<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


genes_for_branchus2<-plot_genes4+plot_genes5+plot_genes6+plot_layout(ncol=1)
genes_for_branchus2















##B1
cds<-readRDS("5.monocle2/8.TAM_B1/TAM_B1_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p1<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 1)')


cluster_p1<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =RColorBrewer::brewer.pal(3,"Set1") )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 1)')

#挑选top10画图展示
Track_genes_sig <- c("F13A1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"F13A1 (Branch 1)"
plot_genes1<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))



Track_genes_sig <- c("STAB1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"STAB1 (Branch 1)"
plot_genes2<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


Track_genes_sig <- c("FOLR2")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"FOLR2 (Branch 1)"
plot_genes3<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))






##B2
cds<-readRDS("5.monocle2/9.TAM_B2//TAM_B2_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 2)')


cluster_p2<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =c(RColorBrewer::brewer.pal(3,"Set1")[1:2],RColorBrewer::brewer.pal(3,"Set2")[1]) )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 2)')


#挑选top10画图展示
Track_genes_sig <- c("F13A1")
cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GC","GCPM","PE")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"F13A1 (Branch 2)"
plot_genes4<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))



Track_genes_sig <- c("STAB1")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"STAB1 (Branch 2)"
plot_genes5<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


Track_genes_sig <- c("FOLR2")
cds_subset <- cds[Track_genes_sig,]

#基因表达趋势图
rownames(cds_subset)<-"FOLR2 (Branch 2)"
plot_genes6<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.1)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


genes_for_branchus3<-plot_genes4+plot_genes5+plot_genes6+plot_layout(ncol=1)
genes_for_branchus3




##Feature_plot
sub_sce<-readRDS("3.Cluster/8.SubAnnotation/2.myeloid/sub_sce_annotation.rds")
Feature_P1<-FeaturePlot(sub_sce,features = c("MARCO","FN1") ,
                        ncol = 1,pt.size=0,cols = viridis::viridis(5,option = "D"),combine = F) 

plot_list<-list()
for(pt in Feature_P1){
  pt<-pt+theme(plot.title = element_text(size = 10,hjust = 0.5),
               axis.title.x = element_text(size = 10),
               axis.title.y = element_text(size = 10),
               axis.text.x = element_text(size = 8),
               axis.text.y = element_text(size = 10),
               legend.title = element_text(size = 10),
               
               legend.text = element_text(size = 10),
               legend.position = "bottom")
  
  plot_list[[length(plot_list)+1]]<-pt
}

Feature_P1<-ggarrange(plotlist = plot_list,common.legend = T,legend = "bottom",ncol=1)




Feature_P2<-FeaturePlot(sub_sce,features = c("F13A1","STAB1") ,
                        ncol = 1,pt.size=0,cols = viridis::viridis(5,option = "D"),combine = F) 

plot_list<-list()
for(pt in Feature_P2){
  pt<-pt+theme(plot.title = element_text(size = 10,hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 10),
          legend.position = "bottom")
  
  plot_list[[length(plot_list)+1]]<-pt
}

Feature_P2<-ggarrange(plotlist = plot_list,common.legend = T,legend = "bottom",ncol=1)







library(monocle3)

##B1
cds<-readRDS("5.monocle2/8.TAM_B1/TAM_B1_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p1<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 1)')


cluster_p1<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =RColorBrewer::brewer.pal(3,"Set1") )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 1)')

#挑选top10画图展示
Track_genes_sig <- c("SPP1")
cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GC","PBMC")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"GC (Branch 1)"
plot_genes1<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GCPM","PBMC")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"GCPM (Branch 1)"
plot_genes2<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


cds_p1+plot_genes1+plot_genes2







##B2
cds<-readRDS("5.monocle2/9.TAM_B2//TAM_B2_monocle3.rds")


colSet<-c(RColorBrewer::brewer.pal(3,"PuBu")[1:2],
          RColorBrewer::brewer.pal(9,"Blues")[2:7])
cds_p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TAM (Branch 2)')


cluster_p2<-plot_cells(cds, color_cells_by = "Type", label_cell_groups = FALSE, 
                       label_leaves = FALSE,  label_branch_points = FALSE)+
  scale_color_manual(values =c(RColorBrewer::brewer.pal(3,"Set1")[1:2],RColorBrewer::brewer.pal(3,"Set2")[1]) )+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell components (Branch 2)')


#挑选top10画图展示
Track_genes_sig <- c("SPP1")
cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GC")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"GC (Branch 2)"
plot_genes3<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F)+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))


cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GCPM","PE")
                    ) %>%
                    row.names
]

#基因表达趋势图
rownames(cds_subset)<-"GCPM (Branch 2)"
plot_genes4<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,label_by_short_name = F)+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5))





total_p<-ggarrange(Feature_P1,genes_for_branchus1,genes_for_branchus3,genes_for_branchus2,Feature_P2,ncol=5,widths = c(1,2,1,1,1))
pdf("Rplot.pdf",width = 12,height = 4)
print(total_p)
dev.off()