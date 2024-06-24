library(Seurat)
library(data.table)
library("sscVis")
library("sscClust")
library("scPip")
library(ggplot2)
library(ggridges)
library(dplyr)
library(cowplot)
library(PCAtools)
library(ggpubr)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/3.GCPM1/")



# 展示特定行名函数
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  library(grid)
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  return(heatmap)
}

test<-readRDS("5.monocle2/3.Mono/monocle2.RDS")
test$Type<-factor(test$Type,levels = c("PE","GCPM","GN","GC"))

library(monocle)
monocle2_time<-plot_cell_trajectory(test,color_by = "Pseudotime",
                                    cell_size =0.4)+
  theme(legend.position = "right")+
  scale_color_gradientn(colors = c(RColorBrewer::brewer.pal(5,"Purples")))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory of TRM-derived Mph')+
  facet_wrap(~Type)
monocle2_time


monocle2_state<-plot_cell_trajectory(test,color_by = "State",
                                    cell_size =0.4)+
  theme(legend.position = "right")+
  #scale_color_gradientn(colors = c(RColorBrewer::brewer.pal(5,"Purples")))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell state')+
  facet_wrap(~Type)+
  guides(colour = guide_legend(ncol=1,override.aes = list(size=3)))
monocle2_state



pseudotime_tab<-data.frame(Cell_ID=colnames(test),
                           Pseudotime=test$Pseudotime,
                           Type=test$Type,
                           Cluster=test$Cluster,
                           State=test$State)

monocle2_density<-ggplot(pseudotime_tab,aes(x=Pseudotime,color=Cluster))+
  geom_density(size=0.8)+
  facet_grid(Cluster~.,switch = "y")+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  ggtitle('Cell density by pseudotime')+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position="none",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))

monocle2_density2<-ggplot(pseudotime_tab,aes(x=Pseudotime,color=Cluster))+
  geom_density(size=0.8)+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  facet_grid(Cluster~Type,switch = "y",scales = "free_y")+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))
monocle2_density2



##TF activity
cg_n<-readRDS("7.scenic/TAM/cg_n.RDS")
seu<-readRDS("7.scenic/TAM/seu_obj.RDS")
Idents(seu)<-seu$Type
seu<-subset(seu,idents=c("GCPM","GC","PE","GN"))
monocle1<-readRDS("5.monocle2/3.Mono/monocle2.RDS")

pseudotime_tab<-data.frame(Cell_ID=colnames(monocle1),
                           Pseudotime=monocle1$Pseudotime,
                           Cluster=monocle1$Cluster)
rownames(pseudotime_tab)<-pseudotime_tab$Cell_ID
cg_n<-AddMetaData(cg_n,metadata = seu@meta.data)
cg_n<-AddMetaData(cg_n,metadata = pseudotime_tab)
seu<-AddMetaData(seu,metadata = pseudotime_tab)

cg_n$Cluster<-factor(cg_n$Cluster,levels=c("F13A1+ TRM","CXCL9+ Mph","TAM"))
seu$Cluster<-factor(seu$Cluster,levels=c("F13A1+ TRM","CXCL9+ Mph","TAM"))
Idents(cg_n)<-cg_n$Cluster


gene_list<-c("KDM5A" , "ELF1","HIC1","HDAC2","MLX" ,"CEBPA")

seu$tf_FOXO3<-as.numeric(cg_n@assays$RNA@counts["FOXO3",])
seu$tf_MAF<-as.numeric(cg_n@assays$RNA@counts["MAF",])
seu$tf_NFATC1<-as.numeric(cg_n@assays$RNA@counts["NFATC1",])
seu$tf_IRF8<-as.numeric(cg_n@assays$RNA@counts["IRF8",])
seu$tf_HIC1<-as.numeric(cg_n@assays$RNA@counts["HIC1",])
seu$tf_CEBPA<-as.numeric(cg_n@assays$RNA@counts["CEBPA",])
seu$tf_MLX<-as.numeric(cg_n@assays$RNA@counts["MLX",])
seu$tf_TCF4<-as.numeric(cg_n@assays$RNA@counts["TCF4",])
seu$tf_ETV5<-as.numeric(cg_n@assays$RNA@counts["ETV5",])
DefaultAssay(seu)<-"RNA"
seu<-ScaleData(seu)
FeatureP1<-FeaturePlot(seu,features = c("tf_MAF" , "tf_FOXO3","tf_TCF4" ,"tf_IRF8","tf_NFATC1","tf_HIC1" ,"tf_MLX","tf_ETV5","tf_CEBPA"),
                       ncol=3,
                       min.cutoff = 0,
                       max.cutoff = 2,combine = F) 

Vln_list<-list()
for (i in 1:length(FeatureP1)){
  tmp_p<-FeatureP1[[i]]
  tmp_p<-tmp_p&
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(5,"BuPu")) &
    theme(axis.text = element_blank(), axis.title = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_blank(), 
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          legend.title  = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          panel.border = element_rect(fill = NA,color="grey"))
  if(i==3){
    tmp_p<-tmp_p+theme(legend.position = "right")
  }
  Vln_list[[i]]<-tmp_p
}

library(patchwork)
FeatureP1_1<-Vln_list[[1]]+Vln_list[[2]]+Vln_list[[3]]+
  Vln_list[[4]]+Vln_list[[5]]+Vln_list[[6]]+
  Vln_list[[7]]+Vln_list[[8]]+Vln_list[[9]]+plot_layout(ncol=3,nrow=3)


FeatureP1


FeatureP2<-VlnPlot(seu,features = c("tf_MAF" , "tf_FOXO3","tf_TCF4" ,"tf_IRF8","tf_NFATC1","tf_HIC1" ,"tf_MLX","tf_ETV5","tf_CEBPA"),
                       ncol=3,pt.size = 0,group.by = "Cluster",combine=F)

Vln_list<-list()
for (i in 1:length(FeatureP2)){
  tmp_p<-FeatureP2[[i]]
  tmp_p<-tmp_p+
    stat_compare_means(label = "p.format",label.y.npc=0.9,label.x.npc=0.1) +
    scale_fill_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                       labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                       values = RColorBrewer::brewer.pal(3,"Set2"))+
    theme(axis.text = element_blank(), axis.title = element_blank(), 
          axis.ticks = element_blank(), axis.line = element_blank(), 
          legend.text = element_text(size = 10),
          axis.text.x = element_blank(), 
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          legend.title  = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          panel.border = element_rect(fill = NA,color="grey"))
  if(i>=7){
    tmp_p<-tmp_p+theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1))
  }
  Vln_list[[i]]<-tmp_p
}
  
FeatureP2_1<-Vln_list[[1]]+Vln_list[[2]]+Vln_list[[3]]+
  Vln_list[[4]]+Vln_list[[5]]+Vln_list[[6]]+
  Vln_list[[7]]+Vln_list[[8]]+Vln_list[[9]]+plot_layout(ncol=3,nrow=3)



DefaultAssay(seu)<-"RNA"
seu<-ScaleData(seu)

gene_list<-c("SPP1","TREM2","CTSD","GPNMB")
gene_tab<-data.frame(Cell_ID=colnames(seu),
                     t(seu@assays$RNA@scale.data[gene_list,]))
gene_tab<-reshape2::melt(gene_tab)
names(gene_tab)<-c("Cell_ID","gene_name","value")
gene_info<-data.frame(Cell_ID=colnames(seu),
                      tf_CEBPA=seu$tf_CEBPA)
merged_tab<-merge(gene_tab,gene_info,by="Cell_ID")
tf_CEBPA<-ggplot(data=merged_tab,aes(tf_CEBPA,value))+
  geom_point(color="grey")+
  ggtitle("Correlation of CEBPA activity with TAM markers")+
  xlab("Activity of CEBPA")+
  ylab("Gene expression")+
  geom_smooth()+
  stat_cor()+
  facet_wrap(~gene_name,scales = "free_y",ncol=2)+
  theme(axis.text = element_text(size = 10), axis.title =element_text(size = 10), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10), 
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        legend.title  = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))
tf_CEBPA




library(SCENIC)
packageVersion("SCENIC")  
library(SCopeLoomR)
scenicLoomPath='7.scenic/TAM/auc_mtx.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)






##T_NK
TF_name="MAF"
Cluster_name="Mph_C2_TRM-F13A1"
combined<-readRDS("7.scenic/TAM/seu_obj.RDS")
combined$Cluster<-factor(combined$Cluster,levels = c("F13A1+ TRM","CXCL9+ Mph","TAM"))
markers<-regulons[[which(names(regulons)==paste(TF_name,"(+)",sep=""))]]

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined<-subset(combined, features =markers )

Exp_table<-data.frame(Cluster=combined$Cluster,t(combined@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))


pheatmapp<-pheatmap::pheatmap(ave_Exp_table,
                              color=colorRampPalette(c("navy" ,'white','#7E2324'), bias=1)(50), border_color=NA,
                              
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              treeheight_row = 10,
                              main = paste("Regulons of",TF_name),
                              angle_col = 90
                              
)

cell_markers<-read.table("3.Cluster/8.SubAnnotation/2.myeloid/sub_sce_Cluster_markers.txt")
top10 <-  top_n(group_by(cell_markers,cluster), n = 30, wt = avg_log2FC)
top10_gene<-top10$gene[which(top10$cluster==Cluster_name)]

gene_name<-intersect(rownames(ave_Exp_table),top10_gene)
pheatmapp1<-add.flag(pheatmapp,kept.labels = gene_name,repel.degree = 0.2)
pheatmapp1<-plot_grid(pheatmapp1)



##T_NK
TF_name="MLX"
Cluster_name="Mph_C4_TAM-SPP1"
combined<-readRDS("7.scenic/TAM/seu_obj.RDS")
combined$Cluster<-factor(combined$Cluster,levels = c("F13A1+ TRM","CXCL9+ Mph","TAM"))
markers<-regulons[[which(names(regulons)==paste(TF_name,"(+)",sep=""))]]

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined<-subset(combined, features =markers )

Exp_table<-data.frame(Cluster=combined$Cluster,t(combined@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))


pheatmapp<-pheatmap::pheatmap(ave_Exp_table,
                              color=colorRampPalette(c("navy" ,'white','#7E2324'), bias=1)(50), border_color=NA,
                              
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              treeheight_row = 10,
                              main = paste("Regulons of",TF_name),
                              angle_col = 90
                              
)

cell_markers<-read.table("3.Cluster/8.SubAnnotation/2.myeloid/sub_sce_Cluster_markers.txt")
top10 <-  top_n(group_by(cell_markers,cluster), n = 30, wt = avg_log2FC)
top10_gene<-top10$gene[which(top10$cluster==Cluster_name)]

gene_name<-intersect(rownames(ave_Exp_table),top10_gene)
pheatmapp2<-add.flag(pheatmapp,kept.labels = gene_name,repel.degree = 0.2)
pheatmapp2<-plot_grid(pheatmapp2)






##T_NK
TF_name="ETV5"
Cluster_name="Mph_C4_TAM-SPP1"
combined<-readRDS("7.scenic/TAM/seu_obj.RDS")
combined$Cluster<-factor(combined$Cluster,levels = c("F13A1+ TRM","CXCL9+ Mph","TAM"))
markers<-regulons[[which(names(regulons)==paste(TF_name,"(+)",sep=""))]]

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined<-subset(combined, features =markers )

Exp_table<-data.frame(Cluster=combined$Cluster,t(combined@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))


pheatmapp<-pheatmap::pheatmap(ave_Exp_table,
                              color=colorRampPalette(c("navy" ,'white','#7E2324'), bias=1)(50), border_color=NA,
                              
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              treeheight_row = 10,
                              main = paste("Regulons of",TF_name),
                              angle_col = 90
                              
)

cell_markers<-read.table("3.Cluster/8.SubAnnotation/2.myeloid/sub_sce_Cluster_markers.txt")
top10 <-  top_n(group_by(cell_markers,cluster), n = 40, wt = avg_log2FC)
top10_gene<-top10$gene[which(top10$cluster==Cluster_name)]
top10_gene<-top10_gene[which(top10_gene!="C3")]

gene_name<-intersect(rownames(ave_Exp_table),top10_gene)
pheatmapp3<-add.flag(pheatmapp,kept.labels = gene_name,repel.degree = 0.2)
pheatmapp3<-plot_grid(pheatmapp3)




##T_NK
TF_name="IRF8"
Cluster_name="Mph_C3-CXCL9"
combined<-readRDS("7.scenic/TAM/seu_obj.RDS")
combined$Cluster<-factor(combined$Cluster,levels = c("F13A1+ TRM","CXCL9+ Mph","TAM"))
markers<-regulons[[which(names(regulons)==paste(TF_name,"(+)",sep=""))]]

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined<-subset(combined, features =markers )

Exp_table<-data.frame(Cluster=combined$Cluster,t(combined@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))


pheatmapp<-pheatmap::pheatmap(ave_Exp_table,
                              color=colorRampPalette(c("navy" ,'white','#7E2324'), bias=1)(50), border_color=NA,
                              
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              treeheight_row = 10,
                              main = paste("Regulons of",TF_name),
                              angle_col = 90
                              
)

cell_markers<-read.table("3.Cluster/8.SubAnnotation/2.myeloid/sub_sce_Cluster_markers.txt")
top10 <-  top_n(group_by(cell_markers,cluster), n = 30, wt = avg_log2FC)
top10_gene<-top10$gene[which(top10$cluster==Cluster_name)]

gene_name<-intersect(rownames(ave_Exp_table),top10_gene)
pheatmapp4<-add.flag(pheatmapp,kept.labels = gene_name,repel.degree = 0.2)
pheatmapp4<-plot_grid(pheatmapp4)



total_p1<-monocle2_density+monocle2_time+monocle2_state+plot_layout(ncol=3,widths = c(0.8,1,1))
total_p3<-ggarrange(FeatureP1_1,FeatureP2_1,ncol=2,widths = c(1,0.8))
total_p4<-ggarrange(pheatmapp1,pheatmapp2,pheatmapp3,pheatmapp4,ncol=4,common.legend = T,legend="right")
  
total_p<-ggarrange(total_p1,monocle2_density2,total_p3,total_p4,nrow = 4,ncol = 1,heights=c(3,3,5,4))
pdf("14.Figure_plot/Figure3_S6.pdf",width = 12,height = 15)
print(total_p)
dev.off()