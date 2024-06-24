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
library(RColorBrewer)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/3.GCPM1/")

# obj: seurat obj
# group_mat: matrix with 2 cols, Cell_ID and Cluster
# markers: gene_list to calculate
Cal_average_exp<-function(obj,group_mat,markers){
  DefaultAssay(obj)<-"RNA"
  obj<-ScaleData(obj)
  obj<-subset(obj,features = markers)
  Exp<-t(obj@assays$RNA@scale.data)
}


##T_NK
combined1<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")

Idents(combined1)<-combined1$Cluster
combined<-subset(combined1,idents=c("Mono","Mph","DC","Neutro" ))

markers<-c(
  "LILRA4","LAMP3","IDO1","CD1C","AREG",
  "MNDA","FCGR3B","CSF3R","CXCR2","NAMPT",
  "IL1B","S100A8","S100A9","CD14",
  
  "FCGR3A","CD68","CD163","C1QA","C1QC",
  
  
  
  "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
  
  "CD274","PDCD1LG2","CD276","VSIG4","LGALS9","HAVCR2",
  
  "CXCL1","CXCL8","CXCL9","CXCL10","SPP1","MMP12"
)  

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined1<-subset(combined, features =markers )
combined1$SubCluster<-factor(combined1$SubCluster,levels=c(    "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                                                               "DC_C4_cDC2-CD1C"   ,  
                                                               
                                                               
                                                               "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                                               "Mph_C1-TIMP1"    ,
                                                               "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                               "Mph_C5_TAM-FBP1",
                                                               "Neutro_C1-CD16B" 
                                                               
))   

Exp_table<-data.frame(Cluster=combined1$SubCluster,t(combined1@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))

colData1<-data.frame(Cluster=c(rep("DC",4),rep("Monocytes",3),rep("TAM",5),rep("Neutrophil",1) ))
rownames(colData1)<-colnames(ave_Exp_table)

rowData1<-data.frame(Markers=c(
  rep("DC markers",5),
  rep("Neutro markers",5),
  rep("Mono markers",4),
  rep("Macro markers",5),
  rep("MHC II molecular",4),
  rep("Inhibitory receptors",6),
  rep("Cytokines",6)
))
rownames(rowData1)<-markers
ave_Exp_table<-ave_Exp_table[markers,]

pheatmap_p1<-pheatmap::pheatmap(ave_Exp_table,
                                color=colorRampPalette(c("navy" ,"navy" ,'white','#D97777','#7E2324'), bias=1)(50), border_color=NA,
                                #color=colorRampPalette(c('#3C5488FF','white','#DC0000FF'), bias=1)(50), border_color=NA,
                                
                                #color=colorRampPalette(c('red',"blue"), bias=1)(50), border_color=NA,
                                annotation_row = rowData1,
                                annotation_col = colData1,
                                cluster_rows = F,
                                cluster_cols = F,
                                scale = "row",
                                gaps_col=c(4,7,12),
                                gaps_row = c(5,10,14,19,23,29)
)
print(pheatmap_p1)
pheatmap_p1<-plot_grid(pheatmap_p1$gtable)






library(GSVA)
combined1<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")

mymatrix <- combined1@assays$RNA@counts
mymatrix <- as.matrix(mymatrix)
mymatrix[1:5,1:5]
##load the gene set 

features=list("M1"=c("TNF","IRF5","CD86","IL1B","KYNU","IRF1","CCR7","CD40","CXCL9","IL23A"),
              "M2"=c("CLEC7A","CCL4","TGFB1","CTSD","CTSC","CTSB","CTSA","FN1","MSR1","TNFSF12","LYVE1","IL4R","MMP19","CCL20","MMP14","CCL18","CD276","SPP1","VEGFB"),
              "Angiogenesis"=c("TYMP","VCAN","CD44","FYN","VEGFA","TNFAIP6","E2F3","MMP9","ITGAV","SPP1","CXCR4","PTK2","CCND2","EZH2"),
              "Phagocytosis"=c("CD163","C1QB","MERTK","MRC1"),
              "Checkpoint"=c("LAIR1","SIRPA","HAVCR2","NECTIN2","PVR","CD274","ADORA2A","IDO1","BTLA")
)

es.dif <- gsva(mymatrix, features, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=20)

melt_tab<-reshape2::melt(es.dif)
colnames(melt_tab)<-c("Pathway","Cell_ID","value")
cell_info<-data.frame(Cell_ID=colnames(combined1),
                      Cluster=combined1$SubCluster,
                      Type=combined1$Type)

merged_tab<-merge(melt_tab,cell_info,by="Cell_ID")
merged_tab<-merged_tab[which(merged_tab$Pathway %in% c("M1","M2","Checkpoint")),]
FSbox_p<-ggplot(merged_tab, aes(Cluster, value, fill = Cluster))+
  geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
  #geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format",label.y.npc = 0.9)+
  theme_classic()+
  #scale_fill_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), 
        legend.position = "none",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,vjust=1,hjust=1),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5,face="bold"),
        strip.background = element_blank())+
  facet_wrap(~Pathway,scales = "free_y")

print(FSbox_p)


FSbox_p1<-ggplot(merged_tab, aes(Cluster, value, fill = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0.8),outlier.shape=NA)+
  #geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format",label.y.npc = 0.9)+
  theme_classic()+
  #scale_fill_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), 
        legend.position = "right",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,vjust=1,hjust=1),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5,face="bold"),
        strip.background = element_blank())+
  facet_wrap(~Pathway,scales = "free_y")

print(FSbox_p1)
VlnPlot(combined1,features = c("M1","M2","Angiogenesis","Phagocytosis","Checkpoint"),group.by = "SubCluster",pt.size = 0)


##T_NK
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$SubCluster
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)

Vln_P1<-VlnPlot(sub_sce,features = c("CXCL9","MMP12","SPP1","TREM2","CD274","HAVCR2") ,ncol = 2,pt.size=0,combine = F)


ps_list<-list()
k=0
for (i in Vln_P1){
  i<-i+theme(plot.title = element_text(size = 10,hjust = 0.5),
             axis.title.x = element_text(size = 10),
             axis.title.y = element_text(size = 10),
             axis.text.x = element_text(size = 9),
             axis.text.y = element_text(size = 10),
             legend.title = element_text(size = 10),
             
             legend.text = element_text(size = 10),
             legend.position = "none")+
    scale_y_continuous(expand = c(0,0),limits = c(-1,ceiling(max(i$data[[1]]))))
  
  
  k=k+1
  
  if (k<=3){
    i<-i+theme(axis.title.x = element_blank(),
               axis.text.x = element_blank())
  }
  
  ps_list[[k]]<-i
}

Vln_P1<-ps_list[[1]]+ps_list[[2]]+ps_list[[3]]+ps_list[[4]]+ps_list[[5]]+ps_list[[6]]+
  plot_layout(ncol=3,nrow = 2)
#Vln_P1<-ggarrange(plotlist = ps_list,heights = c(1,1,1.5),ncol = 5,nrow = 2)




##T_NK
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("GC","GCPM"))
Idents(sub_sce)<-sub_sce$Cluster
sub_sce<-subset(sub_sce,idents=c("Mono","Mph"))

Idents(sub_sce)<-sub_sce$SubCluster
Vln_P2<-VlnPlot(sub_sce,features = c("SPP1") ,
                ncol = 1,pt.size=0,
                split.by = "Type",cols =RColorBrewer::brewer.pal(3,"Set3")[1:2] )+
  stat_compare_means(label = "p.signif",label.y = 9.5)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10),
        legend.position = "right")+
  scale_y_continuous(expand = c(0,0),limits = c(-0.5,11))


##Feature_plot
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("GC","GCPM"))
Feature_P1<-FeaturePlot(sub_sce,features = c("SPP1") ,
                ncol = 1,pt.size=0,
                split.by = "Type",cols = viridis::viridis(3,option = "D")) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10),
        legend.position = "none")





##All Cells
combined<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(combined)<-combined$Type
combined<-subset(combined,idents=c("GC","GCPM"))

Idents(combined)<-combined$SubCluster
combined<-subset(combined,idents=c("Mph_C4_TAM-SPP1"  ))
combined$SubCluster<-as.character(combined$SubCluster)
combined$SubCluster[which(combined$SubCluster=="Mph_C4_TAM-SPP1")]<-"SPP1+ TAM"

tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,Type),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="Type")
sum_table$per<-sum_table$num/sum_table$total_num*100


Cluster<-unique(sum_table$Cluster)
g.colSet <- readRDS("/home/zhengyq/data/single_cell/7.GC_gut//1.data/colSet.list.rds")
g.colSet1<-g.colSet$Set8[1:length(Cluster)]
names(g.colSet1)<-Cluster

sum_table$Type<-factor(sum_table$Type,levels=c("GN" , "GC","GCPM","PE"))
stack_p<-ggplot(data=sum_table, aes(x=Type, y=per, fill=Cluster)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.8)+ 
  theme_classic()+
  scale_fill_manual(values = g.colSet1)+
  labs(x = 'Type', y ="Component",title=paste("Cell proportion")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,101),breaks = seq(0,100,25),labels = c("0","25%","50%","75%","100%"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(size = 10,vjust = 1,hjust=1,angle=45),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
  )
print(stack_p)





library(monocle3)

cds<-readRDS("5.monocle2/7.TAM_monocle3/TAM_monocle3.rds")

colSet<-c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
          brewer.pal(8, "Set1"),
          brewer.pal(8, "Dark2")[1],
          "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(colSet)<-c(   "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                    "DC_C4_cDC2-CD1C"   ,  
                    
                    
                    "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                    "Mph_C1-TIMP1"    ,
                    "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                    "Mph_C5_TAM-FBP1",
                    "Neutro_C1-CD16B"   )

cds_p1<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme_classic()+
  scale_color_gradientn(colors = c(RColorBrewer::brewer.pal(5,"Purples")))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        legend.title  = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory of Mono/TAM')
cds_p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme_classic()+
  xlab("UMAP_1")+ylab("UMAP_2")+
  facet_wrap(~Type,ncol=1,nrow = 5)+
  ggtitle('Pseudotime of Mono/TAM')+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        legend.title  = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  theme(legend.position = "right")

cds_p3<-plot_cells(cds, color_cells_by = "SubCluster", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme_classic()+
  ggtitle('Cell trajectory of Mono/TAM')+
  scale_color_manual(breaks = c("Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                "Mph_C1-TIMP1"    ,
                                "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                "Mph_C5_TAM-FBP1"),
                     labels =c("Mono_C1"  ,
                               "Mono_C2"   ,  "Mono_C3"  , 
                               "Mph_C1"    ,
                               "Mph_C2_TRM" ,"Mph_C3", "Mph_C4_TAM", 
                               "Mph_C5_TAM"),
                     values = c("#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5",
                                "#b3b3b3" ,"#E41A1C" ,"#377EB8", "#4DAF4A", "#984EA3"))+
  facet_wrap(~Type,ncol=1,nrow = 5)+
  xlab("UMAP_1")+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        legend.title  = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  theme(legend.position = "right")


combined <- readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mono","Mph"))

Idents(combined)<-combined$Type
combined<-subset(combined,idents=c("GCPM","PE","PBMC","GC","GN"))

#combined$Cluster<-combined$SubCluster
DefaultAssay(combined)<-"RNA"

combined<-AddMetaData(combined,combined@reductions$umap@cell.embeddings,col.name = colnames(combined@reductions$umap@cell.embeddings))
combined$Branchus<-"Other"
combined$Branchus[which((combined$UMAP_1<3&combined$UMAP_2<2&combined$UMAP_1>(-4.5)&combined$SubCluster=="Mph_C1-TIMP1")|combined$SubCluster %in% c("Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A" ))]<-"Mono-derived"
combined$Branchus[which(combined$UMAP_1<(-0.5)& combined$SubCluster %in% c("Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                                           "Mph_C5_TAM-FBP1"))]<-"TRM-derived"

Idents(combined)<-combined$Branchus
combined<-subset(combined,idents=c("Mono-derived","TRM-derived"))

combined$Pseudotime<-pseudotime(cds)

gene_list<-c("SPP1" , "TREM2","C3AR1","MSR1","FCN1","IL1B","CD44","VCAN")

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

select_TF_mt<-as.matrix(combined@assays$RNA@scale.data[gene_list,])
select_TF_mt<-reshape2::melt(select_TF_mt)
colnames(select_TF_mt)<-c("TF_names","Cell_ID","value")

merged_tab<-merge(combined@meta.data,select_TF_mt,by="Cell_ID")
EXP_time<-ggplot(merged_tab,aes(x=Pseudotime,y=value,color=Branchus))+
  geom_smooth(se = F)+
  ggtitle("Expression by pseudotime")+
  theme_classic2()+
  scale_color_manual(values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "top",
        panel.border = element_rect(fill = NA,color="grey"))+
  facet_wrap(~TF_names,scales="free_y",ncol=2)
EXP_time

melt_tab<-reshape2::melt(es.dif)
colnames(melt_tab)<-c("Pathway","Cell_ID","value")

merged_tab<-merge(merged_tab,melt_tab,by="Cell_ID")
merged_tab<-merged_tab[merged_tab$Pathway %in% c("M1","M2","Angiogenesis","Checkpoint"),]
Pathway_time<-ggplot(merged_tab,aes(x=Pseudotime,y=value.y,color=Branchus))+
  geom_smooth(se = F)+
  ggtitle("Pathway by pseudotime")+
  theme_classic2()+
  scale_color_manual(values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA,color="grey"))+
  facet_wrap(~Pathway,scales="free_y",ncol=1)
Pathway_time


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
Idents(cg_n)<-cg_n$Cluster


gene_list<-c("KDM5A" , "ELF1","HIC1","HDAC2","MLX" ,"CEBPA")

seu$tf_KDM5A<-as.numeric(cg_n@assays$RNA@counts["KDM5A",])
seu$tf_ELF1<-as.numeric(cg_n@assays$RNA@counts["ELF1",])
seu$tf_HIC1<-as.numeric(cg_n@assays$RNA@counts["HIC1",])
seu$tf_HDAC2<-as.numeric(cg_n@assays$RNA@counts["HDAC2",])
seu$tf_MLX<-as.numeric(cg_n@assays$RNA@counts["MLX",])
seu$tf_CEBPA<-as.numeric(cg_n@assays$RNA@counts["CEBPA",])
DefaultAssay(seu)<-"RNA"
seu<-ScaleData(seu)
FeatureP1<-FeaturePlot(seu,features = c("tf_KDM5A" , "tf_ELF1","tf_HIC1","tf_HDAC2","tf_MLX" ,"tf_CEBPA"),
            ncol=6,
            min.cutoff = 0,
            max.cutoff = 2) &
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(5,"BuPu")) &
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        legend.title  = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))&
  theme(legend.position = "none") +
  theme(legend.position = "right") 


FeatureP2<-FeaturePlot(seu,features = c("KDM5A" , "ELF1","HIC1","HDAC2","MLX" ,"CEBPA"),
                      ncol=6,
                      min.cutoff = 0) &
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(5,"BuPu")) &
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        legend.title  = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))&
  theme(legend.position = "none") +
  theme(legend.position = "right") 


tab<-data.frame(SPP1=as.numeric(seu@assays$RNA@scale.data["CEBPA",]),
                seu@meta.data)
ggplot(data=tab,aes(tf_CEBPA,SPP1))+
  geom_point()+
  geom_smooth()+
  stat_cor()


total_p1<-ggarrange(pheatmap_p1,FSbox_p+Vln_P1+plot_layout(ncol=1,nrow=2,heights = c(1,1.5)),ncol=2,widths = c(1,1.4))
total_p2<-ggarrange(Vln_P2,Feature_P1,stack_p,ncol=3,nrow = 1,widths  = c(1.2,1.4,0.7))
total_p3<-cds_p1+cds_p2+plot_layout(ncol=2,nrow = 1)
total_p4<-cluster_p1+cluster_p2+plot_layout(ncol=1,nrow = 2)
total_p5<-ggarrange(cds_p2+cds_p3,EXP_time+Pathway_time+plot_layout(width=c(2,1)),ncol=2,nrow =1,widths  = c(1,1))


total_p<-ggarrange(total_p1,total_p5,FeatureP1 ,FeatureP2,ncol=1,nrow = 4,heights = c(7,6,2,2))
pdf("14.Figure_plot/Figure3_S5.pdf",width = 12,height = 17)
print(total_p)
dev.off()
