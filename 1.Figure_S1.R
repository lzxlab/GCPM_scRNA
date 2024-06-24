### 以下是堆叠小提琴图展示经典marker基因 #############################################
library(tidyverse)
library(RColorBrewer)
library(scales)
library(reshape2)
library(Seurat)
setwd("~/data/single_cell/3.GCPM1/")



sub_sce<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE","PBMC" ))
dir.create("3.Cluster/13.plot/1.total/")
sub_sce$Cluster<-as.character(sub_sce$Cluster)
sub_sce$Cluster[which(sub_sce$Cluster %in% c("iCAF" ,  "mCAF"))]<-"CAF"
sub_sce$Cluster[which(sub_sce$Cluster %in% c("CD4_conv",  "CD4_Treg"))]<-"CD4"

sub_sce$Cluster<-factor(sub_sce$Cluster,levels=c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,     "MAIT"     , "Tgd"     ,
                                                     "Mast"   ,  "DC"  ,              
                                                     "Mono"       ,            
                                                     "Mph"   ,  "Neutro" ,          
                                                     "CAF", "PC" , "Endo"    ))       


Cluster<-c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,     "MAIT"     , "Tgd"     ,
               "Mast"   ,  "DC"  ,              
               "Mono"       ,            
               "Mph"   ,  "Neutro" ,          
               "CAF", "PC" , "Endo"    )
g.colSet <- readRDS("g.colSet.RDS")
g.colSet1<-c("#1F77B4" ,"#FF7F0E" , "#9467BD", "#AEC7E8","#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5" ,"#C49C94",
             "#F7B6D2" ,"#C7C7C7" ,"#DBDB8D", "#9EDAE5")



library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(Cluster) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

umap_p1<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Type))+
  geom_point(aes(color=Type),size=0.1,alpha=0.6) +
  scale_color_manual(breaks = c(levels(sub_sce$Type)),
                     values = g.colSet1)+
  ggtitle("Tissue type")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
umap_p1


umap_p2<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Sample_ID))+
  geom_point(aes(color=Sample_ID),size=0.1,alpha=0.6) +
  ggtitle("Samples")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
umap_p2

sub_sce$Study<-as.character(sub_sce$Study)
sub_sce$Study[which(grepl("SYSUCC",sub_sce$Study))]<-"SYSUCC"
umap_p3<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Study))+
  geom_point(aes(color=Study),size=0.1,alpha=0.6) +
  ggtitle("Datasets")+
  scale_color_manual(values = RColorBrewer::brewer.pal(2,"Set3"))+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = c(0.2,0.2),
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
umap_p3


Idents(sub_sce)<-sub_sce$Efficiency_ICB
sub_sce1<-subset(sub_sce,idents=c("PR","PD","SD"))
umap_p4<-ggplot(sub_sce1@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Efficiency_ICB))+
  geom_point(aes(color=Efficiency_ICB),size=0.1,alpha=0.6) +
  ggtitle("ICB Efficiency")+
  scale_color_manual(values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.15,0.2),
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
umap_p4




plot_gene=c("CD79A", "JCHAIN"  , 
            "CD3D","CD4" ,   "CD8A","GNLY"  ,
            "CD1C" ,"IDO1" ,"CD14","CD68","S100A9", 
            
              
            "CPA3","CSF3R",
            
                
            
            "PECAM1","SPARC","LUM","RGS5","ACTA2"
            
)


feature_plot<-FeaturePlot(sub_sce,features =plot_gene,pt.size = 2,max.cutoff = 3,min.cutoff = 0,
                          combine = F,raster = T) 

feature_list<-list()
k=0
for (i in 1:length(feature_plot)){
  tmp_p<-feature_plot[[i]]
  tmp_p<-tmp_p+theme_classic()+ 
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(5,"BuPu"))+
    theme(axis.title = element_blank(), legend.text = element_text(size = 10),
          axis.text = element_blank(),
          axis.line = element_blank(),axis.ticks = element_blank(),
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          panel.border = element_rect(fill=NA,color="grey"),
          legend.position="none")
  
  
  k=k+1
  if (k %in% c(6,12,18)){
    tmp_p<- tmp_p+theme(legend.position="right")
  }
  
  feature_list[[k]]<-tmp_p
}


feature_ps<-ggarrange(plotlist =feature_list,ncol=6,nrow=3,
                      widths = c(1,1,1,1,1,1.4))







##Stack plot

sub_sce<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE","PBMC" ))
dir.create("3.Cluster/13.plot/1.total/")
sub_sce$Cluster<-as.character(sub_sce$Cluster)
sub_sce$Cluster[which(sub_sce$Cluster %in% c("iCAF" ,  "mCAF"))]<-"CAF"
sub_sce$Cluster[which(sub_sce$Cluster %in% c("CD4_conv",  "CD4_Treg"))]<-"CD4"

sub_sce$Cluster<-factor(sub_sce$Cluster,levels=c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,     "MAIT"     , "Tgd"     ,
                                                     "Mast"   ,  "DC"  ,              
                                                     "Mono"       ,            
                                                     "Mph"   ,  "Neutro" ,          
                                                     "CAF", "PC" , "Endo"    ))   
combined<-sub_sce

##All cell type
tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,Sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,Type,Sample_ID),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,Sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by=c("Type","Sample_ID"))
sum_table$per<-sum_table$num/sum_table$total_num*100


Cluster<-c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,     "MAIT"     , "Tgd"     ,
               "Mast"   ,  "DC"  ,              
               "Mono"       ,            
               "Mph"   ,  "Neutro" ,          
               "CAF", "PC" , "Endo"    )
g.colSet <- readRDS("/home/zhengyq/data/single_cell/7.GC_gut//1.data/colSet.list.rds")
g.colSet1<-c("#1F77B4" ,"#FF7F0E" , "#9467BD", "#AEC7E8","#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5" ,"#C49C94",
             "#F7B6D2" ,"#C7C7C7" ,"#DBDB8D", "#9EDAE5")
names(g.colSet1)<-Cluster

sum_table1<-sum_table[sum_table$Cluster %in% c("B","Plasma"),]
sum_table1<-dplyr::summarize(group_by(sum_table1,Type,Sample_ID),per=sum(per))

sum_table1<-sum_table1[order(sum_table1$per,decreasing = T),]
Sample_ID<-unique(c(sum_table1$Sample_ID,sum_table$Sample_ID))
sum_table$Sample_ID<-factor(sum_table$Sample_ID,levels=Sample_ID)

sum_table$Type<-factor(sum_table$Type,levels=c("GN" , "GC","GCPM","PE","PBMC"))
stack_p2<-ggplot(data=sum_table, aes(x=Sample_ID, y=per, fill=Cluster)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.8)+ 
  theme_classic()+
  scale_fill_manual(values = g.colSet1)+
  labs(x = 'Type', y ="Percentage",title=paste("")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,101),breaks = seq(0,100,25),labels = c("0","25%","50%","75%","100%"))+
  theme(axis.text.y = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "right",
        axis.text.x = element_blank(),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10,vjust = 0.5,hjust=0.5)
  )+
  facet_grid(~Type,scales = "free_x", space = 'free_x')
print(stack_p2)





sub_sce<-readRDS("3.Cluster/13.Annotation/6.Epi_sub_sce_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM" ))
dir.create("3.Cluster/13.plot/1.total/")


sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c( "Epi_C1-LGALS3" ,  "Epi_C2-FTL"   , "Epi_C3-NEAT1" ,
                                                        "Epi_C4-GKN2"  , "Epi_C5-prolif" ,"Epi_C6-KRT20" , "Epi_C7-SPARC" , "Epi_C8-PHGR1"  ,"Epi_C9-CKB",
                                                        "Epi_C10-PGC"  , "Epi_C11-SCG5" , "Epi_C12-PGA5"  ))       


Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "Epi_C1"        ,       
                   "Epi_C2"   ,            
                   "Epi_C3"        ,       "Epi_C4"   ,           
                   "Epi_C5"       ,       "Epi_C6" ,          
                   "Epi_C7"      ,      "Epi_C8" ,  "Epi_C9",
                   "Epi_C10"  , "Epi_C11" , "Epi_C12") 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


Cluster<-unique(sort(sub_sce$SubCluster))
g.colSet <- readRDS("g.colSet.RDS")
g.colSet1<-g.colSet$Set1[1:length(Cluster)]
names(g.colSet1)<-Cluster
#g.colSet1<-list("Cluster"=g.colSet1)



#g.colSet1<-list("SubCluster"=g.colSet1)

library(ggforce)
library(ggpubr)
library(ggrepel)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(SubCluster1) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )


epi_p1<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c( "Epi_C1-LGALS3" ,  "Epi_C2-FTL"   , "Epi_C3-NEAT1" ,
                                "Epi_C4-GKN2"  , "Epi_C5-prolif" ,"Epi_C6-KRT20" , "Epi_C7-SPARC" , "Epi_C8-PHGR1"  ,"Epi_C9-CKB",
                                "Epi_C10-PGC"  , "Epi_C11-SCG5" , "Epi_C12-PGA5"    ),
                     values = g.colSet1)+
  ggtitle("Clustering of epithelial cells")+
  geom_text_repel(aes(label = SubCluster1), data = class_avg,color="black",size=2.5)+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
epi_p1



epi_p2<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Type))+
  geom_point(aes(color=Type),size=0.1) +
  scale_color_manual(values = c("#1F77B4" ,"#FF7F0E" , "#9467BD"))+
  ggtitle("Tissue type")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
epi_p2


epi_p3<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Cluster))+
  geom_point(aes(color=Cluster),size=0.1) +
  scale_color_manual(values = c("#1F77B4" ,"#FF7F0E" , "#9467BD"))+
  ggtitle("CNV")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=1.5)))
epi_p3




total_p1<-umap_p1+umap_p3+umap_p4+umap_p2+plot_layout(ncol=4)

total_p2<-ggarrange(feature_ps,NULL,
                    ncol=2,widths = c(10,0.5))

total_p3<-stack_p2

right_plot<-epi_p1+epi_p2/epi_p3+plot_layout(ncol=2,widths = c(2,1))

total_p4<-ggarrange(right_plot,NULL,ncol=2,widths=c(1.5,1))

total_p<-ggarrange(total_p1,total_p2,stack_p2,total_p4,ncol=1,nrow = 4,heights = c(3,5,4,3))
pdf("14.Figure_plot/Figure_1_S1.pdf",width = 12,height = 15)
print(total_p)
dev.off()


png("14.Figure_plot/Figure_1_S1.png",width = 12,height = 15,units = "in",res = 300)
print(total_p)
dev.off()
