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




CD4_sce<-readRDS("3.Cluster/13.Annotation/1.CD4_sub_sce_annotation.rds")
CD4_sce<-AddMetaData(CD4_sce,CD4_sce@reductions$umap@cell.embeddings,col.name = colnames(CD4_sce@reductions$umap@cell.embeddings))
CD4_sce$SubCluster<-factor(CD4_sce$SubCluster,levels=c("CD4_C1_Tn-CCR7"   ,  "CD4_C2_Tcm-AXNA1" , "CD4_C3_Tfh-CXCL13",
                                                       "CD4_C4_Tc-GZMA"  ,  "CD4_C5_Treg-FOXP3" ,"CD4_C6-EGR1"  ,  "CD4_C7-NEAT1" ,
                                                       "CD4_C8-TIMP1"  , "CD4_C9-ISG"  ,  "CD4_C10-prolif"   ))

CD4_sce$SubCluster1<-as.numeric(CD4_sce$SubCluster)

CD4_sce$TopCluster<-"CD4+ T cells"

CD8_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")
CD8_sce<-AddMetaData(CD8_sce,CD8_sce@reductions$umap@cell.embeddings,col.name = colnames(CD8_sce@reductions$umap@cell.embeddings))
CD8_sce$SubCluster<-factor(CD8_sce$SubCluster,levels=c( "CD8_C1_Tn-TCF7"   ,   "CD8_C2_Tcm-IL7R"   ,  "CD8_C3_Tem-IFNG-AS1",
                                                        "CD8_C4_Tem-GZMK"  ,   "CD8_C5_Tex-HAVCR2" ,  "CD8_C6_Trm-ZNF683"  ,
                                                        "CD8_C7-prolif"   ,    "CD8_C8-NEAT1"      ,  "CD8_C9-ISG"   , "NK_C1-CD16"        ,  "NK_C2-TNFRSF18"    , 
                                                        "NK_C3-KIT"          ,"NKT_C1-KLRD1", "Tgd_C1-TRDV2"  ))

CD8_sce$SubCluster1<-as.numeric(CD8_sce$SubCluster)
CD8_sce$TopCluster<-"CD8+ T/ILC cells"

B_sce<-readRDS("3.Cluster/13.Annotation/4.B_sub_sce_annotation.rds")
B_sce<-AddMetaData(B_sce,B_sce@reductions$umap@cell.embeddings,col.name = colnames(B_sce@reductions$umap@cell.embeddings))
B_sce$SubCluster<-factor(B_sce$SubCluster,levels=c("B_C1-TCL1A"     ,  "B_C2-CXCR4",    "B_C3-CD69"       ,       
                                                   "B_C4-RGS13"       , "B_C5-IL7R"     , "B_C6-prolif", "Plasma_C1-IGHG1" , "Plasma_C2-IGHA1",
                                                   "Plasma_C3-IGHM"  ))

B_sce$SubCluster1<-as.numeric(B_sce$SubCluster)
B_sce$TopCluster<-"B/plasma cells"


sub_sce<-merge(x=CD4_sce,y=c(B_sce,CD8_sce))

sub_sce$TopCluster<-factor(sub_sce$TopCluster,levels=c("CD4+ T cells","CD8+ T/ILC cells","B/plasma cells"))


Cluster<-unique(sort(sub_sce$SubCluster))
g.colSet1<-c("#1F77B4" ,"#FF7F0E" ,"#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#AEC7E8", "#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5", "#C49C94",
             "#F7B6D2", "#C7C7C7" ,"#DBDB8D", "#9EDAE5", "#A787C0", "#4391A9" ,"#8DC594", "#72BF5A",
             "#B49D99" ,"#F6955A" ,"#D9CE99", "#FDAE53", "#408DBF",
             "#D8AC60", "#969D62", "#92CF72",
             "#659E47" ,"#B294C7" ,"#815AA8", "#EBD57C" ,"#C48244", "#B15928", "#E94330", "#E4986B",
             "#F06C45", "#FE9E37" ,"#CCEBC5", "#A6CEE3" ,"#8EC694" ,"#FDB45D", "#CFC099", "#9B9C64",
             "#60B64D", "#EA4833", "#2F82B9", "#8D6B99" ,"#F16667", "#E62E30", "#FE982C")

Cluster<-c("CD4",
           "CD4_C1_Tn-CCR7"   ,  "CD4_C2_Tcm-AXNA1" , "CD4_C3_Tfh-CXCL13",
           "CD4_C4_Tc-GZMA"  ,  "CD4_C5_Treg-FOXP3" ,"CD4_C6-EGR1"  ,  "CD4_C7-NEAT1" ,
           "CD4_C8-TIMP1"  , "CD4_C9-ISG"  ,  "CD4_C10-prolif" ,     
           "CD8",
           "CD8_C1_Tn-TCF7"   ,   "CD8_C2_Tcm-IL7R"   ,  "CD8_C3_Tem-IFNG-AS1",
           "CD8_C4_Tem-GZMK"  ,   "CD8_C5_Tex-HAVCR2" ,  "CD8_C6_Trm-ZNF683"  ,
           "CD8_C7-prolif"   ,    "CD8_C8-NEAT1"      ,  "CD8_C9-ISG"   , 
           
           "NK",
           "NK_C1-CD16"        ,  "NK_C2-TNFRSF18"    , 
           "NK_C3-KIT" ,
           "NKT",
           "NKT_C1-KLRD1",
           "Tgd",
           "Tgd_C1-TRDV2" ,
           "B",
           "B_C1-TCL1A"     ,  "B_C2-CXCR4",    "B_C3-CD69"       ,       
           "B_C4-RGS13"       , "B_C5-IL7R"     , "B_C6-prolif" ,
           "Plasma",
           "Plasma_C1-IGHG1" , "Plasma_C2-IGHA1",
           "Plasma_C3-IGHM"        
           
)



library(ggforce)

class_avg <- sub_sce@meta.data %>%
  group_by(SubCluster1,TopCluster) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )


umap_p2<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  #geom_mark_hull(aes(label=Cluster,fill=Cluster,group=Cluster),alpha=0.1,concavity = 1)+
  stat_ellipse(data=sub_sce@meta.data,aes(group=Cluster,color=Cluster),
               alpha=0,show.legend = F,
               level = 0.95, linetype = 'dashed')+
  scale_color_manual(breaks = Cluster,
                     labels= Cluster,
                     values = g.colSet1)+
  ggtitle("Clustering of lymphocytes")+
  geom_text(aes(label = SubCluster1), data = class_avg,color="black")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))+ 
  guides(colour = guide_legend(ncol=7,override.aes = list(size=5)))+
  facet_wrap(~TopCluster,ncol=3,scales = "free")
umap_p2




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
library(monocle)
library(ggpubr)
library(patchwork)
library(cowplot)


do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  #count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


library(plyr)

combined<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")

combined$Type<-factor(combined$Type,levels=c("GN","GC" ,"GCPM","PE","PBMC" ))
meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb[which(meta.tb$Cluster %in% c("CD8","Tgd","NK","NKT","CD4","B","Plasma")),]

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Type,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[OR.dist.mtx>3]<-3
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 0

OR.dist.mtx<-t(OR.dist.mtx)

pheatmp_p2<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = F,
                               cluster_cols = T,
                               treeheight_col = 5,
                               main="Tissue preference by site"
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p2<-plot_grid(pheatmp_p2$gtable)







sub_sce<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","SD", "PR"))



meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb[which(meta.tb$Cluster %in% c("CD8","Tgd","NK","NKT","CD4","B","Plasma")),]

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Efficiency_ICB,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 10
OR.dist.mtx[OR.dist.mtx>3]<-3

OR.dist.mtx<-t(OR.dist.mtx)
OR.dist.mtx<-OR.dist.mtx[c("PR","SD","PD"),]

pheatmp_p3<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = F,
                               cluster_cols = T,
                               treeheight_col = 10,
                               main="Tissue preference by ICB efficiency for GCPM"
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p3<-plot_grid(pheatmp_p3$gtable)





##BOXplot
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
library(ggpubr)
library(easyGgplot2)
library(ggpubr)
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/3.GCPM1/")
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GN","GC","GCPM"))



##All immune cells
combined$Cluster<-factor(as.character(combined$Cluster))

##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c(    "CD4_C5_Treg-FOXP3" ,"Plasma_C1-IGHG1",
                                                        "CD4_C9-ISG" , "CD8_C6_Tex-HAVCR2"  )),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c(    "CD4_C5_Treg-FOXP3" ,"Plasma_C1-IGHG1",
                                                         "CD4_C9-ISG" , "CD8_C6_Tex-HAVCR2"   ))

color_list<-RColorBrewer::brewer.pal(9,"Blues")[c(2,4,6)]
#names(color_list)<-c( "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c( "GN","GC","GCPM"  )
sum_table$Type<-factor(sum_table$Type,levels=c("GN","GC","GCPM"  ))
select_table<-sum_table
box_p<-ggplot(select_table, aes(Type, per, fill = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
  geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format", method = "wilcox.test",
                     comparisons = list( c("GC","GN"),c("GC","GCPM")))+
  theme_classic()+
  scale_fill_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), 
        legend.position = "none",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5,face="bold"),
        strip.background = element_blank())+
  facet_wrap(~Cluster,ncol=2,scales = "free")

print(box_p)






##CD8/NK
library(patchwork)
sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")


Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD", "PR"))




gene_name<-"TCF7"
combined<-AddMetaData(combined,combined@reductions$umap@cell.embeddings,col.name = colnames(combined@reductions$umap@cell.embeddings))
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined$gene<-combined@assays$RNA@scale.data[gene_name,]
gene_exp <- combined@meta.data 

library(viridis)
CD8_prolif<-ggplot(data=gene_exp,aes(x=UMAP_1,y=UMAP_2,color=gene))+
  geom_point(size=0.2)+
  theme_classic()+
  labs(title="TCF7 in CD8+ T cells")+
  scale_colour_continuous(type ="viridis")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))+
  facet_wrap(~Efficiency_ICB,ncol=2)







##CD4
sub_sce<-readRDS("3.Cluster/13.Annotation/1.CD4_sub_sce_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD", "PR"))



gene_name<-"TCF7"
combined<-AddMetaData(combined,combined@reductions$umap@cell.embeddings,col.name = colnames(combined@reductions$umap@cell.embeddings))
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined$gene<-combined@assays$RNA@scale.data[gene_name,]
gene_exp <- combined@meta.data 

library(viridis)
CD4_prolif<-ggplot(data=gene_exp,aes(x=UMAP_1,y=UMAP_2,color=gene))+
  geom_point(size=0.2)+
  theme_classic()+
  labs(title="TCF7 in CD4+ T cells")+
  scale_colour_continuous(type ="viridis")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))+
  facet_wrap(~Efficiency_ICB,ncol=2)






sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")


Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))
combined$Topluster<-"CTL"


##All cell type
tmp_table<-data.frame(Cluster=combined$Cluster,TopCluster=combined$Topluster,Type=combined$Efficiency_ICB,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,TopCluster,Type),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,TopCluster,Type),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by=c("Type","TopCluster"))
sum_table$per<-sum_table$num/sum_table$total_num*100
#sum_table$per[which(sum_table$Cluster %in% c("Epi",
#                                             "CAF","Endo","SMC"))]<- (-(sum_table$per[which(sum_table$Cluster %in% c("Epi",
#                                                                                                                            "CAF","Endo","SMC"))]))


Cluster<-unique(sum_table$Cluster)
g.colSet <- readRDS("/home/zhengyq/data/single_cell/7.GC_gut//1.data/colSet.list.rds")
g.colSet1<-g.colSet$Set8[1:length(Cluster)]
names(g.colSet1)<-Cluster

#sum_table<-sum_table[sum_table$Cluster=="Prolif",]
sum_table$Type<-factor(sum_table$Type,levels=c("PD","PR" ))

sum_table$per<-round(sum_table$per,1)
stack_p<-ggplot(data=sum_table, aes(x=Type, y=per, fill=Cluster)) + 
  geom_bar(stat= 'identity',width = 0.8) +
  #coord_polar(theta = "y")+ 
  theme_classic()+
  scale_fill_manual(values = (RColorBrewer::brewer.pal(4,"Set3")[c(1,2,4,3)]))+
  labs(y = "Percentage",x="") +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,100),breaks = seq(0,100,20))+
  #geom_text(aes(y = per/2 + c(0, cumsum(per)[-length(per)]), x = sum(per)/100, label = per), size = 5) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        #axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank()
  )
print(stack_p)



##Expresion compare CAF
combined<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("NKT","NK","CD8"))
combined$Cluster<-factor(combined$Cluster,levels=c("CD8","NKT","NK"))

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)


data<-combined@assays$RNA@scale.data
data<-as.matrix(data)

gene_name<- c("CD3D","CD8A","CD8B","KLRD1","GNLY","FCGR3A")
mt<-data[gene_name,]
select_fpkm_matrix<-reshape2::melt(mt)
colnames(select_fpkm_matrix)<-c("gene_name","Cell_ID","gene")

Cell_info<-data.frame(Cell_ID=combined$Cell_ID,
                      Type=factor(combined$Efficiency_ICB),
                      Cluster=factor(as.character(combined$Cluster)))
select_fpkm_matrix<-merge(select_fpkm_matrix,Cell_info,by="Cell_ID")
select_fpkm_matrix$Cluster<-factor(select_fpkm_matrix$Cluster,levels=rev(c("CD8","NKT","NK")))


ridge_p1<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 1.5, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  geom_hline(yintercept = seq(from=1, to=3, by = 1),color="grey")+
  ggtitle("Expression of T and NK markers")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(3,"YlGnBu"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  #ylim(c(-1 ,5))+
  
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.line = element_blank(),
        axis.ticks.x  = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey"),
        strip.background = element_blank())+
  facet_wrap(~gene_name,scales = "free_x",ncol=3)
ridge_p1


##Expresion compare CAF
sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("NKT"))
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
NKT_exp<-VlnPlot(combined,features = c("KLRD1","FCGR3A","FGFBP2"),
        pt.size = 0,ncol = 3,
        group.by = "Efficiency_ICB") &
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y.npc = 0.9,size=3) &
  scale_fill_manual(values = c("#E64B35","#4DBBD5")) &
  theme(plot.title = element_text(size = 10,hjust = 0.5,face = "plain"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 0,vjust=0.5,hjust=0.5),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())





library(scRepertoire)

sample_info<-read.table("~/data/single_cell/3.GCPM/1.data/LHY_cohort/sample_info.txt",sep="\t",header = T)
file_list<-sample_info$TCR_ID

total_tcr<-list()
for (i in 1:length(file_list)){
  #i=8
  file_n<-file_list[i]
  sample_ID<-sample_info$Sample_ID[[i]]
  path<-paste0("~/data/single_cell/3.GCPM/1.cellranger/",file_n,"/outs/filtered_contig_annotations.csv")
  tcr <- read.table(path,sep=",",header = T)
  
  
  #tcr$barcode<-paste(sample_ID,tcr$barcode,sep="_")
  total_tcr[[length(total_tcr)+1]]<-tcr
}


total_tcr<-combineTCR(total_tcr,
                      samples = sample_info$Sample_ID)


sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")

seurat <- combineExpression(total_tcr, sub_sce, 
                            cloneCall="gene",  proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
seurat$cloneType[is.na(seurat$cloneType)]<-"None"
seurat$cloneType[seurat$cloneType=="Single (0 < X <= 1)"]<-"Single"
seurat$cloneType[seurat$cloneType=="Small (1 < X <= 5)"]<-"Small 2-5"
seurat$cloneType[seurat$cloneType=="Medium (5 < X <= 20)"]<-"Medium 6-20"
seurat$cloneType[seurat$cloneType=="Large (20 < X <= 100)"]<-"Large 21-100"
seurat$cloneType[seurat$cloneType=="Hyperexpanded (100 < X <= 500)"]<-"Hyper 101-500"
#DimPlot(seurat,group.by = "cloneType")
seurat<-AddMetaData(seurat,seurat@reductions$umap@cell.embeddings,col.name = colnames(seurat@reductions$umap@cell.embeddings))


seurat$cloneType<-factor(seurat$cloneType,levels=c("None","Single","Small 2-5","Medium 6-20","Large 21-100" ,"Hyper 101-500"))


Idents(seurat)<-seurat$Efficiency_ICB
seurat<-subset(seurat,idents=c("PR","PD"))
Idents(seurat)<-seurat$Type
seurat1<-subset(seurat,idents=c("GCPM","PBMC"))
Idents(seurat1)<-seurat1$Cluster
seurat1<-subset(seurat1,idents=c("CD8","Tgd","NKT","NK"))

CD8_p3<-ggplot(seurat1@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=cloneType))+
  geom_point(aes(color=cloneType),size=0.1) +
  scale_color_manual(values = RColorBrewer::brewer.pal(6,"BuPu"))+
  ggtitle("Clonal expansion")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line  = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        panel.border = element_rect(fill=NA,color ="grey"))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=5)))+
  facet_grid(Efficiency_ICB~Type,switch = "y",scales = "free_x")
CD8_p3

Idents(seurat)<-"Cluster"
seurat1<-subset(seurat,idents="NKT")
Idents(seurat1)<-"Type"
seurat1<-subset(seurat1,idents=c("GCPM","PBMC"))
tmp_table<-data.frame(cloneType=seurat1$cloneType,Type=seurat1$Efficiency_ICB,Cluster=seurat1$Type,num=1)
tmp_table$Type<-as.character(tmp_table$Type)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,cloneType,Type),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Cluster,Type),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by=c("Type","Cluster"))
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-factor(sum_table$Cluster,levels=c("GCPM","PBMC"))

#sum_table$Type<-factor(sum_table$Type,levels=c("Vector","CRISPRa"))
stack_p1<-ggplot(data=sum_table, aes(x=Type, y=per, fill=cloneType)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.8)+ 
  theme_classic()+
  scale_fill_manual(values = RColorBrewer::brewer.pal(6,"BuPu"))+
  labs(x = 'Type', y ="Percentage",title=paste("Clonal expansion of NKT")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,101),breaks = seq(0,100,25),labels = c("0","25%","50%","75%","100%"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank()
  )+
  facet_wrap(~Cluster,ncol=1)
print(stack_p1)


##select_TCR
sample_info<-read.table("~/data/single_cell/3.GCPM/1.data/LHY_cohort/sample_info.txt",sep="\t",header = T)
file_list<-sample_info$TCR_ID

total_tcr<-list()
for (i in 1:length(file_list)){
  #i=8
  file_n<-file_list[i]
  sample_ID<-sample_info$Sample_ID[[i]]
  path<-paste0("~/data/single_cell/3.GCPM/1.cellranger/",file_n,"/outs/filtered_contig_annotations.csv")
  tcr <- read.table(path,sep=",",header = T)
  
  
  #tcr$barcode<-paste(sample_ID,tcr$barcode,sep="_")
  total_tcr[[length(total_tcr)+1]]<-tcr
}


total_tcr<-combineTCR(total_tcr,
                      samples = sample_info$Sample_ID)


sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")
Idents(sub_sce)<-"Study"
sub_sce<-subset(sub_sce,idents=c("LHY_SYSUCC"))

seurat <- combineExpression(total_tcr, sub_sce, 
                            cloneCall="gene",  proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
seurat$cloneType[is.na(seurat$cloneType)]<-"None"
seurat$cloneType[seurat$cloneType=="Single (0 < X <= 1)"]<-"Single"
seurat$cloneType[seurat$cloneType=="Small (1 < X <= 5)"]<-"Small 2-5"
seurat$cloneType[seurat$cloneType=="Medium (5 < X <= 20)"]<-"Medium 6-20"
seurat$cloneType[seurat$cloneType=="Large (20 < X <= 100)"]<-"Large 21-100"
seurat$cloneType[seurat$cloneType=="Hyperexpanded (100 < X <= 500)"]<-"Hyper 101-500"
#DimPlot(seurat,group.by = "cloneType")
seurat<-AddMetaData(seurat,seurat@reductions$umap@cell.embeddings,col.name = colnames(seurat@reductions$umap@cell.embeddings))


seurat$cloneType<-factor(seurat$cloneType,levels=c("None","Single","Small 2-5","Medium 6-20","Large 21-100" ,"Hyper 101-500"))

Idents(seurat)<-seurat$Cluster
NKT_seu<-subset(seurat,idents=c("NKT"))
seurat1<-NKT_seu
Idents(seurat1)<-"Type"
seurat1<-subset(seurat1,idents=c("GCPM","PBMC","GC","GN"))
Idents(seurat1)<-"Sample_ID"
seurat1<-subset(seurat1,idents=c("DD37","DD38","DD39","DD40","DD52","DD53","DD54","DD55","DD58","DD59","DD60","DD61","DD71","DD72","DD73","DD74"))
select_barcodes<-seurat1$Cell_ID
sample_info<-read.table("~/data/single_cell/3.GCPM/1.data/LHY_cohort/sample_info.txt",sep="\t",header = T)
file_list<-sample_info$TCR_ID

total_tcr<-list()
ID_list<-c()
type_list<-c()
for (i in 1:length(file_list)){
  #i=8
  file_n<-file_list[i]
  sample_ID<-sample_info$Sample_ID[[i]]
  Type=sample_info$Type[i]
  if (Type %in% c("GCPM","PBMC","GC","GN")){
    path<-paste0("~/data/single_cell/3.GCPM/1.cellranger/",file_n,"/outs/filtered_contig_annotations.csv")
    tcr <- read.table(path,sep=",",header = T)
    tcr$barcode1<-paste(sample_ID,tcr$barcode,sep="_")
    tcr<-tcr[which(tcr$barcode1 %in% select_barcodes),]
    if (nrow(tcr)>0){
      total_tcr[[length(total_tcr)+1]]<-tcr
      ID_list<-c(ID_list,sample_info$Patiant_ID[[i]])
      type_list<-c(type_list,Type)
    }
    
  }
  
}


total_tcr<-combineTCR(total_tcr,
                      samples = ID_list,
                      ID = type_list,
)


subset1<-subsetContig(total_tcr,name="sample",variables="Pt12")
library(limma)
names(subset1)<-strsplit2(names(subset1),"_")[,2]
mat_melt<-clonalOverlap(subset1,cloneCall = "aa",exportTable = T,
                        method = "overlap")
mat_melt <- suppressMessages(melt(as.matrix(mat_melt[,-ncol(mat_melt)])))
#mat_melt$value<-round(mat_melt$value,1)
mean_value <- mean(na.omit(mat_melt[,"value"]))

pt12<-  ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + 
  geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], 
            fill = NA, 
            lwd= 0.25, 
            color = "grey") +
  xlim(c("GC","GCPM","GN"))+
  ylim(c("GCPM","GN","PBMC"))+
  geom_text(aes(label = round(value, digits = 2), 
                color = ifelse(value <= mean_value,
                               "black", "black")), 
            size=3,
            na.rm = TRUE) +
  #scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
  scale_color_identity() +
  theme_classic() + 
  theme(axis.title = element_blank())+ 
  scale_fill_gradientn(colours = c("#F1EEF6","#BDC9E1","#FC8D59","#D7301F"),
                       limits=c(0,0.7),
                       
                       na.value = "white")+
  ggtitle("Pat12")+
  theme(axis.text = element_text(size = 10),  
        #axis.text.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank()
  )

pt12

subset1<-subsetContig(total_tcr,name="sample",variables="Pt13")
library(limma)
names(subset1)<-strsplit2(names(subset1),"_")[,2]
mat_melt<-clonalOverlap(subset1,cloneCall = "aa",exportTable = T,
                        method = "overlap")
mat_melt <- suppressMessages(melt(as.matrix(mat_melt[,-ncol(mat_melt)])))
#mat_melt$value<-round(mat_melt$value,1)
mean_value <- mean(na.omit(mat_melt[,"value"]))

pt13<-  ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + 
  geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], 
            fill = NA, 
            lwd= 0.25, 
            color = "grey") +
  xlim(c("GC","GCPM","GN"))+
  ylim(c("GCPM","GN","PBMC"))+
  geom_text(aes(label = round(value, digits = 2), 
                color = ifelse(value <= mean_value,
                               "black", "black")), 
            size=3,
            na.rm = TRUE) +
  #scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
  scale_color_identity() +
  theme_classic() + 
  theme(axis.title = element_blank())+ 
  scale_fill_gradientn(colours = c("#F1EEF6","#BDC9E1","#FC8D59","#D7301F"),
                       limits=c(0,0.7),
                       
                       na.value = "white")+
  ggtitle("Pat13")+
  theme(axis.text = element_text(size = 10),  
        #axis.text.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank()
  )

pt13


subset1<-subsetContig(total_tcr,name="sample",variables="Pt14")
library(limma)
names(subset1)<-strsplit2(names(subset1),"_")[,2]
mat_melt<-clonalOverlap(subset1,cloneCall = "aa",exportTable = T,
                        method = "overlap")
mat_melt <- suppressMessages(melt(as.matrix(mat_melt[,-ncol(mat_melt)])))
#mat_melt$value<-round(mat_melt$value,1)
mean_value <- mean(na.omit(mat_melt[,"value"]))

pt14<-  ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + 
  geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], 
            fill = NA, 
            lwd= 0.25, 
            color = "grey") +
  xlim(c("GC","GCPM","GN"))+
  ylim(c("GCPM","GN","PBMC"))+
  geom_text(aes(label = round(value, digits = 2), 
                color = ifelse(value <= mean_value,
                               "black", "black")), 
            size=3,
            na.rm = TRUE) +
  #scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
  scale_color_identity() +
  theme_classic() + 
  theme(axis.title = element_blank())+ 
  scale_fill_gradientn(colours = c("#F1EEF6","#BDC9E1","#FC8D59","#D7301F"),
                       limits=c(0,0.7),
                       
                       na.value = "white")+
  ggtitle("Pat14")+
  theme(axis.text = element_text(size = 10),  
        #axis.text.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank()
  )

pt14

subset1<-subsetContig(total_tcr,name="sample",variables="Pt15")
library(limma)
names(subset1)<-strsplit2(names(subset1),"_")[,2]
mat_melt<-clonalOverlap(subset1,cloneCall = "aa",exportTable = T,
                        method = "overlap")
mat_melt <- suppressMessages(melt(as.matrix(mat_melt[,-ncol(mat_melt)])))
#mat_melt$value<-round(mat_melt$value,1)
mean_value <- mean(na.omit(mat_melt[,"value"]))

pt15<-  ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() + 
  geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], 
            fill = NA, 
            lwd= 0.25, 
            color = "grey") +
  xlim(c("GC","GCPM","GN"))+
  ylim(c("GCPM","GN","PBMC"))+
  geom_text(aes(label = round(value, digits = 2), 
                color = ifelse(value <= mean_value,
                               "black", "black")), 
            size=3,
            na.rm = TRUE) +
  #scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
  scale_color_identity() +
  theme_classic() + 
  theme(axis.title = element_blank())+ 
  scale_fill_gradientn(colours = c("#F1EEF6","#BDC9E1","#FC8D59","#D7301F"),
                       limits=c(0,0.7),
                       
                       na.value = "white")+
  ggtitle("Pat15")+
  theme(axis.text = element_text(size = 10),  
        #axis.text.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank()
  )

pt15

pt12<-pt12+theme(axis.text.x = element_blank(),
                 legend.position = "none")
pt13<-pt13+theme(axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 legend.position = "right")
pt14<-pt14+theme(legend.position = "none")
pt15<-pt15+theme(axis.text.y = element_blank(),
                 legend.position = "none")
                 
tcr_overlap<-pt12+pt13+pt14+pt15+plot_layout(ncol=2)

  



library(patchwork)
total_p1<-umap_p2
total_p2<-ggarrange(stack_p,NKT_exp,ncol=2,nrow=1,widths = c(1,1.5))
total_p2<-ggarrange(total_p2,ridge_p1,ncol=1,nrow=2,heights = c(1,1.5))
total_p2<-ggarrange(total_p2,NULL,ncol=2,nrow=1)


total_p3<-ggarrange(tcr_overlap,CD8_p3,stack_p1,ncol=3,widths = c(1.4,1.4,1.2))
#total_p3<-ggarrange(NKT_exp,tcr_overlap,CD8_p3,ncol=3,widths = c(1.2,1.4,1.4))

total_p<-ggarrange(total_p1,plot_spacer()+pheatmp_p3+plot_layout(ncol=2,widths=c(0.2,10)),total_p2,total_p3,ncol=1,heights = c(6,3.2,4,4))
ggsave("14.Figure_plot/Figure2_1.pdf",total_p,width = 12,height = 17.2)


