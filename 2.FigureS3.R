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




sub_sce<-readRDS("3.Cluster/13.Annotation/4.B_sub_sce_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE","PBMC" ))
dir.create("3.Cluster/13.plot/3.myeloid/")
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|Mϕ_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(  "B_C1-TCL1A"     ,  "B_C2-CXCR4",    "B_C3-CD69"       ,       
                                                         "B_C4-RGS13"       , "B_C5-IL7R"     , "B_C6-prolif", "Plasma_C1-IGHG1" , "Plasma_C2-IGHA1",
                                                         "Plasma_C3-IGHM"    ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c("B_C1"     ,  "B_C2",    "B_C3"       ,       
                  "B_C4"       , "B_C5"     , "B_C6", "Plasma_C1" , "Plasma_C2",
                  "Plasma_C3" ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tsne")<-Embeddings(object = sub_sce, reduction = "tsne")

Cluster<-unique(sort(sce$SubCluster))
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster
#g.colSet1<-list("SubCluster"=g.colSet1)

library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(SubCluster1) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )


B_umap1<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  #geom_mark_hull(aes(label=Cluster,fill=Cluster,group=Cluster),alpha=0.1,concavity = 1)+
  stat_ellipse(data=sub_sce@meta.data,aes(group=Cluster,color=Cluster),
               alpha=0.5,show.legend = F,
               level = 0.95, linetype = 'dashed')+
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c(  "B_C1-TCL1A"     ,  "B_C2-CXCR4",    "B_C3-CD69"       ,       
                                 "B_C4-RGS13"       , "B_C5-IL7R"     , "B_C6-prolif", "Plasma_C1-IGHG1" , "Plasma_C2-IGHA1",
                                 "Plasma_C3-IGHM"   ),
                     values = g.colSet1)+
  ggtitle("Clustering of B/plasma cells")+
  geom_text_repel(aes(label = SubCluster1), data = class_avg,color="black")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8))+ 
  guides(colour = guide_legend(override.aes = list(size=2),ncol = 1))


B_umap1








sub_sce<-readRDS("3.Cluster/13.Annotation/1.CD4_sub_sce_annotation.rds")
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)
sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(  "CD4_C1_Tn-CCR7"   ,  "CD4_C2_Tcm-AXNA1" , "CD4_C3_Tfh-CXCL13",
                                                         "CD4_C4_Tc-GZMA"  ,  "CD4_C5_Treg-FOXP3" ,"CD4_C6-EGR1"  ,  "CD4_C7-NEAT1" ,
                                                         "CD4_C8-TIMP1"  , "CD4_C9-ISG"  ,  "CD4_C10-prolif"    ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c(  "CD4_C1_Tn-CCR7"   ,  "CD4_C2_Tcm-AXNA1" , "CD4_C3_Tfh-CXCL13",
                    "CD4_C4_Tc-GZMA"  ,  "CD4_C5_Treg-FOXP3" ,"CD4_C6-EGR1"  ,  "CD4_C7-NEAT1" ,
                    "CD4_C8-TIMP1"  , "CD4_C9-ISG"  ,  "CD4_C10-prolif"   )


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster1
CD4_p<-DotPlot(sub_sce, features = rev(c("LEF1"  ,   "SELL"   ,  "TCF7"  ,  "CCR7"    , "GADD45B",  "NR4A2"  ,    "ANXA1" ,
                                         "IL7R"  ,   "TXNIP"  ,   "LMNA"   ,  "CXCL13" ,  "PDCD1"  ,  "ICOS"  ,  "CXCR5"  ,  "TOX2"  ,
                                         "LAG3"   , "IFNG", "GZMK"    , "KLRG1"   ,  "EGR1","TNF",  
                                         "FOXP3" ,   "IL2RA" ,   "CTLA4"   ,  "TNFRSF4",  "NEAT1"   ,  "HNRNPH1", "PDE3B" ,  "TIMP1","S100A6","GIMAP7",
                                         "IFI6"  ,   "ISG15"  , "IFIT3"  ,  "STMN1" ,   "HMGB2"   , "TUBB"   ,  "RTKN2" )),
               cols = c("grey","red2"),dot.scale=5)+ggplot2::coord_flip()+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  theme(legend.position = "none")+
  ggtitle("Markers of CD4+ T cell clusters")+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))

CD4_p




sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)
sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(   "CD8_C1_Tn-TCF7"   ,   "CD8_C2_Tcm-IL7R"   ,  "CD8_C3_Tem-IFNG-AS1",
                                                          "CD8_C4_Tem-GZMK"  ,   "CD8_C5_Tex-HAVCR2" ,  "CD8_C6_Trm-ZNF683"  ,
                                                          "CD8_C7-prolif"   ,    "CD8_C8-NEAT1"      ,  "CD8_C9-ISG"   , "NK_C1-CD16"        ,  "NK_C2-TNFRSF18"    , 
                                                          "NK_C3-KIT"          ,"NKT_C1-KLRD1", "Tgd_C1-TRDV2"  ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c(    "CD8_C1_Tn-TCF7"   ,   "CD8_C2_Tcm-IL7R"   ,  "CD8_C3_Tem-IFNG-AS1",
                      "CD8_C4_Tem-GZMK"  ,   "CD8_C5_Tex-HAVCR2" ,  "CD8_C6_Trm-ZNF683"  ,
                      "CD8_C7-prolif"   ,    "CD8_C8-NEAT1"      ,  "CD8_C9-ISG"   , "NK_C1-CD16"        ,  "NK_C2-TNFRSF18"    , 
                      "NK_C3-KIT"          ,"NKT_C1-KLRD1", "Tgd_C1-TRDV2"    )


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster1
CD8_p<-DotPlot(sub_sce, features = rev(c( "CD3D","CD8A","TCF7"    , "SELL"   ,  "LTB"   ,     "IL7R"  ,   "CCR7"   ,  "CXCR4" , 
                                          "GZMK"    , "IFNG-AS1"  ,  "TNFSF9"  ,   "CCL4"   ,  "CCL4L2" ,  "HLA-DRA" ,
                                          "CXCL13" ,  "HAVCR2" ,  "PDCD1"  ,   "ZNF683"  ,  "HOPX",  "STMN1"  ,  "HMGB2"  ,  "MKI67"      ,  "NEAT1"   ,  
                                          "ZEB1"   ,"ISG15"  , 
                                          "IFI6"    , "MX1"    ,   "LPP" ,"GNLY"  ,   "FCGR3A"  ,  "FGFBP2" ,"GZMH"   ,  "KRT86"  , "TNFRSF18", "LST1"   ,   
                                          "KIT" ,      
                                          "CMC1"   ,  "KLRG1"  ,  "TRDV2",  "TRGV9")),
               cols = c("grey","red2"),dot.scale=4)+ggplot2::coord_flip()+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  ggtitle("Markers of CD8+ T/NK cell clusters")+
  theme(legend.position = "none")+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))
CD8_p


sub_sce<-readRDS("3.Cluster/13.Annotation/4.B_sub_sce_annotation.rds")
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)
sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c( "B_C1-TCL1A"     ,  "B_C2-CXCR4",    "B_C3-CD69"       ,       
                                                        "B_C4-RGS13"       , "B_C5-IL7R"     , "B_C6-prolif", "Plasma_C1-IGHG1" , "Plasma_C2-IGHA1",
                                                        "Plasma_C3-IGHM"  ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c("B_C1-TCL1A"     ,  "B_C2-CXCR4",    "B_C3-CD69"       ,       
                  "B_C4-RGS13"       , "B_C5-IL7R"     , "B_C6-prolif", "Plasma_C1-IGHG1" , "Plasma_C2-IGHA1",
                  "Plasma_C3-IGHM")


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster1
B_p<-DotPlot(sub_sce, features = rev(c("CD83",  "TCL1A"   ,  "HLA-DRA" , "CXCR4"     , "CD69"    ,   "GPR183" , "HLA-DRB1",
                                       "RGS13","TOX","SOX5","IL7R"     , "IL32"     , "CCL5"   ,   "STMN1"   ,  "HMGB2"    , "TUBA1B"  , 
                                         "IGHG1"    , "IGHG2"  , "IGHG4"  ,"MZB1",
                                         "IGHM","IGHA1","IGHA2" )),
               cols = c("grey","purple"),dot.scale=5)+ggplot2::coord_flip()+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)+
  ggtitle("Markers of B cell clusters")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))
B_p




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

meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb[which(meta.tb$Cluster %in% c("CD8","MAIT","Tgd","NK","NKT","CD4","B","Plasma")),]

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
                               treeheight_col = 10,
                               main="Tissue preference by site"
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p2<-plot_grid(pheatmp_p2$gtable)






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
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GN","GC","GCPM","PE","PBMC"))



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

sum_table<-sum_table[which(sum_table$Cluster %in% c(    "CD4_C5_Treg-FOXP3","CD4_C7-NEAT1" ,
                                                        "CD4_C9-ISG" ,"CD8_C5_Tex-HAVCR2", "B_C1-TCL1A", "Plasma_C1-IGHG1",  "NK_C1-CD16","NKT_C1-KLRD1")),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c(   "CD4_C5_Treg-FOXP3", "CD4_C7-NEAT1" , 
                                                         "CD4_C9-ISG" , "CD8_C5_Tex-HAVCR2" , "B_C1-TCL1A","Plasma_C1-IGHG1","NK_C1-CD16","NKT_C1-KLRD1"  ))

color_list<-RColorBrewer::brewer.pal(5,"BuPu")

select_table<-sum_table
select_table$Type<-factor(select_table$Type,levels=c("GN","GC","GCPM","PE","PBMC"))

box_p<-ggplot(select_table, aes(Type, per, fill = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
  geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format", method = "wilcox.test",
                     comparisons = list( c("GC","GN"),c("GC","GCPM"),c("GN","GCPM"),c("PBMC","GCPM")))+
  
  theme_classic()+
  scale_fill_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), 
        legend.position = "none",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle=45,vjust=1,hjust=1),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5,face="bold"),
        strip.background = element_blank())+
  facet_wrap(~Cluster,ncol=10,scales = "free")

print(box_p)




##vocano
vocano_plot = function(sorted_matrix, Sample_1 = "A", Sample_2 = "B",
                       lfc = 0, pval = 0.05,gene_list=NA,title="DEG analyses"){
  par(mar = c(5, 6, 5, 5))
  sorted_matrix$p_val_adj[which(sorted_matrix$p_val_adj<10^-200)]<-10^-200
  tab = data.frame(logFC = as.numeric(as.character(sorted_matrix$avg_log2FC)), 
                   negLogPval = -log10(as.numeric(as.character(sorted_matrix$p_val_adj))))
  rownames(tab)=rownames(sorted_matrix)
  
  tab <- na.omit(tab)
  tab <- tab[(abs(tab$logFC) > lfc &  tab$negLogPval > -log10(pval)),]
  
  nosigGene = rownames(tab)[(abs(tab$logFC) <= lfc | tab$negLogPval <= -log10(pval))]
  sigGenes_up = rownames(tab)[(tab$logFC > lfc & tab$negLogPval > -log10(pval))]
  sigGenes_down = rownames(tab)[(tab$logFC < -lfc & tab$negLogPval > -log10(pval))]
  draw_up = rownames(tab)[(tab$logFC > lfc & tab$negLogPval > -log10(pval))]
  draw_down = rownames(tab)[(tab$logFC < -lfc & tab$negLogPval > -log10(pval))]
  up_count = length(sigGenes_up)
  down_count = length(sigGenes_down)
  nosig_count = length(nosigGene)
  gap = max(tab$logFC)/50
  FCrange=ceiling(max(abs(tab$logFC)))
  tab[sigGenes_up,"SigGenes"]=paste("1.up:",up_count,sep="")
  tab[sigGenes_down,"SigGenes"]=paste("2.down:",down_count,sep="")
  tab[nosigGene,"SigGenes"]=paste("3.noSig",sep="")
  tab$name=rownames(tab)
  options(stringsAsFactors = FALSE)  ### NOTICE!!!
  DF=data.frame(name=as.character(tab$name),SigGenes=as.factor(tab$SigGenes),logFC=tab$logFC,negLogPval=tab$negLogPval)
  rownames(DF)=rownames(tab)
  #DF <- DF[sort(DF$logFC,index.return=TRUE, decreasing = TRUE)$ix,]
  tophit=DF[c(draw_up[1:10],draw_down[max((length(draw_down)-10), 0):length(draw_down)]),]
  xmax <- ceiling(max(abs(DF$logFC)))
  ymax <- ceiling(max(abs(DF$negLogPval)))*1.1
  
  
  DF$value<-abs(DF$logFC)
  DF<-DF[order(DF$logFC,decreasing = T),]
  DF$Gene<-""
  if (nrow(DF)>20){
    DF$Gene[1:20]<-rownames(DF)[1:20]
  }else{
    DF$Gene<-rownames(DF)
  }
  DF<-DF[order(DF$logFC,decreasing = F),]
  if (nrow(DF)>20){
    DF$Gene[1:20]<-rownames(DF)[1:20]
  }else{
    DF$Gene<-rownames(DF)
  }
  DF$Gene<-""
  DF$Gene[rownames(DF) %in% gene_list]<-rownames(DF)[rownames(DF) %in% gene_list]
  #DF$Gene<-rownames(DF)
  print(c(head(DF,20),tail(DF,20)))
  p <- ggplot(DF, aes(x = logFC, y = negLogPval, label=DF$Gene)) +
    geom_point(aes(color = SigGenes))+ 
    geom_text_repel(max.overlaps = 2000)+
    xlim(-xmax,xmax) + ylim(0,ymax) +
    scale_color_manual(values = c("#B31B21", "#1465AC","grey")) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom") +
    xlab("Log2FoldChange")+ylab(paste("-log10"," FDR",sep=""))+
    geom_vline(aes(xintercept=-lfc),colour="darkgrey", linetype="dashed")+
    geom_vline(aes(xintercept=lfc),colour="darkgrey", linetype="dashed") +
    geom_hline(aes(yintercept=-log10(pval)),colour="darkgrey", linetype="dashed")+
    ggtitle(title,subtitle = paste(Sample_2,Sample_1,sep=paste(rep(" ",15),collapse=""))) +
    annotate("text", x=-xmax*0.9, y=-log10(pval), label= paste("FDR"," < ",pval,sep=""))+
    #annotate("text", x=0, y=-log10(pval), label= "2fold")+
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

libs <- c("limma","ggpubr","ggrepel","reshape2","cluster","splines","ggplot2","gridExtra")
loaded <- sapply(libs, library, character.only = T)
fc = 2
lfc = log2(fc)
pval = 0.05
sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("NKT"))

Idents(combined)<-combined$Efficiency_ICB
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
C1<-"PR"
C2<-"PD"

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
markers1<-markers$gene[which(markers$cluster==C1)]
markers$avg_log2FC[which(markers$cluster==C2)]<-(-(markers$avg_log2FC[which(markers$cluster==C2)]))

gene_list<-c("HLA-DRB5","CD69","DUSP2","CD8B","CXCR4","CX3CR1",
             "VAV3","KLRF1","GNLY","KLRD1","FCGR3A","FGFBP2","TYROBP" )
vocano_p1<-vocano_plot(markers, Sample_1 = C1, Sample_2 = C2, lfc = 0, pval = pval,
            gene_list=gene_list,title="NKT cells by ICB (GCPM)")+
  theme(legend.position = "none")+
  ylim(c(0,40))+
  theme(axis.text = element_text(size = 10), 
        legend.position = "none",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        plot.subtitle = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())


sub_sce<-readRDS("3.Cluster/13.Annotation/2.CD8_sub_sce_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="PBMC")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("NKT"))

Idents(combined)<-combined$Efficiency_ICB
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
C1<-"PR"
C2<-"PD"

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
markers1<-markers$gene[which(markers$cluster==C1)]
markers$avg_log2FC[which(markers$cluster==C2)]<-(-(markers$avg_log2FC[which(markers$cluster==C2)]))

gene_list<-c("HLA-DRB5","CD69","DUSP2","CD8B","CXCR4","KIR3DL2",
             "IFITM2","IFITM1","GNLY","KLRD1","FCGR3A","FGFBP2","TYROBP" )
vocano_p2<-vocano_plot(markers, Sample_1 = C1, Sample_2 = C2, lfc = 0, pval = pval,
                       gene_list=gene_list,title="Differential exression by ICB (NKT cells, PBMC)")

vocano_p2





##GO_BP CD8
Enrich_type<-"KEGG"
Type="All"
Select_cluster_list<-c("CD8","NKT","NK")
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/3.NKT//",cluster,"/",Enrich_type,"_enrich.xls",sep="")
  enrich_table<-read.table(path,sep="\t",header = T,row.names = 1)
  enrich_table$Cluster<-cluster
  total_table<-rbind(total_table,enrich_table)
  if (nrow(enrich_table)>5){
    enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
    enrich_table<-head(enrich_table,5)
  }
  
  enriched_pathway<-c(enriched_pathway,enrich_table$Description)
  
}
enriched_pathway<-unique(enriched_pathway)


total_table1<-total_table[which(total_table$Cluster %in% Select_cluster_list&total_table$Description %in% enriched_pathway),]
total_table1$Description<-factor(total_table1$Description,levels=rev(unique(total_table1$Description)))
total_table1$Cluster<-factor(total_table1$Cluster,levels=Select_cluster_list)


S0<- ggplot(total_table1, aes(x= Cluster, y=Description, size=Count, color=p.adjust, group=Cluster)) + 
  geom_point(alpha = 0.8) + 
  theme_classic()+
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")+
  scale_size(range = c(2, 6))+
  ggtitle(Enrich_type)+
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) 

S0 



##GO_BP CD8
Enrich_type<-"KEGG"
Type="All"
Select_cluster_list<-c("PD","PR")
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/4.NKT/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
  enrich_table<-read.table(path,sep="\t",header = T,row.names = 1)
  enrich_table$Cluster<-cluster
  total_table<-rbind(total_table,enrich_table)
  if (nrow(enrich_table)>5){
    enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
    enrich_table<-head(enrich_table,5)
  }
  
  enriched_pathway<-c(enriched_pathway,enrich_table$Description)
  
}
enriched_pathway<-unique(enriched_pathway)


total_table1<-total_table[which(total_table$Cluster %in% Select_cluster_list&total_table$Description %in% enriched_pathway),]
total_table1$Description<-factor(total_table1$Description,levels=rev(unique(total_table1$Description)))
total_table1$Cluster<-factor(total_table1$Cluster,levels=Select_cluster_list)


S1<- ggplot(total_table1, aes(x= Cluster, y=Description, size=Count, color=p.adjust, group=Cluster)) + 
  geom_point(alpha = 0.8) + 
  theme_classic()+
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")+
  scale_size(range = c(2, 6))+
  ggtitle(Enrich_type)+
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) 

S1 


##GO_BP CD8
Enrich_type<-"GO_BP"
Type="All"
Select_cluster_list<-c("PD","PR")
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/4.NKT/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
  enrich_table<-read.table(path,sep="\t",header = T,row.names = 1)
  enrich_table$Cluster<-cluster
  total_table<-rbind(total_table,enrich_table)
  if (nrow(enrich_table)>5){
    enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
    enrich_table<-head(enrich_table,5)
  }
  
  enriched_pathway<-c(enriched_pathway,enrich_table$Description)
  
}
enriched_pathway<-unique(enriched_pathway)


total_table1<-total_table[which(total_table$Cluster %in% Select_cluster_list&total_table$Description %in% enriched_pathway),]
total_table1$Description<-factor(total_table1$Description,levels=rev(unique(total_table1$Description)))
total_table1$Cluster<-factor(total_table1$Cluster,levels=Select_cluster_list)


S2<- ggplot(total_table1, aes(x= Cluster, y=Description, size=Count, color=p.adjust, group=Cluster)) + 
  geom_point(alpha = 0.8) + 
  theme_classic()+
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")+
  scale_size(range = c(2, 6))+
  ggtitle(Enrich_type)+
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) 

S2 




total_p2<-CD4_p+CD8_p+B_p+plot_layout(ncol=3,widths = c(16,20,11))
total_p3<-S0+vocano_p1+S1+plot_layout(ncol=3,widths = c(0.5,1.1,0.4))

total_p<-ggarrange(total_p2,pheatmp_p2,box_p,total_p3,heights = c(7,3.5,2.5,3.5),ncol = 1,nrow=4)
pdf("14.Figure_plot/Figure2_S3.pdf",width = 12,height = 16.5)
print(total_p)
dev.off()
