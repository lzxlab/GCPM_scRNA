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
setwd("/home/zhengyq/data/single_cell/3.GCPM/")


sub_sce<-readRDS("3.Cluster/5.SubAnnotation/7.Endo/sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("GN","GC" ,"GCPM","PE" ))
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE" ))
dir.create("3.Cluster/13.plot/3.myeloid/")
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(    "Endo_C1_angio-INSR"      ,   "Endo_C2_angio-PGF"     ,   "Endo_C3_immature-APLNR",  
                                                           "Endo_C4_arterial-CXCL12" , "Endo_C5_venous-SELE"     , "Endo_C6_venous-ACKR1"    ,
                                                           "Endo_C7_capillary-CA4"  ,  "Endo_C8_scavenging-C1QC" , "Endo_C9_lymphatics-PROX1",
                                                           "Endo_C10-MYL9"        ,    "Endo_C11-KRT8"        ,   
                                                           "Endo_C12-other"      ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "Endo_C1"      ,   "Endo_C2"     ,   "Endo_C3",  
                   "Endo_C4" , "Endo_C5"     , "Endo_C6"    ,
                   "Endo_C7"  ,  "Endo_C8" , "Endo_C9",
                   "Endo_C10"        ,    "Endo_C11"        ,   
                   "Endo_C12"   ) 



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


umap_p1<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  #geom_mark_hull(aes(label=Cluster,fill=Cluster,group=Cluster),alpha=0.1,concavity = 1)+
  stat_ellipse(data=sub_sce@meta.data,aes(group=Cluster,color=Cluster),
               alpha=0.5,show.legend = F,
               level = 0.95, linetype = 'dashed')+
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c(   "C1_angio-INSR"      ,   "C2_angio-PGF"     ,   "C3_immature-APLNR",  
                                  "C4_arterial-CXCL12" , "C5_venous-SELE"     , "C6_venous-ACKR1"    ,
                                  "C7_capillary-CA4"  ,  "C8_scavenging-C1QC" , "C9_lymphatics-PROX1",
                                  "C10-MYL9"        ,    "C11-KRT8"        ,   
                                  "C12-other"     ),
                     values = g.colSet1)+
  ggtitle("Clustering of Endothelial cells")+
  geom_text(aes(label = SubCluster1), data = class_avg,color="black")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(override.aes = list(size=3)))


umap_p1




#Merged plot
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster
g.colSet1<-list("SubCluster"=g.colSet1)


umap_p2 <- sscVis::ssc.plot.tsne(sce, columns = "SubCluster", splitBy = "Type",
                            reduced.name = "umap",
                            
                            colSet=g.colSet1,size=0.1,label=0,
                            par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                            #vector.friendly=T,
                            #
                            fun.extra=function(p){ p + guides(color=guide_legend(ncol=1, override.aes=list(size=4))) +
                                ggplot2::facet_wrap(~splitBy, ncol = 5)+
                                theme_classic()+
                                theme(legend.position = "none",
                                      plot.title = element_blank())},
                            par.geom_point = list(scale=0.8),
                            par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                            base_aspect_ratio = 1.35,
                            p.ncol = 1)
print(umap_p2)




sub_sce<-readRDS("3.Cluster/5.SubAnnotation/7.Endo//sub_sce_annotation.rds")

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(   "Endo_C1_angio-INSR"      ,   "Endo_C2_angio-PGF"     ,   "Endo_C3_immature-APLNR",  
                                                          "Endo_C4_arterial-CXCL12" , "Endo_C5_venous-SELE"     , "Endo_C6_venous-ACKR1"    ,
                                                          "Endo_C7_capillary-CA4"  ,  "Endo_C8_scavenging-C1QC" , "Endo_C9_lymphatics-PROX1",
                                                          "Endo_C10-MYL9"        ,    "Endo_C11-KRT8"        ,   
                                                          "Endo_C12-other"      ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "Endo_C1_angio-INSR"      ,   "Endo_C2_angio-PGF"     ,   "Endo_C3_immature-APLNR",  
                   "Endo_C4_arterial-CXCL12" , "Endo_C5_venous-SELE"     , "Endo_C6_venous-ACKR1"    ,
                   "Endo_C7_capillary-CA4"  ,  "Endo_C8_scavenging-C1QC" , "Endo_C9_lymphatics-PROX1",
                   "Endo_C10-MYL9"        ,    "Endo_C11-KRT8"        ,   
                   "Endo_C12-other"   ) 


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster1
dot_p<-DotPlot(sub_sce, features = rev(c("INSR"  , "FLT1"    , "KDR"  , "PGF" ,     "ESM1"    , "ANGPT2"   , "APLNR" , "IGFBP5" , "PRCP", 
                                         "GJA4"  ,   "CXCL12",   "FBLN5",   "SELE" ,  "CSF3"  ,   "IL6"    ,  "ACKR1" ,  
                                         "CLU"   ,   "CCL14"  , "CA4" ,  "RBP7"  ,   "TIMP3" ,     "FCER1G"  ,   "C1QC"   ,   "C1QB" ,   
                                         "PROX1"   ,"CCL21"  ,  "TFF3"   ,   "RGS5"   ,  "MYL9" ,"ACTA2"  ,  "KRT8"  ,  "KRT5" ,   "CXCL11" ,  "ISG15"  ,  "CXCL10")),
               cols = c("#007896FF","#FDE725FF"),dot.scale=5)+ggplot2::coord_flip()+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(legend.position = "none")+
  scale_color_viridis(discrete=F, option = "B", begin = 0, end=1, direction=1)+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))

dot_p



##All Cells
combined<-readRDS("3.Cluster/5.SubAnnotation/7.Endo/sub_sce_annotation.rds")


tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,Type),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="Type")
sum_table$per<-sum_table$num/sum_table$total_num*100
#sum_table$per[which(sum_table$Cluster %in% c("Epi",
#                                             "CAF","Endo","SMC"))]<- (-(sum_table$per[which(sum_table$Cluster %in% c("Epi",
#                                                                                                                            "CAF","Endo","SMC"))]))


Cluster<-levels(sum_table$Cluster)
Cluster<-unique(sort(sce$SubCluster))
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster

sum_table$Type<-factor(sum_table$Type,levels=c("GN" , "GC","GCPM","PE","PBMC"))
stack_p<-ggplot(data=sum_table, aes(x=Type, y=per, fill=Cluster)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.8)+ 
  theme_classic()+
  scale_fill_manual(values = g.colSet1)+
  labs(x = 'Type', y ="Percentage",title=paste("")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 10)) +
  #scale_x_continuous(breaks = seq(0,15,5),labels = c("0","5","10","15+"))+
  scale_y_continuous(limits = c(0,101),breaks = seq(0,100,25),labels = c("0","25%","50%","75%","100%"))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust=1),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
  )
print(stack_p)










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
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


library(plyr)

combined<-readRDS("3.Cluster/12.Annotation/combined_annotation.rds")

Idents(combined)<-combined$Type
combined<-subset(combined,idents=c("GN","GC","GCPM","PE"))

meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb[which(meta.tb$Cluster %in% c("Endo")),]

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Type,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[OR.dist.mtx>3]<-3
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 0

#OR.dist.mtx<-t(OR.dist.mtx)
OR.dist.mtx<-OR.dist.mtx[,c("GN","GC","GCPM","PE")]
pheatmp_p2<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = T,
                               cluster_cols = F,
                               
                               main="Tissue distribution",
                               angle_col = 45
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
setwd("/home/zhengyq/data/single_cell/3.GCPM/")
combined1<-readRDS("3.Cluster/9.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GN","GC","GCPM"))
if (F){
  Idents(combined)<-combined$Cluster
  combined<-subset(combined,idents=c( "B"    ,    "Plasma"  ,     "CD4_conv"   ,     "CD4_Treg"       , "CD8"     ,    "MAIT" ,  "NK"   ,  "Tgd",   
                                      "Mast"   , "Neutro" ,  "DC"  ,              
                                      
                                      "TAM"   , "Mono"  ,    "Endo"   ,        
                                      "iCAF" ,  "mCAF" ,   "vCAF",  "Chief_cells"  , "Mucoid_cells"   ,    "Epi"     ))
}



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

sum_table<-sum_table[which(sum_table$Cluster %in% c(  "Endo_C1_angio-INSR"      ,   "Endo_C2_angio-PGF"     ,   "Endo_C3_immature-APLNR",  
                                                      "Endo_C4_arterial-CXCL12" , "Endo_C5_venous-SELE"     , "Endo_C6_venous-ACKR1"    ,
                                                      "Endo_C7_capillary-CA4"  ,  "Endo_C8_scavenging-C1QC" , "Endo_C9_lymphatics-PROX1",
                                                      "Endo_C10-MYL9"        ,    "Endo_C11-KRT8"        ,   
                                                      "Endo_C12-other" )),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c(    "Endo_C1_angio-INSR"      ,   "Endo_C2_angio-PGF"     ,   "Endo_C3_immature-APLNR",  
                                                         "Endo_C4_arterial-CXCL12" , "Endo_C5_venous-SELE"     , "Endo_C6_venous-ACKR1"    ,
                                                         "Endo_C7_capillary-CA4"  ,  "Endo_C8_scavenging-C1QC" , "Endo_C9_lymphatics-PROX1",
                                                         "Endo_C10-MYL9"        ,    "Endo_C11-KRT8"        ,   
                                                         "Endo_C12-other"  ))

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
                     comparisons = list( c("GC","GN"),c("GC","GCPM"),c("GN","GCPM")))+
  theme_classic()+
  scale_fill_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  facet_wrap(~Cluster,ncol=4,scales = "free")

print(box_p)





library(CellChat)


cellchat1<-readRDS("6.cell_chat/1.CAF/GCPM.RDS")
cellchat2<-readRDS("6.cell_chat/1.CAF/GC.RDS")

cellchat <- mergeCellChat(list(cellchat1, cellchat2), add.names = c("GCPM", "GC"))

groupSize <- as.numeric(table(cellchat1@idents))
levels(cellchat1@idents) 
vertex.receiver = seq(1,4) # a numeric vector


pathways.show <- "COMPLEMENT"
select_list<-"PTN|IL|CSF|TEK|VEGF"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])

cell_chat_p2<-netVisual_bubble(cellchat1, pairLR.use = pairLR.use,  sources.use = c(      "DC"    ,      "Endo" ,      
                                                                                                            "Epi"    ,     "iCAF"      ,  "MAIT"     ,   "Mast"  ,      "mCAF"   ,     "Mono" ,       "Neutro"  ,   
                                                                                                            "NK"     ,     "Plasma"    ,  "TAM"      ,   "Tgd"   ,      "vCAF" ), targets.use = 7, remove.isolate = FALSE)







##correlation
combined1<-readRDS("3.Cluster/12.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))

combined$Cluster<-factor(combined$Cluster,levels=c( "B"          ,  "CD4_conv"   ,  "CD4_Treg"  ,   "CD8"      ,    "Chief_cells",  "DC"      ,     "Endo"  ,      
                                                    "Epi"         , "iCAF"      ,   "MAIT"     ,    "Mast"       ,  "mCAF"    ,     "Mono" ,       
                                                    "Mucoid_cells", "Neutro"     ,  "NK"        ,      "Plasma"     ,  "TAM"     ,     "Tgd"  ,       
                                                    "vCAF"))

tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)




library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p


cell_pair_list<-list(c("vCAF","Endo"),
                     c("Endo","CD8"))

for (i in 1:length(cell_pair_list)){
  cell_pair<-cell_pair_list[[i]]
  
  immune_mt<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
  
  Cluster1<-cell_pair[[1]]
  Cluster2<-cell_pair[[2]]
  
  
  CD8_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster1,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster1,-1]))
  
  
  Neutro_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster2,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster2,-1]))
  
  
  
  
  merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")
  
  cor.value<-cor(merged_table$per.x,merged_table$per.y)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
  p.value<-signif(p.value,digits = 2)
  
  p <- ggplot(merged_table, aes(x=per.x, y=per.y))+
    geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=2,color="#374E55")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste(Cluster1,"and",Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 11)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p2<-p1+p2+plot_layout(ncol=1,nrow=2)







##correlation of subcluster
combined<-subset(combined1,idents=c("GC","GCPM","GN"))
combined$Cluster<-as.character(combined$Cluster)
combined$SubCluster<-as.character(combined$SubCluster)
combined$Cluster[which(combined$Cluster=="Endo")]<-combined$SubCluster[which(combined$Cluster=="Endo")]

tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
tmp_table<-tmp_table[which(!(tmp_table$sample_ID %in% c("GSM5573482_sample17"))),]
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)




library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p




cell_pair_list<-list(c("vCAF","Endo_C1_angio-INSR"),
                     c("vCAF","Endo_C2_angio-PGF" ),
                     c("vCAF","Endo_C3_immature-APLNR"),
                     c("vCAF","Endo_C4_arterial-CXCL12" ),
                     c("vCAF","Endo_C5_venous-SELE"      ),
                     c("vCAF" , "Endo_C6_venous-ACKR1" ))

for (i in 1:length(cell_pair_list)){
  cell_pair<-cell_pair_list[[i]]
  
  immune_mt<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
  
  Cluster1<-cell_pair[[1]]
  Cluster2<-cell_pair[[2]]
  
  
  CD8_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster1,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster1,-1]))
  
  
  Neutro_tab<-data.frame(
    sample_ID=names(immune_mt[immune_mt$Cluster==Cluster2,-1]),
    per=as.numeric(immune_mt[immune_mt$Cluster==Cluster2,-1]))
  
  
  
  
  merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")
  
  cor.value<-cor(merged_table$per.x,merged_table$per.y)
  cor.value<-round(cor.value,2)
  p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
  p.value<-signif(p.value,digits = 2)
  
  p <- ggplot(merged_table, aes(x=per.x, y=per.y))+
    geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=2,color="#374E55")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste(Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 11)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p3<-p1+p2+p3+p4+p5+p6+plot_layout(ncol=3)








total_p1<-ggarrange(umap_p1,stack_p,
                    ncol = 2,nrow = 1,
                    widths  = c(3,1))

total_p2<-ggarrange(total_p1,umap_p2,
                    ncol = 1,nrow = 2,
                    heights  = c(2,1))

total_p3<-ggarrange(total_p2,dot_p,
                    ncol = 2,nrow = 1,
                    widths = c(2,1))

total_p4<-ggarrange(pheatmp_p2,box_p,
                    ncol = 2,widths=c(1,1.8))

total_p5<-ggarrange(cor_p2,cell_chat_p2,cor_p3,
                    ncol = 3,widths=c(0.7,1.3,2))

total_p<-ggarrange(total_p3,total_p4,total_p5,
                   ncol = 1,nrow = 3,
                   heights = c(7,5,4.5))

#print(total_p3)


ggsave("14.Figure_plot/Figure6_S6_Endo.pdf",total_p,width = 12,height = 16.5)



