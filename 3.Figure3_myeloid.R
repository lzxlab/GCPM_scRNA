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


sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE","PBMC" ))
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(  "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                                                         "DC_C4_cDC2-CD1C"   ,  
                                                         
                                                         
                                                         "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                                         "Mph_C1-TIMP1"    ,
                                                         "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                         "Mph_C5_TAM-FBP1",
                                                         "Neutro_C1-CD16B"  ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c("DC_C1"  ,   "DC_C2"  ,   "DC_C3" ,  
                  "DC_C4"   ,  
                  
                  "Mono_C1"  ,
                  "Mono_C2"   ,  "Mono_C3"  ,
                  "Mph_C1"    ,
                  "Mph_C2" ,"Mph_C3", "Mph_C4", 
                  "Mph_C5",  
                  "Neutro_C1" ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tsne")<-Embeddings(object = sub_sce, reduction = "tsne")

Cluster<-unique(sort(sce$SubCluster))
library(RColorBrewer)
g.colSet1 <- c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
                        brewer.pal(8, "Set1"),
                        brewer.pal(8, "Dark2")[1],
                        "#fc4e2a","#fb9a99","#f781bf","#e7298a")
g.colSet1<-g.colSet1[1:length(Cluster)]
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


p1<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c(   "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                                  "DC_C4_cDC2-CD1C"   ,  
                                  
                                  
                                  "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                  "Mph_C1-TIMP1"    ,
                                  "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                  "Mph_C5_TAM-FBP1",
                                  "Neutro_C1-CD16B"   ),
                     values = g.colSet1)+
  ggtitle("Clustering of myeloid cells")+
  geom_text(aes(label = SubCluster1), data = class_avg,color="black",size=3)+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(), 
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 9))+ 
  guides(colour = guide_legend(override.aes = list(size=3)))


p1




#Merged plot
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("GC" ,"GCPM","PE","PBMC" ))
sub_sce$Type<-factor(sub_sce$Type,levels=c("GC" ,"GCPM","PE","PBMC" ))
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(    "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                                                           "DC_C4_cDC2-CD1C"   ,  
                                                           
                                                           
                                                           "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                                           "Mph_C1-TIMP1"    ,
                                                           "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                           "Mph_C5_TAM-FBP1",
                                                           "Neutro_C1-CD16B"  ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "DC_C1"  ,   "DC_C2"  ,   "DC_C3" ,  
                   "DC_C4"   ,  
                   
                   "Mono_C1"  ,
                   "Mono_C2"   ,  "Mono_C3"  , 
                   "Mph_C1"    ,
                   "Mph_C2" ,"Mph_C3", "Mph_C4", 
                   "Mph_C5",  
                   "Neutro_C1"   ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tsne")<-Embeddings(object = sub_sce, reduction = "tsne")

g.colSet1<-c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
             brewer.pal(8, "Set1"),
             brewer.pal(8, "Dark2")[1],
             "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(g.colSet1)<-Cluster
g.colSet1<-list("SubCluster"=g.colSet1)


p2 <- sscVis::ssc.plot.tsne(sce, columns = "SubCluster", splitBy = "Type",
                            reduced.name = "umap",
                            
                            colSet=g.colSet1,size=0.1,label=0,
                            par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                            #vector.friendly=T,
                            #
                            fun.extra=function(p){ p + guides(color=guide_legend(ncol=1, override.aes=list(size=4))) +
                                ggplot2::facet_wrap(~splitBy, ncol = 2)+
                                theme_classic()+
                                theme(legend.position = "none",
                                      plot.title = element_blank())},
                            par.geom_point = list(scale=0.8),
                            par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                            base_aspect_ratio = 1.35,
                            p.ncol = 2)
print(p2)






##MIlo
library(miloR)
milo.obj<-readRDS("4.characteristics/7.4.Milo_myeloid.milo.obj.RDS")
milo.res<-readRDS("4.characteristics/7.4.Milo_myeloid.milo.res.RDS")
nh_graph_myeloid <- plotNhoodGraphDA(milo.obj, milo.res, layout="umap",alpha=0.1) +
  theme_classic()+
  labs(title = "Nhoods of myeloid cells",x="UMAP_1",y="UMAP_2")+
  scale_fill_gradientn(colors = c('navy',"white","red3"))+
  theme(axis.title = element_text(size = 10),
        axis.text = element_blank(), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 12,vjust = 0.5,hjust=0.5))

milo.res$SubCluster<-as.character(milo.res$SubCluster)

order_tab<-dplyr::summarize(group_by(milo.res,SubCluster),FC=mean(logFC))
order_tab<-order_tab[order(order_tab$FC),]
milo.res$SubCluster<-factor(milo.res$SubCluster,levels=order_tab$SubCluster)

beeswarm_myeloid<-plotDAbeeswarm(milo.res, group.by = "SubCluster")+
  theme_classic()+
  scale_color_gradientn(colors = c('navy',"white","red3"))+
  geom_hline(yintercept = c(0),linetype=c("solid"),size=0.5,color="grey")+
  theme(axis.text.y = element_text(size = 9), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5))
nh_graph_myeloid+beeswarm_myeloid





sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(    "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                                                           "DC_C4_cDC2-CD1C"   ,  
                                                           
                                                           
                                                           "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                                           "Mph_C1-TIMP1"    ,
                                                           "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                           "Mph_C5_TAM-FBP1",
                                                           "Neutro_C1-CD16B" ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "DC_C1"  ,   "DC_C2"  ,   "DC_C3" ,  
                   "DC_C4"   ,  
                   
                   "Mono_C1"  ,
                   "Mono_C2"   ,  "Mono_C3"  , 
                   "Mph_C1"    ,
                   "Mph_C2" ,"Mph_C3", "Mph_C4", 
                   "Mph_C5",  
                   "Neutro_C1") 


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster

features=list("M1"=c("TNF","IRF5","CD86","IL1B","KYNU","IRF1","CCR7","CD40","CXCL9","IL23A"),
              "M2"=c("CLEC7A","CCL4","TGFB1","CTSD","CTSC","CTSB","CTSA","FN1","MSR1","TNFSF12","LYVE1","IL4R","MMP19","CCL20","MMP14","CCL18","CD276","VEGFB","CCL22","IRF4","CCL17"),
              "Angiogenesis"=c("TYMP","VCAN","CD44","FYN","VEGFA","TNFAIP6","E2F3","MMP9","ITGAV","SPP1","CXCR4","PTK2","CCND2","EZH2"),
              "Phagocytosis"=c("CD163","C1QB","MERTK","MRC1"),
              "Checkpoint"=c("LAIR1","SIRPA","HAVCR2","NECTIN2","PVR","CD274","ADORA2A","IDO1","BTLA")
              )

dot_p<-DotPlot(sub_sce, features = features,dot.scale = 4,
               #cols = RColorBrewer::brewer.pal(3,"BuPu")
                 )+
  theme_bw()+
  RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_distiller(palette = "BuPu",direction = 1)+
  #scale_color_brewer()+
  #scale_color_viridis(discrete=F, option = "E", begin = 0, end=1, direction=1)+
  theme(legend.position = "right")+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 90),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))
dot_p




##All Cells
combined<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")


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
g.colSet1 <- c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
               brewer.pal(8, "Set1"),
               brewer.pal(8, "Dark2")[1],
               "#fc4e2a","#fb9a99","#f781bf","#e7298a")
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

sample_list<-unique(sort(combined$Sample_ID))
sample_list<-sample_list[which(!(sample_list %in% c("GSM5573482_sample17")))]
Idents(combined)<-combined$Sample_ID
combined<-subset(combined,idents=sample_list)

meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb[which(meta.tb$Cluster %in% c("Mph","Mono","Neutro","DC")),]

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
OR.dist.mtx<-OR.dist.mtx[,c("GN","GC","GCPM","PE","PBMC")]
pheatmp_p2<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = T,
                               cluster_cols = F,
                               treeheight_row = 15,
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
setwd("/home/zhengyq/data/single_cell/3.GCPM1/")
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))

##All immune cells
combined$Cluster<-factor(as.character(combined$Cluster))

##CD4
tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

cluster_list1<-c(   "DC_C1_pDC-GPR183"  , "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9",
                    "DC_C2_cDC1-LAMP3"   ,  "Mph_C4_TAM-SPP1", 
                    "Neutro_C1-CD16B")
cluster_list2<-c(   "GPR183+ GPR183"  , "F13A1+ TRM" ,"CXCL9+ Mph",
                    "LAMP3+ cDC1"   ,  "SPP1+ TAM", 
                    "CD16B+ Neutro")
for (i in 1:length(cluster_list1)){
  sum_table$Cluster[which(sum_table$Cluster==cluster_list1[[i]])]<-cluster_list2[[i]]
  
}
sum_table<-sum_table[which(sum_table$Cluster %in% c(   "F13A1+ TRM" ,"CXCL9+ Mph",
                                                      "LAMP3+ cDC1"   ,  "SPP1+ TAM")),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c( "F13A1+ TRM" ,"CXCL9+ Mph",
                                                      "LAMP3+ cDC1"   ,  "SPP1+ TAM"))

color_list<-c("#7876B1","#4DBBD5","#E64B35")
#names(color_list)<-c( "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c( "GN","GC","GCPM"  )
sum_table$Type<-factor(sum_table$Type,levels=c("GN","GC","GCPM"  ))
select_table<-sum_table
box_p<-ggplot(select_table, aes(Type, per, color = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
  geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format", method = "wilcox.test",
                     comparisons = list( c("GC","GN"),c("GC","GCPM"),c("GN","GCPM")))+
  theme_classic()+
  scale_color_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Percentage (%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  facet_wrap(~Cluster,ncol=2,scales = "free")

print(box_p)




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
  facet_wrap(~Type,ncol=1,nrow = 3)+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  theme(legend.position = "none")

cds_p3<-plot_cells(cds_subset, color_cells_by = "SubCluster", label_cell_groups = FALSE, 
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
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  theme(legend.position = "right")







##CD8/NK
library(patchwork)
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents=c("GC","GCPM"))



combined<-AddMetaData(combined,combined@reductions$umap@cell.embeddings,col.name = colnames(combined@reductions$umap@cell.embeddings))
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

gene_name<-"SPP1"
combined$gene<-combined@assays$RNA@scale.data[gene_name,]
gene_exp <- combined@meta.data 

library(viridis)
TAM_SPP1<-ggplot(data=gene_exp,aes(x=UMAP_1,y=UMAP_2,color=gene))+
  geom_point(size=0.2)+
  theme_classic()+
  labs(title="SPP1 in myeloid cells")+
  scale_colour_continuous(type ="viridis")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))+
  facet_wrap(~Type,ncol=2)



##Expresion compare CAF
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents=c("GCPM","GC"))



Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mph"))

library(limma)
combined$SubCluster<-as.character(combined$SubCluster)
combined$SubCluster<-strsplit2(combined$SubCluster,"-")[,1]
combined$SubCluster<-paste("Mph",strsplit2(combined$SubCluster,"_")[,2],sep="_")

color_list<-c("#4DBBD5","#E64B35")
#names(color_list)<-c( "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c( "GC","GCPM"  )

box_p1<-VlnPlot(combined,features = c("SPP1","CD81","TGFB1","MSR1","TREM2","HAVCR2"),split.by = "Type",
                group.by = "SubCluster",pt.size = 0,combine = F)
ps_list<-list()
for(i in 1:length(box_p1)){
  tmp_p<-box_p1[[i]]
  if (i==1){
    tmp_p<-tmp_p+scale_fill_manual(values = color_list)+
      stat_compare_means(label = "p.signif",label.y.npc = 0.9)+
      #ylim(c(-2,12))+
      theme(plot.title = element_text(size = 10,hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "none",
            legend.text = element_text(size = 10),
            strip.background = element_blank())
  }else{
    tmp_p<-tmp_p+scale_fill_manual(values = color_list)+
      stat_compare_means(label = "p.signif",label.y.npc = 0.9)+
      #ylim(c(-2,12))+
      theme(plot.title = element_text(size = 10,hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 10),
            legend.position = "none",
            legend.text = element_text(size = 10),
            strip.background = element_blank())
  }
  
  if (i<4){
    tmp_p<-tmp_p+theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank() )
  }
  
  
  ps_list[[length(ps_list)+1]]<-tmp_p
}

ps_list[[3]]<-ps_list[[3]]+theme(legend.position = "right")
box_pp<-ps_list[[1]]+ps_list[[2]]+ps_list[[3]]+
  ps_list[[4]]+ps_list[[5]]+ps_list[[6]]+plot_layout(ncol=3,nrow=2)


box_pp


test<-readRDS("5.monocle2/3.Mono/monocle2.RDS")
detach("package:monocle3", unload = TRUE)
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
  ggtitle('Cell trajectory of TRM-derived Mph')
monocle2_time

test$Type<-factor(test$Type,levels = c("PE","GCPM","GN","GC"))
monocle2_cluster<-plot_cell_trajectory(test,color_by = "Cluster",
                                       cell_size =0.4)+
  theme(legend.position = "none")+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory of TRM-derived Mph')+
  guides(colour = guide_legend(ncol=1,override.aes = list(size=3)))
monocle2_cluster

monocle2_split<-plot_cell_trajectory(test,color_by = "Cluster",show_branch_points = F,
                                     cell_size =0.4)+
  facet_wrap(~Type,ncol=2)+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(legend.position = "bottom")+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by tissues')+
  guides(colour = guide_legend(ncol=3,override.aes = list(size=3)))
monocle2_split

test_genes=c("CXCL9","SPP1")
plot_genes<-plot_genes_branched_pseudotime(test[test_genes,],
                               branch_point = 2,
                               color_by = "Cluster",
                               cell_size=0.2,
                               ncol = 2)+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(legend.position = "none")+
  theme(axis.text = element_blank(), axis.title.x = element_text(size = 10),axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  #ggtitle('Expression by pseudotime')+
  guides(colour = guide_legend(ncol=1,override.aes = list(size=3)))


##TF activity
cg_n<-readRDS("7.scenic/TAM/cg_n.RDS")
seu<-readRDS("7.scenic/TAM/seu_obj.RDS")
monocle1<-readRDS("5.monocle2/3.Mono/monocle2.RDS")

pseudotime_tab<-data.frame(Cell_ID=colnames(monocle1),
                           Pseudotime=monocle1$Pseudotime,
                           Cluster=monocle1$Cluster)
rownames(pseudotime_tab)<-pseudotime_tab$Cell_ID
cg_n<-AddMetaData(cg_n,metadata = pseudotime_tab)
seu<-AddMetaData(seu,metadata = pseudotime_tab)

cg_n$Cluster<-factor(cg_n$Cluster,levels=c("F13A1+ TRM","CXCL9+ Mph","TAM"))
Idents(cg_n)<-cg_n$Cluster
gene_list<-c( "MAF"  ,"FOXO3"  ,"TCF4" ,"KDM5A" ,"ELF1"  , "IRF8", "NFATC1"  , "HIC1" ,"SNAI1" , "RUNX3" ,
              "MLX" ,  "CEBPA", "ETV5" ,  "SREBF1" , "SPI1"   )

cg_n$Cluster1<-paste0("G",as.numeric(cg_n$Cluster))
cg_n$Cluster1<-factor(cg_n$Cluster1,levels = rev(unique(sort(cg_n$Cluster1))))
TF_dot<-DotPlot(cg_n, features = unique(gene_list),group.by = "Cluster1",
        cols = c("grey","blue"))+RotatedAxis()+
  ggtitle("Inferred TF activity")+
  scale_x_discrete("")+scale_y_discrete("")+
  theme_bw()+
  RotatedAxis()+
  #xlim(c("G1","G2","G3"))+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_distiller(palette = "BuPu",direction = 1)+
  #scale_color_brewer()+
  #scale_color_viridis(discrete=F, option = "E", begin = 0, end=1, direction=1)+
  theme(legend.position = "right")+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10,colour = RColorBrewer::brewer.pal(3,"Set2")),
        legend.title = element_blank(),
        
        legend.text = element_text(size = 10))
TF_dot

gene_list<-c("MAF" , "TCF4","IRF8","HIC1","MLX" ,"CEBPA")
select_TF_mt<-as.matrix(cg_n@assays$RNA[gene_list,])
select_TF_mt<-reshape2::melt(select_TF_mt)
colnames(select_TF_mt)<-c("TF_names","Cell_ID","value")

merged_tab<-merge(pseudotime_tab,select_TF_mt,by="Cell_ID")
TF_time<-ggplot(merged_tab,aes(x=Pseudotime,y=value,color=Cluster))+
  geom_smooth(se = F)+
  ggtitle("Inferred TF activity")+
  theme_classic2()+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA,color="grey"))+
  facet_wrap(~TF_names,scales="free_y",ncol=1)
TF_time


DefaultAssay(seu)<-"RNA"
seu<-ScaleData(seu)

select_TF_mt<-as.matrix(seu@assays$RNA@scale.data[gene_list,])
select_TF_mt<-reshape2::melt(select_TF_mt)
colnames(select_TF_mt)<-c("TF_names","Cell_ID","value")

merged_tab<-merge(pseudotime_tab,select_TF_mt,by="Cell_ID")
EXP_time<-ggplot(merged_tab,aes(x=Pseudotime,y=value,color=Cluster))+
  geom_smooth(se = F)+
  ggtitle("Expression of TF")+
  theme_classic2()+
  scale_color_manual(breaks = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     labels = c("F13A1+ TRM","CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA,color="grey"))+
  facet_wrap(~TF_names,scales="free_y",ncol=1)
EXP_time

TF_plot<-ggarrange(TF_dot,TF_time+EXP_time,ncol=1,nrow=2,heights = c(1,3))

total_p1<-p1+nh_graph_myeloid +beeswarm_myeloid+plot_layout(ncol=3,widths = c(2,2,1))
  

total_p2<-ggarrange(p2,pheatmp_p2,
                    ncol = 2,nrow = 1,
                    heights  = c(1,1))

total_p3<-ggarrange(box_p,box_pp,
                    ncol = 2,nrow = 1,
                    widths = c(1,2.3))


total_p4<-ggarrange(dot_p,total_p3,
                    ncol = 1)



SPP1_plot<-ggarrange(box_p1,plot_genes10,
                       ncol=1,nrow =2,heights = c(0.9,1))

total_p5<-ggarrange(cds_p1,cds_p3,plot_genes,box_pp,ncol = 4,
                    widths = c(0.9,0.5,1.3,0.7))


mid_p<-ggarrange(plot_genes+monocle2_split+plot_layout(ncol=1,heights = c(1,2.5)),TF_dot,ncol=1,nrow=2,heights = c(2,1))
right_p<-TF_time+EXP_time+plot_layout(ncol=2)
total_p5<-ggarrange(cds_p1+monocle2_time+monocle2_cluster+plot_layout(ncol=1,nrow=3,heights = c(2,1,1)),mid_p,right_p,ncol = 3,
                    widths = c(2.5,3,3))

total_p<-ggarrange(total_p1,total_p4,total_p5,
                   ncol = 1,nrow = 3,
                   heights = c(3.5,7,6))



ggsave("14.Figure_plot/Figure3.pdf",total_p,width = 12,height = 16.5 )



