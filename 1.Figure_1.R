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





setwd("/home/zhengyq/data/single_cell/3.GCPM1")


sub_sce<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE","PBMC" ))
dir.create("3.Cluster/13.plot/1.total/")
sub_sce$Cluster<-as.character(sub_sce$Cluster)
sub_sce$Cluster[which(sub_sce$Cluster %in% c("iCAF" ,  "mCAF"))]<-"CAF"
sub_sce$Cluster[which(sub_sce$Cluster %in% c("CD4_conv",  "CD4_Treg"))]<-"CD4"

sub_sce$Cluster<-factor(sub_sce$Cluster,levels=c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,  "NKT"   ,    "MAIT"     , "Tgd"     ,
                                                     "Mast"   ,  "DC"  ,              
                                                     "Mono"       ,            
                                                     "Mph"   ,  "Neutro" ,          
                                                     "CAF", "PC" , "EC"    ))       


Cluster<-c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"  ,  "NKT"  ,     "MAIT"     , "Tgd"     ,
               "Mast"   ,  "DC"  ,              
               "Mono"       ,            
               "Mph"   ,  "Neutro" ,          
               "CAF", "PC" , "EC"    )
g.colSet <- readRDS("g.colSet.RDS")
g.colSet1<-c("#1F77B4" ,"#FF7F0E" , "#9467BD", "#AEC7E8","#E377C2", "#7F7F7F",
             "#BCBD22" ,"#17BECF" ,"#FFBB78", "#98DF8A" ,"#FF9896", "#C5B0D5" ,"#C49C94",
             "#F7B6D2" ,"#C7C7C7" ,"#DBDB8D", "#9EDAE5")
names(g.colSet1)<-Cluster
#g.colSet1<-list("Cluster"=g.colSet1)



#g.colSet1<-list("SubCluster"=g.colSet1)

library(ggforce)
sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
class_avg <- sub_sce@meta.data %>%
  group_by(Cluster) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )


p1<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Cluster))+
  geom_point(aes(color=Cluster),size=0.1) +
  #geom_mark_hull(aes(label=Cluster,fill=Cluster,group=Cluster),alpha=0.1,concavity = 1)+
  stat_ellipse(data=sub_sce@meta.data,aes(group=Cluster,color=Cluster),
               alpha=0.5,show.legend = F,
               level = 0.95, linetype = 'dashed')+
  scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                     labels= c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   , "NKT"   ,    "MAIT"     , "Tgd"     ,
                                   "Mast"   ,  "DC"  ,              
                                   "Mono"       ,            
                                   "Mph"   ,  "Neutro" ,          
                                   "CAF", "PC" , "EC"    ),
                     values = g.colSet1)+
  ggtitle("Clustering of TME cells")+
  geom_text_repel(aes(label = Cluster), data = class_avg,color="black")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(ncol=1,override.aes = list(size=4)))
p1


p2<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=Cluster))+
  geom_point(aes(color=Cluster),size=0.1,alpha=0.6) +
  scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                     labels= c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,   "NKT"   ,     "MAIT"     , "Tgd"     ,
                                   "Mast"   ,  "DC"  ,              
                                   "Mono"       ,            
                                   "Mph"   ,  "Neutro" ,          
                                   "CAF", "PC" , "EC"    ),
                     values = g.colSet1)+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        ##axis.title = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.text  = element_text(size = 12),
        strip.background = element_blank())+ 
  facet_wrap(~Type,ncol=5)
p2




library(ggdensity)
library(viridis)


Idents(sub_sce)<-"Type"
type_list<-c("GCPM","GN","GC","PBMC","PE")
for (i in 1:length(type_list)){
  sub_sce1<-subset(sub_sce,idents=type_list[i])
  sample_ID<-runif(20000,1,ncol(sub_sce1))
  cell_ID<-colnames(sub_sce1)[sample_ID]
  if (i==1){
    total_lst<-cell_ID
  }else{
    total_lst<-c(total_lst,cell_ID)
  }
}

sub_sce$merge_var<-"Drop"
sub_sce$merge_var[colnames(sub_sce) %in% total_lst]<-"Keep"
Idents(sub_sce)<-sub_sce$merge_var
combined<-subset(sub_sce,idents="Keep")

sample_ID<-colnames(combined)[combined$Type=="PE"&combined$Cluster=="Neutro"]
num1<-runif(2200,1,length(sample_ID))
cell_ID<-sample_ID[num1]

keep_id_list<-setdiff(colnames(combined),cell_ID)

combined$merge_var<-"Drop"
combined$merge_var[colnames(combined) %in% keep_id_list]<-"Keep"
Idents(combined)<-combined$merge_var
combined<-subset(combined,idents="Keep")


umap_data<-combined@meta.data
p3<-ggplot(data=umap_data,aes(UMAP_1,UMAP_2))+
  geom_hex(bins=60) + 
  scale_fill_gradientn(colors=c("#3300FF" ,"#001AFF" ,"#0066FF" ,"#00B3FF" ,"#00FFFF" ,
                                "#00FF66" ,"#33FF00", "#80FF00" ,"#CCFF00","#FFE500" ,"#FFE500" ,"#FFE500" ,"#FFE500" ,"#FF9900" ,"#FF9900" ,"#FF9900" ,"#FF9900" ,"#FF4D00","#FF4D00","#FF4D00", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000"))+
  #scale_fill_gradient2(low = "#4900FF", mid = "#00FF92", high = "#FF0000",
  #                       midpoint = 70)+
  #geom_point(size=0)+

  #scale_fill_distiller(palette= "Spectral", direction=1) +
  #scale_fill_gradient2(low = "blue",mid = "green",high = "red4",midpoint = 50)+
  #scale_fill_viridis(option = "A",begin = 0,end = 1 )+
  
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        ##axis.title = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey",size=0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 10),
        strip.text  = element_blank(),
        strip.background = element_blank())+ 
  facet_wrap(~Type,ncol=5)
p3




Cal_average_exp<-function(obj,group_mat,markers){
  DefaultAssay(obj)<-"RNA"
  obj<-ScaleData(obj)
  obj<-subset(obj,features = markers)
  Exp<-t(obj@assays$RNA@scale.data)
}


##T_NK
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Cluster
combined<-combined1

combined$Cluster<-as.character(combined$Cluster)
combined$Cluster[which(combined$Cluster %in% c("MSC","apCAF", "iCAF" ,  "PC",  "mCAF"))]<-"CAF"
combined$Cluster[which(combined$Cluster %in% c("CD4_conv",  "CD4_Treg"))]<-"CD4"


combined$Cluster<-factor(combined$Cluster,levels=c(    "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,"NKT"   ,       "MAIT"     , "Tgd"     ,
                                                       "Mast"   ,  "DC"  ,              
                                                       "Mono"       ,            
                                                       "Mph"   ,  "Neutro" ,          
                                                       "CAF", "PC" , "EC"    ))


markers<-c("CD79A","MS4A1",  "JCHAIN"  , "MZB1",
           "CD3D","CD4" ,     
           "CD8A","GZMK"  ,     
           "GNLY"  ,   "NKG7"  ,  
           "TPSAB1","CPA3",
           
           "LILRA4","IDO1" ,"CD1C", 
           "S100A9",  "S100A8",  "CD14", 
            "CD68",    "C1QC",
           "FCGR3B","CSF3R",
           "COL1A1","LUM","ACTA2","RGS5","PECAM1","FLT1"
           
)


if (T){
  sample_ID<-runif(20000,1,ncol(combined))
  combined$merge_var<-"Drop"
  combined$merge_var[sample_ID]<-"Keep"
  Idents(combined)<-combined$merge_var
  combined<-subset(combined,idents="Keep")
}


DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)



combined1<-subset(combined, features =markers )
Exp_table<-data.frame(Cluster=combined1$Cluster,t(combined1@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))


rowData1<-data.frame(Markers=c(
  rep("B cell markers",4),
  rep("CD4 markers",2),
  rep("CD8 markers",2),
  rep("NK markers",2),
  rep("Mast markers",2),
 
  rep("DC markers",3),
  rep("Mph/Mono markers",5),
  rep("Neutro markers",2),
  
  rep("Stromal markers",6)
))
rownames(rowData1)<-markers
rowData1$Markers<-factor(rowData1$Markers,levels=c(unique(rowData1$Markers)))
ave_Exp_table<-data.frame(ave_Exp_table)
ave_Exp_table<-ave_Exp_table[markers,]



ave_Exp_table<-t(ave_Exp_table)

pheat_p1<-pheatmap::pheatmap(ave_Exp_table,
                             color=colorRampPalette(c("navy" ,"navy" ,'white','#D97777','#7E2324'), bias=1)(50), border_color="grey",
                             #color=colorRampPalette(c('#3C5488FF','white','#DC0000FF'), bias=1)(50), border_color=NA,
                             
                             #color=colorRampPalette(c('red',"blue"), bias=1)(50), border_color=NA,
                             #annotation_row = rowData1,
                             annotation_col = rowData1,
                             cluster_rows = F,
                             cluster_cols = F,
                             scale = "column",
                             gaps_row=c(2,8,13),
                             gaps_col = c(4,6,8,10,12,15,20,22), angle_col = "45",
                             annotation_legend = F,
                             fontsize_row = 9,fontsize_col= 9
                             
)
pheat_p1<-plot_grid(pheat_p1$gtable)













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

combined<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")

combined$Cluster<-as.character(combined$Cluster)
combined$Cluster[which(combined$Cluster %in% c("MSC","apCAF", "iCAF" ,  "PC",  "mCAF"))]<-"CAF"
combined$Cluster[which(combined$Cluster %in% c("CD4_conv",  "CD4_Treg"))]<-"CD4"
combined$Cluster[which(combined$Cluster %in% c("CD8",  "MAIT",  "Tgd"))]<-"CD8"


combined$Type<-factor(combined$Type,levels=c("GN" ,  "GC", "GCPM", "PE" ,"PBMC"  ))
meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$Cluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Type,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[OR.dist.mtx>3]<-3
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 0
pheatmp_p2<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               legend_title ="O/e",
                               cluster_rows = T,
                               cluster_cols = F,
                               treeheight_row = 10,
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
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-combined1
combined$Cluster<-as.character(combined$Cluster)



##All immune cells
combined$Cluster<-factor(as.character(combined$Cluster))

##CD4
tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c("Plasma",   "Mast",    
                                                    "Mph"   ,     "EC")),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c("Plasma",   "Mast",    
                                                     "Mph"   ,     "EC"))

color_list<-c("#BC3C29",  "#0072B5",
              "#E18727",  "#20854E",
               "#7876B1")
#names(color_list)<-c( "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c( "GN","GC","GCPM","PE","PBMC"  )
sum_table$Type<-factor(sum_table$Type,levels=c("GN","GC","GCPM","PE","PBMC"  ))
select_table<-sum_table
box_p<-ggplot(select_table, aes(Type, per, color = Type))+
  geom_boxplot(width=0.8,position=position_dodge(0),outlier.shape=NA)+
  geom_jitter( size=0.8, alpha=0.9,width = 0.1)+
  #ggrepel::geom_label_repel(aes(label=sample_ID))+
  stat_compare_means(label = "p.format", method = "t.test",
                     comparisons = list( c("GC","GN"),c("GC","GCPM"),c("GN","GCPM")))+
  theme_classic()+
  scale_color_manual(values = color_list)+
  labs(x = 'Type of tissue', y = 'Relative Abundance(%)') +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5,face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))+
  facet_wrap(~Cluster,ncol=2,scales = "free")

print(box_p)






##GO_BP CD8
Enrich_type<-"KEGG"
Type="All"
Select_cluster_list<-c("GCPM","GC" )
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/2.scRNA_all/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
  enrich_table<-read.table(path,sep="\t",header = T,row.names = 1)
  enrich_table$Cluster<-cluster
  total_table<-rbind(total_table,enrich_table)
  if (nrow(enrich_table)>10){
    enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
    enrich_table<-head(enrich_table,10)
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
  labs(y="Pathways")+
  theme_bw()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) 

S1 




library(patchwork)
total_p1<-ggarrange(p1,pheat_p1,widths = c(1,1.2))
  

total_p2<-p2+p3+plot_layout(ncol=1,nrow=2)

total_p4<-ggarrange(pheatmp_p2,box_p,S1,ncol=3,widths = c(0.7,1.1,1.2))

total_p<-ggarrange(total_p1,total_p2,total_p4,
                   ncol=1,nrow=3,heights = c(4,5,4))


pdf("14.Figure_plot/Figure1.pdf",width = 12,height = 13 )
print(total_p)
dev.off()
