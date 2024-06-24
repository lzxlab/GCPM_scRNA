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

  
Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","SD", "PR"))



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


meta.tb <- combined@meta.data
cellInfo.tb <- meta.tb

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Efficiency_ICB,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 10
OR.dist.mtx[OR.dist.mtx>5]<-5

pheatmp_p1<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = T,
                               cluster_cols = F,
                               
                               main="Distribution of all cells by ICB"
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p1<-plot_grid(pheatmp_p1$gtable)




cellInfo.tb <- meta.tb[meta.tb$Cluster %in% c(  "CAF" ,"EC",       
                                               "PC" ),]

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Efficiency_ICB,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 10
OR.dist.mtx[OR.dist.mtx>3]<-3
OR.dist.mtx<-OR.dist.mtx[,c("PR","SD","PD")]
pheatmp_p2<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(option = "C", n=7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = T,
                               cluster_cols = F,
                               treeheight_row = 10,
                               main="CAF by ICB",
                               angle_col = 90,
                               fontsize_row = 8,
                               fontsize_col = 8,
                               legend = T
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p2<-plot_grid(pheatmp_p2$gtable)





cellInfo.tb <- meta.tb[meta.tb$Cluster %in% c(  "Mph","Mono","DC","Neutro" ),]

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SubCluster,
                             colname.patient = "Sample_ID",
                             loc = cellInfo.tb$Efficiency_ICB,
                             out.prefix="Rplot",
                             pdf.width=4,pdf.height=6,verbose=1)

OR.dist.mtx<-OR.all.list$OR.dist.mtx
OR.dist.mtx[sapply(OR.dist.mtx, is.infinite)] <- 10
OR.dist.mtx[OR.dist.mtx>3]<-3
OR.dist.mtx<-OR.dist.mtx[,c("PR","SD","PD")]

pheatmp_p3<-pheatmap::pheatmap(OR.dist.mtx,
                               color = viridis::viridis(option = "C", n=7),
                               show_rownames = T,
                               show_colnames = T,
                               legend_labels = "OR",
                               cluster_rows = T,
                               cluster_cols = F,
                               treeheight_row = 10,
                               angle_col = 90,
                               main="Myeloid cells by ICB",
                               fontsize_row = 8,
                               fontsize_col = 8,
                               legend = F
)
#pheatmp_p2$gtable<-pheatmp_p2$gtable+scale_color_viridis(discrete=F, option = "D", begin = 0, end=1, direction=1)
pheatmp_p3<-plot_grid(pheatmp_p3$gtable)




##Dimplot
sub_sce<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")


Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents="GCPM")

Idents(sub_sce)<-sub_sce$Efficiency_ICB
sub_sce<-subset(sub_sce,idents=c("PR","PD"))



#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(   "EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
                                                          "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"   ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
                   "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"   ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tsne")<-Embeddings(object = sub_sce, reduction = "tsne")


#Merged plot
Cluster<-unique(sort(sub_sce$SubCluster))
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster
g.colSet1<-list("SubCluster"=g.colSet1)


umap_p2 <- sscVis::ssc.plot.tsne(sce, columns = "SubCluster", splitBy = "Efficiency_ICB",
                                 reduced.name = "umap",
                                 
                                 colSet=g.colSet1,size=0.1,label=0,
                                 par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                                 #vector.friendly=T,
                                 #
                                 fun.extra=function(p){ p + guides(color=guide_legend(ncol=1, override.aes=list(size=3))) +
                                     ggplot2::facet_wrap(~splitBy, ncol = 1,scales = "free_x")+
                                     theme_classic()+
                                     theme(legend.position = "right",
                                           plot.title = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title.x = element_text(size = 10),
                                           axis.title.y = element_text(size = 10),
                                           legend.title = element_text(size = 10),
                                           legend.text = element_text(size = 8),
                                           strip.background = element_blank())},
                                 par.geom_point = list(scale=0.8),
                                 par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                                 base_aspect_ratio = 1.35,
                                 p.ncol = 1)
print(umap_p2)





##Expresion compare CAF
Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c(  "CAF"  ))



data<-combined@assays$integrated@scale.data
data<-as.matrix(data)

gene_name<- c("C3")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Efficiency_ICB),
                               Cluster=factor(as.character(combined$Cluster)),
                               gene=as.numeric(data[gene_name,]))

#select_fpkm_matrix<-merge(select_fpkm_matrix,data1,by="Cell_ID")



Feature_p1<-FeaturePlot(sub_sce,features = "C3",split.by = "Efficiency_ICB",pt.size = 0.2)+
  plot_layout(ncol=1,nrow=2) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))


gene_name<-"C3"
combined<-sub_sce
combined<-AddMetaData(combined,combined@reductions$umap@cell.embeddings,col.name = colnames(combined@reductions$umap@cell.embeddings))
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined$gene<-combined@assays$RNA@scale.data[gene_name,]
gene_exp <- combined@meta.data 

library(viridis)
C3_exp<-ggplot(data=gene_exp,aes(x=UMAP_1,y=UMAP_2,color=gene))+
  geom_point(size=0.2)+
  theme_classic()+
  labs(title="C3 in stromal cells")+
  scale_colour_continuous(type ="viridis")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey"),
        strip.background = element_blank())+
  facet_wrap(~Efficiency_ICB,ncol=1)






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
  
  
  DF$Gene<-""
  DF$Gene[rownames(DF) %in% gene_list]<-rownames(DF)[rownames(DF) %in% gene_list]
  
  print(c(head(DF,20),tail(DF,20)))
  p <- ggplot(DF, aes(x = logFC, y = negLogPval)) +
    geom_point(aes(color = SigGenes))+ 
    geom_text_repel(aes(label=DF$Gene),max.overlaps = 2000)+
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
sub_sce<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("CAF"))

Idents(combined)<-combined$Efficiency_ICB
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
C1<-"PR"
C2<-"PD"

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
markers1<-markers$gene[which(markers$cluster==C1)]
markers$avg_log2FC[which(markers$cluster==C2)]<-(-(markers$avg_log2FC[which(markers$cluster==C2)]))

gene_list<-c("RPS4Y1","GSTM1","DCN","SPP1","C3","KRT18","GGH","LUM","THBS2","IGF1")
vocano_p1<-vocano_plot(markers, Sample_1 = C1, Sample_2 = C2, lfc = 0, pval = pval,
                       gene_list=gene_list,title="CAF by ICB (GCPM)")+
  theme(legend.position = "none")+
  #ylim(c(0,40))+
  theme(axis.text = element_text(size = 10), 
        legend.position = "none",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        plot.subtitle = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())
vocano_p1




##KEGG CAF
Enrich_type<-"KEGG"
Type="All"
Select_cluster_list<-c("PD","PR")
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/7.CAF_ICB/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
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

pathways<-c("Ribosome"                              ,      "Complement and coagulation cascades"  ,      
            "Prostate cancer"                       ,      "Osteoclast differentiation"     ,            
            "Pathways in cancer"                    ,      "Focal adhesion"               ,              
            "Rheumatoid arthritis"                  ,      "MAPK signaling pathway"      ,               
            "mTOR signaling pathway"                 ,     "Bladder cancer"          ,                   
            "Pathogenic Escherichia coli infection" ,      "Protein processing in endoplasmic reticulum",
            "Proteasome"                      ,           
            "Regulation of actin cytoskeleton" ,          
            "Protein export"                    ,                        
            "Phagosome"                                    )
total_table1<-total_table[which(total_table$Cluster %in% Select_cluster_list&total_table$Description %in% pathways),]
total_table1$Description<-factor(total_table1$Description,levels=rev(unique(total_table1$Description)))
total_table1$Cluster<-factor(total_table1$Cluster,levels=Select_cluster_list)


S1<- ggplot(total_table1, aes(x= Cluster, y=Description, size=Count, color=p.adjust, group=Cluster)) + 
  geom_point(alpha = 0.8) + 
  theme_classic()+
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")+
  scale_size(range = c(2, 6))+
  ggtitle(paste(Enrich_type,"in CAF"))+
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) 

S1 




sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")

Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mph"))

Idents(combined)<-combined$Efficiency_ICB
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
C1<-"PR"
C2<-"PD"

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
markers1<-markers$gene[which(markers$cluster==C1)]
markers$avg_log2FC[which(markers$cluster==C2)]<-(-(markers$avg_log2FC[which(markers$cluster==C2)]))

gene_list<-c("RPS4Y1","RGS1","PLAUR","SPP1","CD9","GPNMB","HIF1A","C5AR1","C3AR1","GGH")
vocano_p2<-vocano_plot(markers, Sample_1 = C1, Sample_2 = C2, lfc = 0, pval = pval,
                       gene_list=gene_list,title="TAM by ICB (GCPM)")+
  theme(legend.position = "none")+
  #ylim(c(0,40))+
  theme(axis.text = element_text(size = 10), 
        legend.position = "none",
        axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        plot.subtitle = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())
vocano_p2




##KEGG CAF
Enrich_type<-"KEGG"
Type="All"
Select_cluster_list<-c("PD","PR")
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/6.Mph_ICB/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
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

pathways<-c( "Rheumatoid arthritis"          ,  "Pathways in cancer"  ,   "Osteoclast differentiation"   , 
             "Leishmaniasis"                ,       "Graft-versus-host disease"      ,    
             "Type I diabetes mellitus"     ,       "Allograft rejection"          ,      
             "Phagosome"                   ,        "Toll-like receptor signaling pathway"            ,      
             "Autoimmune thyroid disease"   ,       "Antigen processing and presentation",
             "Cardiac muscle contraction"   ,       "Proteasome"   ,                      
            "Fc gamma R-mediated phagocytosis")

total_table1<-total_table[which(total_table$Cluster %in% Select_cluster_list&total_table$Description %in% pathways),]
total_table1$Description<-factor(total_table1$Description,levels=rev(unique(total_table1$Description)))
total_table1$Cluster<-factor(total_table1$Cluster,levels=Select_cluster_list)


S2<- ggplot(total_table1, aes(x= Cluster, y=Description, size=Count, color=p.adjust, group=Cluster)) + 
  geom_point(alpha = 0.8) + 
  theme_classic()+
  scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")+
  scale_size(range = c(2, 6))+
  ggtitle(paste(Enrich_type,"in Mph"))+
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) 

S2 



##Myeloid
sub_sce<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")


Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents="GCPM")

Idents(sub_sce)<-sub_sce$Efficiency_ICB
sub_sce<-subset(sub_sce,idents=c("PR","PD"))


dir.create("3.Cluster/13.plot/3.myeloid/")
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

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
                   "Neutro_C1"  ) 



names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)


sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tsne")<-Embeddings(object = sub_sce, reduction = "tsne")

Cluster<-unique(sort(sce$SubCluster))


#Merged plot
library(RColorBrewer)
g.colSet1 <- c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
               brewer.pal(8, "Set1"),
               brewer.pal(8, "Dark2")[1],
               "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(g.colSet1)<-Cluster
g.colSet1<-list("SubCluster"=g.colSet1)


umap_p3 <- sscVis::ssc.plot.tsne(sce, columns = "SubCluster", splitBy = "Efficiency_ICB",
                                 reduced.name = "umap",
                                 
                                 colSet=g.colSet1,size=0.1,label=0,
                                 par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                                 #vector.friendly=T,
                                 #
                                 fun.extra=function(p){ p + guides(color=guide_legend(ncol=1, override.aes=list(size=3))) +
                                     ggplot2::facet_wrap(~splitBy, ncol = 1,scales = "free_x")+
                                     theme_classic()+
                                     theme(legend.position = "right",
                                           plot.title = element_blank(),
                                           axis.text = element_blank(),
                                           axis.title.x = element_text(size = 10),
                                           axis.title.y = element_text(size = 10),
                                           legend.title = element_text(size = 10),
                                           legend.text = element_text(size = 8),
                                           strip.background = element_blank())},
                                 par.geom_point = list(scale=0.8),
                                 par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                                 base_aspect_ratio = 1.35,
                                 p.ncol = 1)
print(umap_p3)






##Expresion compare CAF
Idents(sub_sce)<-sub_sce$Type
combined<-subset(sub_sce,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$SubCluster
combined<-subset(combined,idents=c( "Mph_C1-TIMP1"    ,
                                    "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                    "Mph_C5_TAM-FBP1"  ))

combined$SubCluster<-as.character(combined$SubCluster)
combined$SubCluster<-gsub("_TAM|_TRM","",combined$SubCluster)

data<-combined@assays$integrated@scale.data
data<-as.matrix(data)

gene_name<- c("SPP1")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Efficiency_ICB),
                               Cluster=factor(as.character(combined$Cluster)),
                               gene=as.numeric(data[gene_name,]))

#select_fpkm_matrix<-merge(select_fpkm_matrix,data1,by="Cell_ID")



Feature_p3<-FeaturePlot(sub_sce,features = "SPP1",split.by = "Efficiency_ICB",pt.size = 0.2)+
  plot_layout(ncol=1,nrow=2) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))


gene_name<-"SPP1"
combined<-sub_sce
combined<-AddMetaData(combined,combined@reductions$umap@cell.embeddings,col.name = colnames(combined@reductions$umap@cell.embeddings))
DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined$gene<-combined@assays$RNA@scale.data[gene_name,]
gene_exp <- combined@meta.data 

library(viridis)
SPP1_exp<-ggplot(data=gene_exp,aes(x=UMAP_1,y=UMAP_2,color=gene))+
  geom_point(size=0.2)+
  theme_classic()+
  labs(title="SPP1 in TAM")+
  scale_colour_continuous(type ="viridis")+ 
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey"),
        strip.background = element_blank())+
  facet_wrap(~Efficiency_ICB,ncol=1)






##Expresion compare CAF
combined<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")
Idents(combined)<-combined$Type
combined<-subset(combined,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("CAF"))

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Vln_p1<-VlnPlot(combined,features = c("C3","THBS2","RPS4Y1","GGH"),
                group.by = "Efficiency_ICB",
                cols = c("#E64B35","#4DBBD5"),
                pt.size = 0) +
  plot_layout(ncol=2) &
  stat_compare_means(label="p.format",label.y.npc = 0.9) &
  #scale_fill_manual(values = g.colSet1)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10,angle=0,vjust=0.5,hjust=0.5),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
Vln_p1



##Expresion compare myeloid
combined<-readRDS("3.Cluster/13.Annotation/3.myeloid_sub_sce_annotation.rds")
Idents(combined)<-combined$Type
combined<-subset(combined,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mph"))


DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Vln_p2<-VlnPlot(combined,features = c("SPP1","GPNMB","RPS4Y1","GGH"),
                group.by = "Efficiency_ICB",
                cols = c("#E64B35","#4DBBD5"),
                pt.size = 0) +
  plot_layout(ncol=2) &
  stat_compare_means(label="p.format",label.y.npc = 0.9) &
  #scale_fill_manual(values = g.colSet1)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10,angle=0,vjust=0.5,hjust=0.5),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
Vln_p2





library(CellChat)
library(patchwork)

cellchat1<-readRDS("6.cell_chat/2.CAF_Eff/PD.RDS")
cellchat2<-readRDS("6.cell_chat/2.CAF_Eff/PR.RDS")



cellchat <- mergeCellChat(list(cellchat1, cellchat2), add.names = c("PD", "PR"))



#cellchat@idents<-factor(cellchat@idents,levels=c( "B" ,"Plasma"   ,    "CD4_conv"   ,   "CD4_Treg",  "CD8", "MAIT"  ,   "NK"   ,    "DC"  , "Macro"  ,    "Mono"  ,  
#                                                  "Mast" ,    "Endo", "APC","pAD"  ,  "SMC"            )) 

groupSize <- as.numeric(table(cellchat1@idents))
levels(cellchat1@idents) 
vertex.receiver = seq(1,4) 




pathways.show <- "COMPLEMENT"
select_list<-"CXCL|CCL|MDK|PTN|IL|CSF|THBS|C3"
select_list<-"CXCL|CCL|C3"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])
cell_chat_p1<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,  comparison = c(1, 2) , sources.use = c( "B"       ,     "CD4",     "CD8"  ,        "DC"         ,  "EC"  ,      
                                                                                                                "MAIT" ,        "Mast"      ,   "CAF"  ,       "Mono"     ,   
                                                                                                           "Neutro"    ,   "NK"     ,  "NKT"     ,         "Mph"  ,           
                                                                                                           "PC" ), targets.use = "Mph" , remove.isolate = FALSE)



pathways.show <- "COMPLEMENT"
select_list<-"CXCL|CCL|MDK|PTN|IL|CSF|THBS|C3|SPP1"
select_list<-"CXCL|SPP1"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])
cell_chat_p2<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,  comparison = c(1, 2) , sources.use = "Mph", targets.use = c( "B"       ,     "CD4",     "CD8"  ,        "DC"         ,  "EC"  ,      
                                                                                                                                "MAIT" ,        "Mast"      ,   "CAF"  ,       "Mono"     ,   
                                                                                                                                "Neutro"    ,   "NK"     ,  "NKT"     ,         "Mph"  ,           
                                                                                                                                "PC" ), remove.isolate = FALSE)




library(monocle3)
cds<-readRDS("5.monocle2/9.TAM_B2/TAM_B2_monocle3.rds")

colSet<-c(brewer.pal(12, "Set3")[-c(2,3,9,12)],"#b3b3b3",
          brewer.pal(8, "Set1"),
          brewer.pal(8, "Dark2")[1],
          "#fc4e2a","#fb9a99","#f781bf","#e7298a")
names(colSet)<-c(    "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                     "DC_C4_cDC2-CD1C"   ,  
                     
                     "Mono_C1-S100A9"  ,
                     "Mono_C2-FCN1"   ,  "Mono_C3-CX3CR1"  , "Mono_C4-HSPA1A"      , 
                     "Mph_C1-TIMP1"    ,
                     "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                     "Mph_C5_TAM-TREM2",  "Mph_C6_TAM-VSIG4" ,
                     "Neutro_C1-CD16B" ,
                     "Neutro_C2-PROK2"     )

cds_subset <- cds[, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GCPM") &
                        Efficiency_ICB %in% c("PR","PD")
                    ) %>%
                    row.names
]


pseudotime_TAM<-plot_cells(cds_subset, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  ggtitle('Cell trajectory of TRM-derived TAM')+
  facet_wrap(~Efficiency_ICB,ncol=2,scales = "free_y")



Track_genes_sig <- c("SPP1")
cds_subset <- cds[Track_genes_sig, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GCPM") &
                        Efficiency_ICB %in% c("PR","PD")
                    ) %>%
                    row.names
]
#基因表达趋势图
rownames(cds_subset)<-"SPP1 (Branch 2)"
plot_genes10<-plot_genes_in_pseudotime(cds_subset, color_cells_by="SubCluster", 
                                       ncol = 1,cell_size = 0.3,label_by_short_name = F,min_expr = 0.01)+
  theme_classic()+
  scale_color_manual(values = colSet)+
  ggtitle("SPP1 (Branch 2)")+
  #scale_y_continuous(limits = c(0,3),breaks = seq(0,3),labels = c(0.1,1,10,100))+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.background = element_blank(),
        plot.title = element_text(size=10,hjust = 0.5))+
  facet_wrap(~Efficiency_ICB,ncol=1,scales = "free_x")







##pheudotime_CAF
cds<-readRDS("5.monocle2/5.CAF_monocle3/CAF_monocle3.rds")
cds_subset <- cds[, 
                  colData(cds) %>%
                    subset(
                      Type %in% c("GCPM") &
                        Efficiency_ICB %in% c("PR","PD")
                    ) %>%
                    row.names
]


library(monocle3)
#library(monocle)
pseudotime_CAF<-plot_cells(cds_subset, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                          label_leaves = FALSE,  label_branch_points = FALSE)+ 
  ggtitle('Cell trajectory of CAF')+
  theme_classic()+
  theme(legend.position = "right")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10),
        strip.background = element_blank())+
  facet_wrap(~Efficiency_ICB,ncol=2,scales = "free_y")



#挑选top10画图展示
Track_genes_sig <- c("ACTA2","FAP","CD34")
Cluster<-c( "Endo_C1-ACKR1" , "Endo_C2-PLVAP" , "Endo_C3-MAGI1",
            "Endo_C4-RGS5"  , "iCAF_C1-CXCL14", "mCAF_C1-THBS2" ,
            "mCAF_C2-TGFBI",  "mCAF_C3-SLIT2" ,
            "PC_C1-MYH11"  ,  "PC_C2-RGS5" )
colSet <- c(RColorBrewer::brewer.pal(8,"Set3"),
            RColorBrewer::brewer.pal(8,"Set2"))
names(colSet)<-Cluster


#基因表达趋势图
plot_genes1<-plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="SubCluster", 
                                      ncol = 1,cell_size = 0.3,min_expr = 0.1)+
  facet_wrap(~f_id,ncol=1,scales = "free_x")+
  theme(legend.position = "none")+
  theme_classic()+
  scale_color_manual(values = colSet)+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        strip.background = element_blank(),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size=12,hjust = 0.5))





##Expresion compare myeloid
combined<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined)<-combined$Type
combined<-subset(combined,idents="GCPM")
Idents(combined)<-combined$Efficiency_ICB
combined<-subset(combined,idents=c("PD","PR"))


DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)

Vln_p3<-VlnPlot(combined,features = c("SPP1","CD44"),
                split.by = "Efficiency_ICB",group.by = "Cluster",
                cols = c("#E64B35","#4DBBD5"),
                pt.size = 0) +
  plot_layout(ncol=2) &
  stat_compare_means(label="p.signif",label.y.npc = 0.9) &
  #scale_fill_manual(values = g.colSet1)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
Vln_p3



total_p1<-pheatmp_p3+umap_p3+pheatmp_p2+umap_p2+plot_layout(ncol=4,nrow = 1,widths = c(0.65,1.1,0.75,1.1))

total_p3<-cell_chat_p1+cell_chat_p2+plot_layout(ncol=2,nrow = 1,widths = c(1,1))

right_p<-vocano_p1+vocano_p2+plot_layout(nrow=2,ncol = 1)
midt_p<-S1+S2+plot_layout(nrow=2,ncol = 1)
left_r<-ggarrange(Vln_p1,Vln_p2,ncol = 1,nrow=2)
left_p<-ggarrange(C3_exp,SPP1_exp,ncol = 1,nrow=2,common.legend = T,legend = "right")
total_p5<-ggarrange(right_p,left_r,left_p,midt_p,ncol=4,widths = c(0.9,0.8,0.7,1.6))


total_p4<-ggarrange(total_p1,total_p5,total_p3,Vln_p3,ncol=1,nrow = 4,heights=c(3.5,6,3,2.5))


total_p<-total_p4
dir.create("14.Figure_plot/")
ggsave("14.Figure_plot/Figure6.pdf",total_p,width = 12,height = 15)
