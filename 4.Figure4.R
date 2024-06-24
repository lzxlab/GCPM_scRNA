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
theme_set(theme_minimal())
setwd("/home/zhengyq/data/single_cell/3.GCPM1/")



sub_sce<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC" ,"GCPM","PE"))
#sub_sce$SubCluster<-as.character(sub_sce$SubCluster)
#sub_sce$SubCluster<-gsub("DC_|Neutro_|TAM_|Mono_","",sub_sce$SubCluster)

sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c(  "EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
                                                         "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"    ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c( "EC_C1" , "EC_C2" , "EC_C3",  "EC_C4"  ,"EC_C5",
                   "iCAF_C1",
                   "mCAF_C1" , "mCAF_C2" ,
                   "PC_C1"  ,  "PC_C2" ) 



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
                     labels= c("EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
                               "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"    ),
                     values = g.colSet1)+
  ggtitle("Clustering of stromal cells")+
  geom_text(aes(label = SubCluster1), data = class_avg,color="black",size=3)+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))+ 
  guides(colour = guide_legend(override.aes = list(size=3)))


umap_p1


Idents(sub_sce)<-sub_sce$Type
sub_sce<-subset(sub_sce,idents=c("GCPM","PE","GC","GN"))
sub_sce$Type<-factor(sub_sce$Type,levels=c("GN","GC","PE","GCPM"))

umap_p2<-ggplot(sub_sce@meta.data ,aes(x=UMAP_1,y=UMAP_2,color=SubCluster))+
  geom_point(aes(color=SubCluster),size=0.1) +
  
  scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                     labels= c( "EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
                                "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"   ),
                     values = g.colSet1)+
  ggtitle("Clustering by tissue types")+
  theme_classic()+
  #theme(text=element_text(family="Arial",size=10)) +
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        legend.position="none")+ 
  facet_wrap(~Type,ncol=2,scales = "free")+
  guides(colour = guide_legend(override.aes = list(size=3)))


umap_p2






#MIlo
library(miloR)
milo.obj<-readRDS("4.characteristics/7.5.Milo_stromal.milo.obj.RDS")
milo.res<-readRDS("4.characteristics/7.5.Milo_stromal.milo.res.RDS")
nh_graph_stromal <- plotNhoodGraphDA(milo.obj, milo.res, layout="umap",alpha=0.1) +
  scale_fill_gradientn(colors = c('navy',"white","red3"))+
  theme_classic()+
  labs(title = "Nhoods of stromal cells",x="UMAP_1",y="UMAP_2")+
  theme(axis.title = element_text(size = 10),
        axis.text = element_blank(), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 10,hjust=0.5))

milo.res$SubCluster<-as.character(milo.res$SubCluster)

order_tab<-dplyr::summarize(group_by(milo.res,SubCluster),FC=mean(logFC))
order_tab<-order_tab[order(order_tab$FC),]
milo.res$SubCluster<-factor(milo.res$SubCluster,levels=order_tab$SubCluster)

beeswarm_stromal<-plotDAbeeswarm(milo.res, group.by = "SubCluster")+
  theme_classic()+
  scale_color_gradientn(colors = c('navy',"white","red3"))+
  geom_hline(yintercept = c(0),linetype=c("solid"),size=0.5,color="grey")+
  labs(title = "Differential Nhoods")+
  annotate("text",y = 1,x=10.5, label= "GC      GCPM")+
  theme(axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 10,hjust=0.5),
        plot.subtitle = element_text(size = 10,hjust=0.5))




sub_sce<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")
DefaultAssay(sub_sce)<-"RNA"
sub_sce<-ScaleData(sub_sce)
sub_sce$SubCluster<-factor(sub_sce$SubCluster,levels=c( "EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
                                                        "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"    ))       

Idents(sub_sce)<-sub_sce$SubCluster
new.cluster.id<-c("EC_C1" , "EC_C2" , "EC_C3",  "EC_C4"  ,"EC_C5",
                  "iCAF_C1",
                  "mCAF_C1" , "mCAF_C2" , 
                  "PC_C1"  ,  "PC_C2"  ) 


library(viridis)
names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)
sub_sce$SubCluster1<-Idents(sub_sce)

Idents(sub_sce)<-sub_sce$SubCluster1
dot_p<-DotPlot(sub_sce, features = c("CD320","PECAM1","ACKR1" ,"PLVAP","INSR","FLT1","MAGI1","CXCL14","PDGFRA"  , "CCL11", "CXCL1" ,"NFATC2", "NEAT1" , 
                                     
                                     "THBS2","IGF1" , "FAP",
                                     "COL1A1" ,"KRT8"   ,  "SLPI"  , "KRT19" , "SLIT2",  "ROR1"  , "KIF26B" ,"TUBA1B",
                                     "STMN1", "MYH11" , "ADIRF" , "SORBS2","CD36"  , "RGS5" ,  "THY1"),
               cols = c("grey","red2"),dot.scale=5)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  theme(legend.position = "none")+
  scale_color_viridis(discrete=F, option = "E", begin = 0, end=1, direction=1)+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10))
dot_p




##All Cells
combined<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")


tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,Type),num=sum(num))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="Type")
sum_table$per<-sum_table$num/sum_table$total_num*100
#sum_table$per[which(sum_table$Cluster %in% c("Epi",
#                                             "CAF","EC","SMC"))]<- (-(sum_table$per[which(sum_table$Cluster %in% c("Epi",
#                                                                                                                            "CAF","EC","SMC"))]))


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
combined$SubCluster<-factor(as.character(combined$SubCluster))

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

sum_table<-sum_table[which(sum_table$Cluster %in% c("iCAF_C1-CXCL14",
                                                    "mCAF_C1-THBS2" , "mCAF_C2-KRT8" ,
                                                    "PC_C1-MYH11"   )),]

sum_table$Cluster<-factor(sum_table$Cluster,levels=c("iCAF_C1-CXCL14",
                                                     "mCAF_C1-THBS2" , "mCAF_C2-KRT8" , 
                                                     "PC_C1-MYH11"   ))

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
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10), legend.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5,face="bold"),
        strip.background = element_blank())+
  facet_wrap(~Cluster,ncol=2,scales = "free")

print(box_p)






##Expresion compare CAF
combined<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")


Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("CAF"))
library(limma)
combined$SubCluster<-as.character(combined$SubCluster)
combined$SubCluster<-strsplit2(combined$SubCluster,"-")[,1]

combined$SubCluster<-factor(combined$SubCluster,levels=c( "iCAF_C1",
                                                          "mCAF_C1" , "mCAF_C2" ) )


Cluster<-c( "EC_C1" , "EC_C2" , "EC_C3",  "EC_C4"  ,
            "iCAF_C1",
            "mCAF_C1" , "mCAF_C2" , 
            "PC_C1"  ,  "PC_C2" ) 
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2"))
names(g.colSet1)<-Cluster

box_p1<-VlnPlot(combined,features = c("CXCL14","CCL11","ACTA2","COL1A1","FAP","THBS2","IGF1"),
                group.by = "SubCluster",
                cols = g.colSet1,
                pt.size = 0) +
  plot_layout(ncol=7) &
  #scale_fill_manual(values = g.colSet1)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p1







library(ggpubr)
library(magrittr)
library(ggsignif)
library(ggplot2)
library(ggsci)
library(RColorBrewer)




Gene_name="THBS2_CAF"
select_gene=Gene_name
outdir="3.Cluster/13.Plot/"
setwd("/home/zhengyq/data/single_cell/3.GCPM1/")

#Gene_name="SLFN5"
#outdir="pancancer_exp/"


RNA_matrix<-data.table::fread("~/data/TCGA/pancancer_expression/GDC-PANCAN.htseq_fpkm-uq.tsv",header = T,stringsAsFactors = F)
RNA_matrix<-data.frame(RNA_matrix)
gene_probe<-data.table::fread("~/data/TCGA/pancancer_expression/gencode.v22.annotation.gene.probeMap",header = T,stringsAsFactors = F)
phenotype_matrix<-data.table::fread("~/data/TCGA/pancancer_expression/GDC-PANCAN.basic_phenotype.tsv",header = T,stringsAsFactors = F)
colnames(gene_probe)[1]<-"Ensembl_ID"
colnames(RNA_matrix)[1]<-"Ensembl_ID"
gene_probe<-gene_probe[,c(1,2)]
RNA_matrix<-merge(gene_probe,RNA_matrix,by="Ensembl_ID")
RNA_matrix<-data.frame(RNA_matrix)
RNA_matrix<-RNA_matrix[which(!(duplicated(RNA_matrix$gene))),]
rownames(RNA_matrix)<-RNA_matrix$gene

RNA_matrix<-RNA_matrix[,3:ncol(RNA_matrix)]
colnames(RNA_matrix)<-gsub("\\.","-",colnames(RNA_matrix))

Gene_name1="THBS2_CAF"
Gene_name2="SPP1_TAM"

mymatrix<-as.matrix(RNA_matrix)
Tip_Markers<-c("SFRP4","COL8A1","COL1A1","MGP","COL1A2","COL3A1","CCDC80","CTHRC1","THBS2","IGF1","VCAN","C7","C3","FBLN1","CYP1B1","LXN","FNDC1","HOPX","DCN")
Endo_Markers<-c("SPP1","APOC1","APOE","GPNMB","LGMN","CTSD","TREM2","FABP5","CTSB","CD9","CTSL","LIPA","MSR1","ACP5","FTL","NUPR1","CD81","C1QC","CD68","CTSZ")
mysymbol1<-data.frame(Gene_set="Tip_ECs",Gene_symbol=Tip_Markers)
mysymbol2<-data.frame(Gene_set="Endo",Gene_symbol=Endo_Markers)
mysymbol<-rbind(mysymbol1,mysymbol2)

colnames(mysymbol)<-c("Gene_set","Gene_symbol")
head(mysymbol)
table(mysymbol$Gene_set)



type <- unique(mysymbol$Gene_set)
type
gs <- list()
for (i in type){
  tmp <- mysymbol$Gene_symbol[which(mysymbol$Gene_set == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs

library(GSVA)
es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=20)




##survplot
mt<-es.dif[1,]/es.dif[2,]
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=es.dif[1,],
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)


select_fpkm_matrix<-merge_matrix
colnames(select_fpkm_matrix)[1]<-"Tumor_ID"


library(limma)
select_fpkm_matrix$Tumor_ID<-gsub("\\.","-",select_fpkm_matrix$Tumor_ID)
Type_list<-strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,4]
select_fpkm_matrix<-select_fpkm_matrix[which(!(grepl("^1",Type_list))),]

select_fpkm_matrix$Tumor_ID<-paste(strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,1],strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,2],strsplit2(select_fpkm_matrix$Tumor_ID,"-")[,3],sep="-")



##survival_analysis
surv_table<-read.table("~/data/TCGA/survival/TCGA_survival.txt",header = T,stringsAsFactors = F,sep="\t")
surv_table<-na.omit(surv_table)

library(survival)
library(survminer)
library(patchwork)
##Cut off by OS
total_matrix<-merge(select_fpkm_matrix,surv_table,by="Tumor_ID")



Study_name<-"TCGA-STAD"
merged_matrix<-total_matrix[total_matrix$project_id==Study_name,]
res.cut <- surv_cutpoint(merged_matrix, #数据集
                         time = "PFI.Time", #生存状态
                         event = "PFI", #生存时间
                         variables = c("gene1") #需要计算的数据列名
)
merged_matrix$gene_level<-"Low"
merged_matrix$gene_level[merged_matrix$gene1>res.cut$cutpoint$cutpoint]<-"High"

surv_fit<-survfit(Surv(OS.Time , OS) ~ gene_level,data= merged_matrix)

gg_surv1<-ggsurvplot(surv_fit,
                     conf.int = F,
                     #fun = "cumhaz",
                     linetype =  1, # Change line type by groups
                     size=0.5,
                     censor = F,
                     #surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_classic(),# Change ggplot2 theme
                     
                     palette = c("#E72C19", "#224FA2"),
                     title = paste(Study_name,"(OS)"),
                     #font.family = "Arial",
                     #axis
                     xscale = "d_y",
                     pval = T,
                     surv.scale = "percent",
                     xlim = c(0, 2920), 
                     break.time.by=730.5,
                     xlab = "Time from diagnosis (years)",
                     #ylim = c(0, 0.05),
                     break.y.by=NULL,
                     ylab = "OS rate (%)",
                     #legend
                     legend = c(0.8,0.75),
                     legend.title = Gene_name,
                     legend.labs = c("High","Low"),
                     #risk table
                     risk.table = F,# Add risk table
                     
                     #字体
                     font.tickslab = c(10, "black"),
                     font.x = c(10, "black"),
                     font.y = c(10, "black"),
                     font.main = c(10, "black"),
                     font.legend = c(10, "black"),
)
gg_surv1$plot<-gg_surv1$plot+theme(plot.title = element_text(hjust = 0.5)) 


surv_fit<-survfit(Surv(PFI.Time , PFI) ~ gene_level,data= merged_matrix)
gg_surv2<-ggsurvplot(surv_fit,
                     conf.int = F,
                     #fun = "cumhaz",
                     linetype =  1, # Change line type by groups
                     size=0.5,
                     censor = F,
                     #surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_classic(),# Change ggplot2 theme
                     
                     palette = c("#E72C19", "#224FA2"),
                     title = paste(Study_name,"(PFS)"),
                     #font.family = "Arial",
                     #axis
                     xscale = "d_y",
                     pval = T,
                     surv.scale = "percent",
                     xlim = c(0, 2920), 
                     break.time.by=730.5,
                     xlab = "Time from diagnosis (years)",
                     #ylim = c(0, 0.05),
                     break.y.by=NULL,
                     ylab = "OS rate (%)",
                     #legend
                     legend = c(0.8,0.75),
                     legend.title = Gene_name,
                     legend.labs = c("High","Low"),
                     #risk table
                     risk.table = F,# Add risk table
                     
                     #字体
                     font.tickslab = c(10, "black"),
                     font.x = c(10, "black"),
                     font.y = c(10, "black"),
                     font.main = c(10, "black"),
                     font.legend = c(10, "black"),
)
gg_surv2$plot<-gg_surv2$plot+theme(plot.title = element_text(hjust = 0.5)) 





##
RNA_data<-read.table("/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/Exp_TPM.mat",sep="\t",row.names=1,header=T)


RNA_data<-RNA_data[c("THBS2","FAP","IGF1"),]
immune_mt<-data.frame(RNA_data)
immune_mt$Cell.Type<-rownames(immune_mt)
immune_mt_melt<-melt(immune_mt,value.name = "Per",id.vars="Cell.Type")
immune_mt_melt$Type<-"GN"
immune_mt_melt$Type[which(grepl("P",immune_mt_melt$variable))]<-"GCPM"
immune_mt_melt$Type[which(grepl("T",immune_mt_melt$variable))]<-"GC"

Immune_cell_list<-c("THBS2","FAP","IGF1")

color_list<-c(RColorBrewer::brewer.pal(7,"Set2")[1],c("#56B4E9", "#E69F00"))
#names(color_list)<-c("Normal", "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c("GN", "GC" ,"GCPM" )
immune_mt_melt$Type<-factor(immune_mt_melt$Type,levels=c("GN", "GC" ,"GCPM" ))
immune_mt_melt$variable<-gsub("N|P|T","",immune_mt_melt$variable)


for (i in 1:length(Immune_cell_list)){
  #i=43
  Immune_cell<-Immune_cell_list[i]
  #Immune_cell="SMS"
  select_table<-immune_mt_melt[which(immune_mt_melt$Cell.Type==Immune_cell),]
  #my_comparisons <- list( c("Normal", "Intestinal"), c("Normal", "Diffuse"), c("Normal", "Mixed"),c("Normal", "Metastatic") )
  p<-ggpaired(select_table, x="Type", y="Per", fill="Type",id = "variable",
              add="jitter",line.color = "gray", line.size = 0.5,
              palette=color_list,
              xlab="Region", 
              ylab="TPM", title = paste(Immune_cell,"in RNA-seq"),width = 0.8,
              legend.title=" ",show.legend = F) + 
    theme_classic()+
    stat_compare_means( label = "p.format",method = "wilcox.test",paired = F,
                        comparisons=list(c("GN", "GC"  ),c( "GC" ,"GCPM" ))) +#配对t检验
    theme(legend.position = 'none',
          axis.text = element_text(size = 10), 
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
          strip.background = element_blank()
    )
  assign(paste("tmp_p", i, sep = ""), p)
  print(p)
}

box_p3<-tmp_p1+tmp_p2+tmp_p3+plot_layout(ncol=3,nrow=1) 





test<-readRDS("5.monocle2/4.Endo/monocle2.RDS")
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
  ggtitle('Cell trajectory of PCs and ECs')
monocle2_time




Cluster<-c(  "EC_C1-CD320"  ,  "EC_C2-ACKR1" ,   "EC_C3-MAGI1"   , "EC_C4-INSR"   ,  "EC_C5-PLVAP",    "iCAF_C1-CXCL14",
             "mCAF_C1-THBS2",   "mCAF_C2-KRT8"  , "PC_C1-MYH11" ,   "PC_C2-RGS5"   ) 
g.colSet1 <- c(RColorBrewer::brewer.pal(8,"Set3"),
               RColorBrewer::brewer.pal(8,"Set2")[3:9])
names(g.colSet1)<-Cluster

test$Cluster<-as.character(test$Cluster)
cluster_list1<-c(  "PC_C1-MYH11" ,   "PC_C2-RGS5"   )
cluster_list2<-c(   "PC"  ,  "PC")
for (i in 1:length(cluster_list1)){
  test$SubCluster[which(test$SubCluster==cluster_list1[[i]])]<-cluster_list2[[i]]
  
}

test$SubCluster<-factor(test$SubCluster,levels=c("PC",
                                                 "CD320+ EC" ,"EC" ))
test$Type<-factor(test$Type,levels = c("PE","GCPM","GN","GC"))
monocle2_cluster<-plot_cell_trajectory(test,color_by = "SubCluster",
                                       cell_size =0.4)+
  theme(legend.position = "none")+
  scale_color_manual(values = RColorBrewer::brewer.pal(4,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))+
  ggtitle('Cell trajectory by clusters')+
  guides(colour = guide_legend(ncol=1,override.aes = list(size=3)))
monocle2_cluster


pseudotime_tab<-data.frame(Cell_ID=colnames(test),
                           Pseudotime=test$Pseudotime,
                           Type=test$Type,
                           Cluster=test$SubCluster,
                           State=test$State)

monocle2_density<-ggplot(pseudotime_tab,aes(x=Pseudotime,color=Cluster))+
  geom_density(size=0.8)+
  theme_classic2()+
  facet_grid(Cluster~.,switch = "y",scales = "free_y")+
  scale_color_manual(breaks = c("PC","CD320+ EC" , "EC" ),
                     labels = c("PC","CD320+ EC" , "EC" ),
                     values = RColorBrewer::brewer.pal(3,"Set2"))+
  ggtitle('Cell density by pseudotime')+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10), 
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position="none",
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))
monocle2_density





##TF activity
seu<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")
Idents(seu)<-seu$Cluster
seu<-subset(seu,idents=c("PC","EC"))
pseudotime_tab<-data.frame(Cell_ID=colnames(test),
                           Pseudotime=test$Pseudotime,
                           Cluster=test$SubCluster)
rownames(pseudotime_tab)<-pseudotime_tab$Cell_ID
seu<-AddMetaData(seu,metadata = pseudotime_tab)
DefaultAssay(seu)<-"RNA"
seu<-ScaleData(seu)

gene_list<-c("PDGFRB" , "MYH11","CD320","RGS5","FLT1","PLVAP","PECAM1","ACKR1")

select_TF_mt<-as.matrix(seu@assays$RNA@scale.data[gene_list,])
select_TF_mt<-reshape2::melt(select_TF_mt)
colnames(select_TF_mt)<-c("TF_names","Cell_ID","value")

merged_tab<-merge(pseudotime_tab,select_TF_mt,by="Cell_ID")
EXP_time<-ggplot(merged_tab,aes(x=Pseudotime,y=value,color=Cluster))+
  geom_smooth(se = F)+
  ggtitle("Expression of TF")+
  theme_classic2()+
  scale_color_manual(values = RColorBrewer::brewer.pal(3,"Set2"))+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA,color="grey"))+
  facet_wrap(~TF_names,scales="free_y",ncol=4)
EXP_time



total_p2<-umap_p2+nh_graph_stromal+beeswarm_stromal+box_p+plot_layout(ncol = 4,nrow = 1,
                                                                      widths  = c(1.8,1.8,1,2))


total_p1<-ggarrange(umap_p1,dot_p,ncol=2,widths = c(1,1.2))




total_p4<-ggarrange(box_p1,box_p,
                    nrow = 2,ncol = 1,heights = c(1,2))


total_p5<-ggarrange(dot_p,total_p4,
                    ncol = 2,nrow = 1,
                    widths = c(1,2))


surv_p<-gg_surv1$plot+gg_surv2$plot+plot_layout(ncol = 2,widths = c(1,1))
total_p9<-ggarrange(box_p3,surv_p,ncol = 2,widths = c(1,1)) 


letf_p<-monocle2_time+monocle2_cluster+plot_layout(ncol=1,nrow = 2,heights = c(1,1))

total_p10<-ggarrange(letf_p,monocle2_density,EXP_time,ncol=3,widths = c(0.7,0.45,1.25))

total_p<-ggarrange(total_p1,box_p1,total_p2,total_p9,total_p10,
                   ncol=1,nrow = 5,
                   heights = c(3.5,2,3,2.5,2.5))


dir.create("14.Figure_plot/")
ggsave("14.Figure_plot/Figure4.pdf",total_p,width = 12,height = 13.5)
