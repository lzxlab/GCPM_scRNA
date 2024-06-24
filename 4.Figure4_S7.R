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
combined<-readRDS("3.Cluster/13.Annotation/5.Stromal_sub_sce_annotation.rds")
Idents(combined)<-combined$Cluster

markers<-c(
  "PECAM1","FLT1","ACKR1","INSR",
  
  "CXCL1","CXCL2","CXCL9","CXCL10","CXCL14","CCL5","CCL11","CCL21",
  
  "MGP","COL1A1","COL2A1","LUM","CLU","IGFBP4","FAP","YAP1","ACTA2",
  
  "CD34","PI16","APOD","APOE",
  "ACTC1","ACTG2","MYL9","MYH11","THY1","RGS5",
  
  "EPCAM","KRT18","KRT19","KRT7"
  
  
)  

DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
combined1<-subset(combined, features =markers )
combined1$SubCluster<-factor(combined1$SubCluster,levels=c(   "Endo_C1-ACKR1" , "Endo_C2-PLVAP" , "Endo_C3-MAGI1",  "Endo_C4-RGS5"  ,
                                                              "iCAF_C1-CXCL14",
                                                              "mCAF_C1-THBS2" , "mCAF_C2-TGFBI" , "mCAF_C3-SLIT2" ,
                                                              "PC_C1-MYH11"  ,  "PC_C2-RGS5" 
                                                              
))   

Exp_table<-data.frame(Cluster=combined1$SubCluster,t(combined1@assays$RNA@counts))
ave_Exp_table<-aggregate(.~Cluster,data=Exp_table,mean)
rownames(ave_Exp_table)<-ave_Exp_table$Cluster
ave_Exp_table<-ave_Exp_table[,-1]
ave_Exp_table<-t(data.frame(ave_Exp_table))
rownames(ave_Exp_table)<-gsub("\\.","-",rownames(ave_Exp_table))

colData1<-data.frame(Cluster=c(rep("Endo",4),rep("iCAF",4),rep("PC",2)))
rownames(colData1)<-colnames(ave_Exp_table)

rowData1<-data.frame(Markers=c(
  rep("ENdothelial",4),
  rep("Cytokines",8),
  rep("CAF markers",9),
  rep("MSC markers",4),
  rep("Peicyte/SMC markers",6),
  rep("EMT markers",4)
  
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
                       gaps_col=c(4,8),
                       gaps_row = c(4,12,21,25,31)
)
print(pheatmap_p1)
pheatmap_p1<-plot_grid(pheatmap_p1$gtable)

cellchat<-readRDS("6.cell_chat/4.CAF_IGF/mCAF_THBS2.RDS")
pathways.show <- "THBS"
select_list<-"THBS2"
groupSize <- as.numeric(table(cellchat@idents))
levels(cellchat@idents) 
vertex.receiver = seq(1,4) # a numeric vector
pairLR.use = data.frame(interaction_name=cellchat@LR$LRsig$interaction_name[grepl(select_list,cellchat@LR$LRsig$interaction_name)])


cell_chat1<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,   sources.use = 2, targets.use = c(1:16), remove.isolate = FALSE)




pathways.show <- "IGF"
select_list<-"IGF1|IGF2"
groupSize <- as.numeric(table(cellchat@idents))
levels(cellchat@idents) 
vertex.receiver = seq(1,4) # a numeric vector
pairLR.use = data.frame(interaction_name=cellchat@LR$LRsig$interaction_name[grepl(select_list,cellchat@LR$LRsig$interaction_name)])

cell_chat2<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,   sources.use = 2, targets.use = c(1:16), remove.isolate = FALSE)




##Enrichment

##GO_BP CD8
Enrich_type<-"GO_BP"
Type="Stromal"
Select_cluster_list<-c( "iCAF_C1-CXCL14",
                        "mCAF_C1-THBS2" , "mCAF_C2-TGFBI" , "mCAF_C3-SLIT2"  )
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/4.stromal/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
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
Enrich_type<-"KEGG"
Type="Stromal"
Select_cluster_list<-c( "iCAF_C1-CXCL14",
                        "mCAF_C1-THBS2" , "mCAF_C2-TGFBI" , "mCAF_C3-SLIT2"  )
total_table<-data.frame()
enriched_pathway<-c()
for (cluster in Select_cluster_list){
  path<-paste("4.characteristics/6.Enrichment/4.stromal/",cluster,"/",Enrich_type,"_enrich.xls",sep="")
  enrich_table<-data.table::fread(path)
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





##survival

library(ggpubr)
library(magrittr)
library(ggsignif)
library(ggplot2)
library(ggsci)
library(RColorBrewer)




Gene_name="THBS2_CAF"
select_gene=Gene_name
outdir="3.Cluster/13.Plot/"
setwd("/home/zhengyq/data/single_cell/3.GCPM/")

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

project_list<-c("TCGA-KIRP","TCGA-KIRC","TCGA-BLCA","TCGA-LGG","TCGA-PAAD","TCGA-THCA","TCGA-HNSC","TCGA-LUSC","TCGA-SKCM","TCGA-OV")


sur_plot_list<-list()
for (i in 1:length(project_list)){
  Study_name<-project_list[[i]]
  merged_matrix<-total_matrix[total_matrix$project_id==Study_name,]
  res.cut <- surv_cutpoint(merged_matrix, #数据集
                           time = "OS.Time", #生存状态
                           event = "OS", #生存时间
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
                       title = paste(Study_name),
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
                       legend = c(0.2,0.35),
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
  
  
  sur_plot_list[[length(sur_plot_list)+1]]<-gg_surv1$plot
  print(gg_surv1)
}






total_p3<-ggarrange(plotlist =sur_plot_list,ncol = 5,nrow = 2,
                    common.legend = T,legend = "bottom")





right_P<-cell_chat1+cell_chat2+plot_layout(ncol=1,nrow=2)
total_p1<-ggarrange(pheatmap_p1,right_P,ncol=2)
total_p2<-S1+S2+plot_layout(ncol=2,nrow=1)

total_p<-ggarrange(total_p1,total_p2,total_p3,nrow=3,ncol = 1,heights = c(6,6,4.5))
pdf("14.Figure_plot/Figure4_S7.pdf",width = 12,height = 16.5)
print(total_p)
dev.off()
