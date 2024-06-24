
library(CellChat)
library(patchwork)
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



##correlation
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))
combined$Cluster<-as.character(combined$Cluster)
combined$Cluster[grepl("iCAF",combined$SubCluster)]<-"iCAF"
combined$Cluster[grepl("mCAF",combined$SubCluster)]<-"mCAF"
combined$Cluster<-factor(combined$Cluster,levels=c( "B"    ,    "Plasma"  ,     "CD4"   ,   "CD8"     ,   "NK"   ,   "NKT",  "MAIT"     , "Tgd"     ,
                                                    "Mast"   ,  "DC"  ,              
                                                    "Mono"       ,            
                                                    "Mph"   ,  "Neutro" ,          
                                                    "iCAF","mCAF", "PC" , "EC" ))

tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

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


library(ggcorrplot2)
cor_pp1<-ggcorrplot(cor_value, method = "ellipse",type = "lower",p.mat = cor_test_mat,col = c("#839EDB", "white", "#FF8D8D"),
                    insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001))

cor_pp1<-cor_pp1+
  theme(axis.text = element_text(size=8))
cor_pp1 






total_p1<-cor_pp1

pdf("14.Figure_plot/Figure5_S8_correlation.pdf",width = 12,height = 13)
print(total_p1)
dev.off()
