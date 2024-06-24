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

tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
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
  theme(axis.text = element_text(size=10))
cor_pp1 


tmp_table<-data.frame(Cluster=combined$SubCluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)

sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)

sum_table<-sum_table[which(sum_table$Cluster %in% c("Mph_C1-TIMP1"    ,
                                                    "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                    "Mph_C5_TAM-FBP1" ,      "iCAF_C1-CXCL14",
                                                    "mCAF_C1-THBS2",  "mCAF_C2-KRT8" )),]

cluster_list1<-c(   "Mph_C1-TIMP1"    ,
                    "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                    "Mph_C5_TAM-FBP1" ,      "iCAF_C1-CXCL14",
                    "mCAF_C1-THBS2",  "mCAF_C2-KRT8")
cluster_list2<-c(   "TIMP1+ Mph"    ,
                    "F13A1+ TRM" ,"CXCL9+ Mph", "SPP1+ TAM", 
                    "FBP1+ TAM" ,      "CXCL14+ iCAF",
                    "THBS2+ mCAF",  "KRT8+ mCAF")
for (i in 1:length(cluster_list1)){
  sum_table$Cluster[which(sum_table$Cluster==cluster_list1[[i]])]<-cluster_list2[[i]]
  
}
sum_table$Cluster<-factor(sum_table$Cluster,levels=cluster_list2)




library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p


library(ggcorrplot2)
cor_pp2<-ggcorrplot(cor_value, method = "ellipse",type = "lower",p.mat = cor_test_mat,col = c("#839EDB", "white", "#FF8D8D"),
                   insig = "label_sig", sig.lvl = c(0.05, 0.01, 0.001))

cor_pp2<-cor_pp2+
  theme(axis.text = element_text(size=10))
cor_pp2 





##Cor
tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
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
cell_pair_list<-list(c("mCAF","Mph"),
                     c("PC","EC"),
                     c("PC","Plasma"),
                     c("mCAF","Plasma"))

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
    geom_point(data = merged_table,aes(x=per.x, y=per.y),position=position_dodge(0),size=2,color="grey")+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(paste("Correlation of",Cluster1,"and",Cluster2)))+
    labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 10)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p2<-p1+p2+p3+p4+plot_layout(ncol=2)
cor_p2<-p1





library(CellChat)
library(patchwork)

cellchat1<-readRDS("6.cell_chat/1.CAF/GCPM.RDS")
cellchat2<-readRDS("6.cell_chat/1.CAF/GC.RDS")



cellchat <- mergeCellChat(list(cellchat1, cellchat2), add.names = c("GCPM", "GC"))



#cellchat@idents<-factor(cellchat@idents,levels=c( "B" ,"Plasma"   ,    "CD4_conv"   ,   "CD4_Treg",  "CD8", "MAIT"  ,   "NK"   ,    "DC"  , "Macro"  ,    "Mono"  ,  
#                                                  "Mast" ,    "EC", "APC","pAD"  ,  "SMC"            )) 

groupSize <- as.numeric(table(cellchat1@idents))
levels(cellchat1@idents) 
vertex.receiver = seq(1,4) 




pathways.show <- "COMPLEMENT"
select_list<-"CXCL|MDK|CSF|C3_"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])
cell_chat_p1<-netVisual_bubble(cellchat1, pairLR.use = pairLR.use,  sources.use = c(  "B"     ,      "CD4" ,   "CD8"  ,    "DC"    ,      "EC" ,      
                                                                                                            "Epi"    ,     "iCAF"      ,  "MAIT"     ,   "Mast"  ,      "mCAF"   ,     "Mono" ,       "Neutro"  ,   
                                                                                                            "NK"     ,   "NKT"  ,     "Mph"      ,    "PC" ), targets.use = "Mph", remove.isolate = FALSE)+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 8))


pathways.show <- "COMPLEMENT"
select_list<-"CXCL|MDK|PTN|IL|CSF|THBS|THY|C3|TEK|VEGF"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])

cell_chat_p2<-netVisual_bubble(cellchat, pairLR.use = pairLR.use,  comparison = c(1, 2) , sources.use = c(1:19), targets.use = 7, remove.isolate = FALSE)






library(CellChat)
library(patchwork)

cellchat1<-readRDS("6.cell_chat/3.mCAF_THBS2/mCAF_THBS2.RDS")

pathways.show <- "COMPLEMENT"
select_list<-"C3_"
pairLR.use = data.frame(interaction_name=cellchat1@LR$LRsig$interaction_name[grepl(select_list,cellchat1@LR$LRsig$interaction_name)])

cell_chat_p3<-netVisual_bubble(cellchat1, pairLR.use = pairLR.use,  sources.use = 1, targets.use = c(1:7), remove.isolate = FALSE)+
  coord_flip()+
  ggtitle("Cell interference")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,vjust=1,hjust=1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10))




#correlation
combined1<-readRDS("3.Cluster/13.Annotation/combined_annotation.rds")
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))


Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("CAF","Mph"))
combined$SubCluster<-as.character(combined$SubCluster)

cluster_list1<-c(   "Mph_C1-TIMP1"    ,
                    "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                    "Mph_C5_TAM-FBP1" ,      "iCAF_C1-CXCL14",
                    "mCAF_C1-THBS2",  "mCAF_C2-KRT8")
cluster_list2<-c(   "TIMP1+ Mph"    ,
                    "F13A1+ TRM" ,"CXCL9+ Mph", "SPP1+ TAM", 
                    "FBP1+ TAM" ,      "CXCL14+ iCAF",
                    "THBS2+ mCAF",  "KRT8+ mCAF")
for (i in 1:length(cluster_list1)){
  combined$SubCluster[which(combined$SubCluster==cluster_list1[[i]])]<-cluster_list2[[i]]
  
}
combined$SubCluster<-factor(combined$SubCluster,levels=rev(cluster_list2))



DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
data<-combined@assays$RNA@scale.data
data<-as.matrix(data)

gene_name<- c("C3")


select_fpkm_matrix1<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Efficiency_ICB),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]),
                               gene_name=gene_name)
gene_name<- c("C3AR1")
select_fpkm_matrix2<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Efficiency_ICB),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]),
                               gene_name=gene_name)
select_fpkm_matrix<-rbind(select_fpkm_matrix1,select_fpkm_matrix2)

box_p1<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0, #尾部修剪，数值越大修剪程度越高
                      scale= 2, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  ggtitle("Expression of C3 and C3AR1")+
  xlab("Expression level")+
  scale_fill_manual(values = RColorBrewer::brewer.pal(8,"GnBu"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  geom_hline(yintercept = seq(from=1, to=8, by = 1),color="grey")+
  xlim(c(-1 ,4))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.line = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        panel.border = element_rect(fill=NA,color="grey"),
        strip.background = element_blank())+
  facet_wrap(~gene_name,ncol=2,scales="free_x")
box_p1





Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))

Idents(combined)<-combined$Cluster

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mph","Mono","DC"))
combined$SubCluster<-factor(combined$SubCluster,levels=c( "DC_C1_pDC-GPR183"  ,   "DC_C2_cDC1-LAMP3"  ,   "DC_C3_cDC1-IRF8" ,  
                                                          "DC_C4_cDC2-CD1C"   ,  
                                                          
                                                          
                                                          "Mono_C1-FCN1"   ,  "Mono_C2-CX3CR1"  , "Mono_C3-HSPA1A"      , 
                                                          "Mph_C1-TIMP1"    ,
                                                          "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                          "Mph_C5_TAM-FBP1" ))


DefaultAssay(combined)<-"RNA"
combined<-ScaleData(combined)
data<-combined@assays$RNA@scale.data
data<-as.matrix(data)

gene_name<- c("C3AR1")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=factor(combined$Efficiency_ICB),
                               Cluster=combined$SubCluster,
                               gene=as.numeric(data[gene_name,]))


box_p2<-ggplot(data=select_fpkm_matrix,aes(x=gene,y=Cluster,fill=Cluster))+
  geom_density_ridges(alpha = 0.8,
                      #color= 'white',
                      rel_min_height= 0.01, #尾部修剪，数值越大修剪程度越高
                      scale= 1.8, #山脊重叠程度调整，scale = 1时刚好触及基线，数值越大重叠度越高
                      quantile_lines= TRUE, #显示分位数线
                      quantiles= 2 ) + 
  ggtitle("Expression of component C3AR1")+
  xlab("Expression level")+
  scale_fill_manual(values = viridis::viridis(12,option = "D"))+
  #stat_compare_means(method = "wilcox.test",label = "p.format",size=3)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlim(c(-1 ,4))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
box_p2



##Expression
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

markers<-c("C3","C3AR1")
Vln_P1<-VlnPlot(combined,features = markers,pt.size = 0,combine = F,group.by = "Cluster")

ps_list<-list()
k=0
for (i in Vln_P1){
  i<-i+theme(plot.title = element_text(size = 10,hjust = 0.5),
             axis.title.x = element_text(size = 10),
             axis.title.y = element_text(size = 10),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             legend.title = element_text(size = 10),
             
             legend.text = element_text(size = 10),
             legend.position = "none")
  
  
  k=k+1
  
  if (k<=1){
    i<-i+theme(axis.title.x = element_blank(),
               axis.text.x = element_blank())
  }
  
  ps_list[[k]]<-i
}
Vln_P1<-ps_list[[1]]+ps_list[[2]]+plot_layout(ncol=1,nrow = 2)





markers<-c("C3","C3AR1")
Feature_P1<-FeaturePlot(combined,features = markers,raster=F,ncol=1,max.cutoff =2) &
  scale_color_gradientn(colors=c("white",RColorBrewer::brewer.pal(8,"Purples"))) &
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 10),
        panel.border=element_rect(fill=NA,color="grey"),
        legend.text = element_text(size = 10),
        legend.position = "none")

Feature_P1



Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("CAF","Mph"))
combined$SubCluster<-factor(combined$SubCluster,levels=c("Mph_C1-TIMP1"    ,
                                                         "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9", "Mph_C4_TAM-SPP1", 
                                                         "Mph_C5_TAM-FBP1" ,          "iCAF_C1-CXCL14",
                                                         "mCAF_C1-THBS2",  "mCAF_C2-KRT8"    ))



markers<-c("C3")
Vln_P1<-VlnPlot(combined,features = markers,pt.size = 0,combine = F,group.by = "SubCluster",same.y.lims = T)

#ps_list<-list()
for (i in Vln_P1){
  i<-i+theme(plot.title = element_text(size = 10,hjust = 0.5),
             axis.title.x = element_text(size = 10),
             axis.title.y = element_text(size = 10),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             legend.title = element_text(size = 10),
             
             legend.text = element_text(size = 10),
             legend.position = "none")
  
  
  k=k+1
  
  if (k<=3){
    i<-i+theme(axis.title.x = element_blank(),
               axis.text.x = element_blank())
  }
  
  ps_list[[k]]<-i
}
Vln_P2<-ps_list[[3]]



Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mph","Mono"))
markers<-c("C3AR1","SPP1")
Vln_P1<-VlnPlot(combined,features = markers,pt.size = 0,combine = F,group.by = "SubCluster",same.y.lims = T)
k=3
#ps_list<-list()
for (i in Vln_P1){
  i<-i+theme(plot.title = element_text(size = 10,hjust = 0.5),
             axis.title.x = element_text(size = 10),
             axis.title.y = element_text(size = 10),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             legend.title = element_text(size = 10),
             
             legend.text = element_text(size = 10),
             legend.position = "none")+
    coord_flip()
  
  
  k=k+1
  
  if (k>=5){
    i<-i+theme(axis.title.y = element_blank(),
               axis.text.y = element_blank())
  }
  
  ps_list[[k]]<-i
}
Vln_P2<-ps_list[[4]]+ps_list[[5]]+plot_layout(ncol=2,nrow = 1,heights = c(1.1,1))

Vln_P1<-ps_list[[1]]+ps_list[[2]]+plot_layout(ncol=1,nrow = 2)









##correlation of subcluster
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM","GN"))
combined$Cluster<-as.character(combined$Cluster)
combined$SubCluster<-as.character(combined$SubCluster)
combined$Cluster[which(combined$Cluster %in% c("CAF","Mph","Mono"))]<-combined$SubCluster[which(combined$Cluster %in% c("CAF","Mph","Mono"))]

tmp_table<-data.frame(Cluster=combined$Cluster,Type=combined$Type,sample_ID=combined$Sample_ID,num=1)
tmp_table<-na.omit(tmp_table)
sum_table1<-dplyr::summarise(group_by(tmp_table,Cluster,sample_ID),num=sum(num))
sum_table1 <-data.frame(sum_table1)
sum_table1 <- tidyr::complete(sum_table1,Cluster,sample_ID, fill = list(num = 0))
sum_table2<-dplyr::summarise(group_by(tmp_table,Type,sample_ID),total_num=sum(num))
sum_table<-merge(sum_table1,sum_table2,by="sample_ID")
sum_table$per<-sum_table$num/sum_table$total_num*100

sum_table$Cluster<-gsub("\\/","_",sum_table$Cluster)


cluster_list1<-c(   "mCAF_C1-THBS2"  , "iCAF_C1-CXCL14" ,"Mph_C1-TIMP1"    ,
                    "Mph_C2_TRM-F13A1" ,"Mph_C3-CXCL9",
                     "Mph_C4_TAM-SPP1")
cluster_list2<-c(   "THBS2+ mCAF"  , "CXCL14+ iCAF" ,"TIMP1+ Mph",
                    "F13A1+ TRM" ,"CXCL9+ Mph",
                    "SPP1+ TAM")
for (i in 1:length(cluster_list1)){
  sum_table$Cluster[which(sum_table$Cluster==cluster_list1[[i]])]<-cluster_list2[[i]]
  
}


library(reshape2)
cor_data<-dcast(sum_table,Cluster~sample_ID,value.var = "per")
rownames(cor_data)<-cor_data$Cluster
cor_data<-cor_data[,-1]
cor_data<-t(cor_data)
cor_value<-cor(cor_data)

library("psych")
cor_test_mat <- corr.test(cor_data)$p



cell_pair_list<-list(c("THBS2+ mCAF","TIMP1+ Mph"),
                     c("THBS2+ mCAF", "F13A1+ TRM"),
                     c("THBS2+ mCAF","CXCL9+ Mph"),
                     c("THBS2+ mCAF", "SPP1+ TAM"))

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
    geom_point(data = merged_table,aes(x=per.x, y=per.y),position=position_dodge(0),size=2,color="grey")+
    stat_cor()+
    geom_smooth(size=0.5,method="lm",se=F)+
    scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
    ylim(c(0,25))+
    xlim(c(0,25))+
    ylab(Cluster2)+
    xlab(Cluster1)+
    #expand_limits(y = c(0,100))+
    ggtitle(paste(Cluster1,"~",Cluster2))+
    theme_classic()+
    theme(legend.position="none",
          plot.title = element_text(size = 10,hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 10)) +
    scale_y_continuous(expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0))
  assign(paste("p", i, sep = ""), p)
  print(p)
}
cor_p4<-p1+p2+p3+p4+plot_layout(ncol=2)












vocano_plot = function(sorted_matrix, Sample_1 = "A", Sample_2 = "B", lfc = 0, pval = 0.05,gene_list=NA){
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
  xmax <- max(abs(DF$logFC))+0.3
  ymax <- ceiling(max(abs(DF$negLogPval)))*1.1
  
  
  DF$value<-abs(DF$logFC)
  DF<-DF[order(DF$negLogPval,decreasing = T),]
  DF$Gene<-""
  if (nrow(DF)>10){
    DF$Gene[DF$logFC>0][1:5]<-rownames(DF)[DF$logFC>0][1:5]
    DF$Gene[DF$logFC<0][1:5]<-rownames(DF)[DF$logFC<0][1:5]
  }else{
    DF$Gene<-rownames(DF)
  }
  DF$Gene[rownames(DF) %in% gene_list]<-rownames(DF)[rownames(DF) %in% gene_list]
  
  
  p <- ggplot(DF, aes(x = logFC, y = negLogPval, label=DF$Gene)) +
    geom_point(aes(color = SigGenes))+ 
    geom_text_repel(max.overlaps = 2000)+
    xlim(-xmax,xmax) + ylim(0,ymax) +
    scale_color_manual(values = c("#B31B21", "#1465AC","grey")) +
    theme_classic(base_size = 10) + theme(legend.position = "bottom") +
    xlab("Log2FoldChange")+ylab(paste("-log10"," FDR",sep=""))+
    geom_vline(aes(xintercept=-lfc),colour="darkgrey", linetype="dashed")+
    geom_vline(aes(xintercept=lfc),colour="darkgrey", linetype="dashed") +
    geom_hline(aes(yintercept=-log10(pval)),colour="darkgrey", linetype="dashed")+
    ggtitle(paste(Sample_2,Sample_1,sep=paste(rep(" ",15),collapse=""))) +
    #expression("Group 1" %->% "Group 2"),
    annotate("text", x=-xmax*0.9, y=-log10(pval), label= paste("FDR"," < ",pval,sep=""))+
    #annotate("text", x=0, y=-log10(pval), label= "2fold")+
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}

libs <- c("limma","ggpubr","ggrepel","reshape2","cluster","splines","ggplot2","gridExtra")
loaded <- sapply(libs, library, character.only = T)
fc = 2
lfc = log2(fc)
pval = 0.05
#vocano_plot(sorted_ori_de, Sample_1 = 'shKDM4C', Sample_2 = 'shControl', lfc = lfc, pval = pval)



##PT  TVAT

Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GCPM","GC"))

Idents(combined)<-combined$Cluster
combined<-subset(combined,idents="CAF")


C1<-"GCPM"
C2<-"GC"
Idents(combined)<-combined$Type
sub_sce<-subset(combined,idents=c(C1,C2))
markers <- FindAllMarkers(sub_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
markers1<-markers$gene[which(markers$cluster==C1)]
markers$avg_log2FC[which(markers$cluster==C2)]<-(-(markers$avg_log2FC[which(markers$cluster==C2)]))
vanc_p1<-vocano_plot(markers, Sample_1 = C1, Sample_2 = C2, lfc = 0, pval = pval,gene_list=c("C3","PI16","POSTN"))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10))






##Expresion compare
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM"))
Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("CAF"))

data<-combined@assays$integrated@scale.data
data<-as.matrix(data)

gene_name<- c("C3")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=combined$Type,
                               Cluster=factor(as.character(combined$Cluster)),
                               gene=as.numeric(data[gene_name,]))

#select_fpkm_matrix<-merge(select_fpkm_matrix,data1,by="Cell_ID")


C3_p1<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_violin(data=select_fpkm_matrix,aes(x=Type,group=Type,fill=Type),trim=FALSE,color="white") + 
  geom_boxplot(data=select_fpkm_matrix,aes(x=Type,group=Type,fill=Type),width=0.2,color="black")+
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  ggtitle("C3 in mCAF")+
  ylab("Expression level")+
  stat_compare_means(method = "t.test",label = "p.format",size=3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  #ylim(c(-1 ,5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
C3_p1








##Expresion compare
Idents(combined1)<-combined1$Type
combined<-subset(combined1,idents=c("GC","GCPM"))
Idents(combined)<-combined$Cluster
combined<-subset(combined,idents=c("Mph"))

data<-combined@assays$integrated@scale.data
data<-as.matrix(data)

gene_name<- c("C3AR1")


select_fpkm_matrix<-data.frame(Cell_ID=combined$Cell_ID,
                               Type=combined$Type,
                               Cluster=factor(as.character(combined$Cluster)),
                               gene=as.numeric(data[gene_name,]))

#select_fpkm_matrix<-merge(select_fpkm_matrix,data1,by="Cell_ID")


C3AR1_p1<-ggplot(data=select_fpkm_matrix,aes(x=Type,y=gene,group=Type,color=Type))+
  geom_violin(data=select_fpkm_matrix,aes(x=Type,group=Type,fill=Type),trim=FALSE,color="white") + 
  geom_boxplot(data=select_fpkm_matrix,aes(x=Type,group=Type,fill=Type),width=0.2,color="black")+
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  ggtitle("C3AR1 in TAM")+
  ylab("Expression level")+
  stat_compare_means(method = "t.test",label = "p.format",size=3)+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  #ylim(c(-1 ,5))+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 10),
        strip.background = element_blank())
C3AR1_p1




##
RNA_data<-read.table("/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/Exp_TPM.mat",sep="\t",row.names=1,header=T)


RNA_data<-RNA_data[c("C3","C3AR1"),]
immune_mt<-data.frame(RNA_data)
immune_mt$Cell.Type<-rownames(immune_mt)
immune_mt_melt<-melt(immune_mt,value.name = "Per",id.vars="Cell.Type")
immune_mt_melt$Type<-"GN"
immune_mt_melt$Type[which(grepl("P",immune_mt_melt$variable))]<-"GCPM"
immune_mt_melt$Type[which(grepl("T",immune_mt_melt$variable))]<-"GC"

Immune_cell_list<-unique(sort(immune_mt_melt$Cell.Type))

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
              ylab="TPM", title = paste(Immune_cell,"in RNA-seq"),
              legend.title=" ",show.legend = F,width=0.8) + 
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
  assign(paste("p", i, sep = ""), p)
  print(p)
}

box_p3<-p1+p2+plot_layout(ncol=1,nrow=2) 


RNA_data<-read.table("/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/Exp_TPM.mat",sep="\t",row.names=1,header=T)
RNA_data<-RNA_data[c("C3","C3AR1","CXCL1","CXCL8","CXCL9","CXCL10","CXCL14","CCL1","CCL2","CCL11","CCL19","CCL20","CCL22","IL1B","SPP1","TGFB1"),]
RNA_data<-data.frame(RNA_data[,grepl("N",colnames(RNA_data))],
                     RNA_data[,grepl("T",colnames(RNA_data))],
                     RNA_data[,grepl("P",colnames(RNA_data))])

colData1<-data.frame(Type=c(rep("GN",12),rep("GC",12),rep("GCPM",12)))
rownames(colData1)<-colnames(RNA_data)

library(pheatmap)
pheat_p1<-pheatmap(RNA_data,
                   annotation_col = colData1,
                   cluster_rows = F,
                   cluster_cols = F,
                   show_colnames = F,
                   scale = "row",
                   color=colorRampPalette(c('#3C5488FF','#3C5488FF','white','#DC0000FF','#DC0000FF'), bias=1)(50),border_color = NA,
                   main = "Expression in bulk RNA data",
                   legend.position = "bottom",
                   fontsize = 10)
pheat_p1<-plot_grid(pheat_p1$gtable)









immune_mt<-read.table("/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/Exp_TPM.mat",sep="\t",row.names=1,header=T)
#immune_mt<-data.frame(Cell.type=rownames(immune_mt),immune_mt)



##VEGFR
Cluster1<-"C3"
Cluster2<-"C3AR1"


CD8_tab<-data.frame(
  sample_ID=names(immune_mt),
  per=as.numeric(immune_mt[Cluster1,]))


Neutro_tab<-data.frame(
  sample_ID=names(immune_mt),
  per=as.numeric(immune_mt[Cluster2,]))




merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")

cor.value<-cor(merged_table$per.x,merged_table$per.y)
cor.value<-round(cor.value,2)
p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
if (p.value==0){
  p.value<-"< 10e-16"
}else if (p.value<0.001){
  p.value<-signif(p.value,digits = 2)
}else{
  p.value<-round(p.value,3)
}

RNA_cor_p1 <- ggplot(merged_table, aes(x=per.x, y=per.y))+
  geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=2,color="#374E55")+
  geom_smooth(size=0.5,method="lm",se=F)+
  scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
  ylim(c(0,25))+
  xlim(c(0,25))+
  ylab(Cluster2)+
  xlab(Cluster1)+
  #expand_limits(y = c(0,100))+
  ggtitle(paste(paste("Cor of",Cluster1,"and",Cluster2)))+
  labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
  theme_classic()+
  theme(legend.position="none",
        plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
print(RNA_cor_p1)


##VEGFR
Cluster1<-"C3AR1"
Cluster2<-"SPP1"


CD8_tab<-data.frame(
  sample_ID=names(immune_mt),
  per=as.numeric(immune_mt[Cluster1,]))


Neutro_tab<-data.frame(
  sample_ID=names(immune_mt),
  per=as.numeric(immune_mt[Cluster2,]))




merged_table<-merge(CD8_tab,Neutro_tab,by="sample_ID")

cor.value<-cor(merged_table$per.x,merged_table$per.y)
cor.value<-round(cor.value,2)
p.value<-cor.test(merged_table$per.x,merged_table$per.y)$p.value
if (p.value==0){
  p.value<-"< 10e-16"
}else if (p.value<0.001){
  p.value<-signif(p.value,digits = 2)
}else{
  p.value<-round(p.value,3)
}

RNA_cor_p2 <- ggplot(merged_table, aes(x=per.x, y=per.y))+
  geom_point(data = merged_table,aes(x=per.x, y=per.y),pch=15,position=position_dodge(0),size=2,color="#374E55")+
  geom_smooth(size=0.5,method="lm",se=F)+
  scale_color_manual(values = c("#55752F","#90343B","#1A476F"))+
  ylim(c(0,25))+
  xlim(c(0,25))+
  ylab(Cluster2)+
  xlab(Cluster1)+
  #expand_limits(y = c(0,100))+
  ggtitle(paste(paste("Cor of",Cluster1,"and",Cluster2)))+
  labs(subtitle = paste("R = ",cor.value,", p = ",p.value,sep=""))+
  theme_classic()+
  theme(legend.position="none",
        plot.title = element_text(size = 10,hjust = 0.5,face = "bold"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        
        legend.text = element_text(size = 10)) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))
print(RNA_cor_p2)


RNA_cor<-RNA_cor_p1+RNA_cor_p2+plot_layout(ncol=1)






##TF activity
library(monocle3)
cds<-readRDS("5.monocle2/7.TAM_monocle3/TAM_monocle3.rds")
cds_subset <- cds[, 
                  colData(cds) %>%
                    subset(
                      SubCluster %in% c("Mph_C4_TAM-SPP1","Mph_C3-CXCL9") &
                      Type %in% c("GCPM","GC")
                    ) %>%
                    row.names
]

cds_p3<-plot_cells(cds_subset, genes = "C3AR1", label_cell_groups = FALSE, scale_to_range = F,trajectory_graph_color = "grey",
                   label_leaves = FALSE,  label_branch_points = FALSE)+
  theme_classic()+
  xlim(c(-10,0))+
  ylim(c(-5,5))+
  ggtitle("Expression of C3AR1")+
  scale_color_gradientn(colors =viridis::viridis(5,option="D") )+
  facet_grid(SubCluster~Type,switch = "y")+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        panel.border = element_rect(fill = NA,color="grey"))
cds_p3


##TF activity
library(monocle)
seu<-readRDS("7.scenic/TAM/seu_obj.RDS")
monocle1<-readRDS("5.monocle2/3.Mono/monocle2.RDS")

pseudotime_tab<-data.frame(Cell_ID=colnames(monocle1),
                           Pseudotime=monocle1$Pseudotime,
                           Cluster=monocle1$Cluster,
                           Type=monocle1$Type)
rownames(pseudotime_tab)<-pseudotime_tab$Cell_ID
seu<-AddMetaData(seu,metadata = pseudotime_tab)



gene_list<-c("C3AR1","MAF")

DefaultAssay(seu)<-"RNA"
seu<-ScaleData(seu)

select_TF_mt<-as.matrix(seu@assays$RNA@scale.data[gene_list,])
select_TF_mt<-reshape2::melt(select_TF_mt)
colnames(select_TF_mt)<-c("TF_names","Cell_ID","value")

merged_tab<-merge(pseudotime_tab,select_TF_mt,by="Cell_ID")
merged_tab<-merged_tab[merged_tab$Cluster %in% c("TAM" ,"CXCL9+ Mph"),]
merged_tab<-merged_tab[merged_tab$TF_names=="C3AR1",]
EXP_time<-ggplot(merged_tab,aes(x=Pseudotime,y=value,color=Cluster))+
  geom_smooth(se = F)+
  ggtitle("C3AR1 in TRM-derived Mph trajectory")+
  theme_classic2()+
  scale_color_manual(breaks = c("CXCL9+ Mph","TAM"),
                     labels = c("CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2")[2:3])+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = NA,color="grey"))
EXP_time

merged_tab<-merged_tab[merged_tab$Type %in% c("GCPM","GC"),]
EXP_time1<-ggplot(merged_tab,aes(x=Pseudotime,y=value,color=Cluster))+
  #geom_point(alpha=0.5,size=0.1)+
  geom_smooth(se = F)+
  theme_classic2()+
  scale_color_manual(breaks = c("CXCL9+ Mph","TAM"),
                     labels = c("CXCL9+ Mph","TAM"),
                     values = RColorBrewer::brewer.pal(3,"Set2")[2:3])+
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.title.x = element_text(size = 10),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.line.y = element_blank(), 
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = NA,color="grey"))+
  facet_wrap(~Type,scales="free_y",ncol=2)
EXP_time1

C3AR1_time<-ggarrange(EXP_time,EXP_time1,ncol=1,nrow = 2)



##cell_migration
data1<-read.table("cell_migration/Cell_num_migration.txt",sep="\t",header = T)
data1$mean.var<-apply(data1[,3:7], 1, mean)
data1$sd<-apply(data1[,3:7], 1, sd)
data1$lower.ci<-data1$mean.var-data1$sd
data1$upper.ci<-data1$mean.var+data1$sd
data1$concentration<-gsub("/ml","",data1$concentration)
data1$concentration<-factor(data1$concentration,levels=c("0ng" , "5ng" , "10ng" ,"20ng"))
HR_p1 <- ggplot(data1, aes(x=concentration, y=mean.var,group=1))+
  geom_bar(stat = "identity",size = 0.1,color="grey",fill="#4DBBD5")+
  #geom_line(size=0.5,position=position_dodge(0),color="#374E55")+
  geom_errorbar(aes(ymin = lower.ci, ymax=upper.ci), #误差条表示均值±标准差
                width=0.1, #误差条末端短横线的宽度
                #position=position_dodge(0), 
                color="#374E55",
                alpha = 0.7,
                size=0.5) +
  
  ylab("Cell number")+
  xlab("Concentration of C3 (ng/ml)")+
  theme_classic2()+
  theme(
    plot.title = element_text(size = 12,hjust = 0.5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    strip.background = element_blank(),
    legend.text = element_text(size = 10))+
  scale_y_continuous(limits=c(0,NA))+
  facet_wrap(~Cell_line,scales = "free",ncol=1)
HR_p1


low_p<-ggarrange(Feature_P1,Vln_P1,ncol=2,widths=c(1,2.5))
cell_chat_p<-ggarrange(cell_chat_p1,low_p,ncol=1,nrow=2,heights = c(1,1.5))
  

total_p1<-ggarrange(cor_pp1,cell_chat_p,ncol = 2,widths = c(1,1))
total_p2<-cell_chat_p1+cell_chat_p2+plot_layout(ncol=2)
total_p3<-ggarrange(cor_pp2,box_p1+
                      theme(axis.title.y = element_blank(),
                            axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)),
                    cor_p4,
                    ncol = 3,widths = c(1,1,1))

total_p5<-ggarrange(cell_chat_p3,cds_p3,C3AR1_time,ncol=3,widths = c(1,1,1))


total_p4<-ggarrange(RNA_cor,pheat_p1,HR_p1,nrow=1,ncol=4,widths = c(0.6,0.8,0.5,1.1))


total_p<-ggarrange(total_p1,total_p3,total_p5,total_p4,ncol=1,heights = c(5.5,4,3,4))

dir.create("14.Figure_plot/")
ggsave("14.Figure_plot/Figure5.pdf",total_p,width = 12,height = 16.5)

ggsave("14.Figure_plot/Figure5_1.pdf",cor_p2,width = 2.8,height = 2.5)


