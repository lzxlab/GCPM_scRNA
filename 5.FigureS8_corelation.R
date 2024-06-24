library(ggpubr)
library(magrittr)
library(ggsignif)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

setwd("/home/zhengyq/data/single_cell/3.GCPM/")






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





##correlation
mt<-es.dif[1,]/es.dif[2,]
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=es.dif[1,],
                         gene2=es.dif[2,])


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"|merge_matrix$sample_type=="Solid Tissue Normal"),]
merge_matrix<-data.frame(merge_matrix)


select_matrix<-merge_matrix
project_id_list<-unique(select_matrix$project_id)
project_id_list<-project_id_list[which(!(project_id_list %in% c("TCGA-CHOL","TCGA-DLBC")))]
select_table<-select_matrix[select_matrix$project_id %in% project_id_list, ]


#绘图
TCGA_cor_p1 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor(size=3)+
  #ylim(c(0,5))+
  ylab(paste0("ssGSEA score (",Gene_name2,")"))+
  xlab(paste0("ssGSEA score (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of THBS2+ CAF and SPP1+ TAM, TCGA pan-caner dataset (N = 10,327; RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free",ncol=8)
print(TCGA_cor_p1)




##C3 C3AR1
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


Gene_name1="C3"
Gene_name2="C3AR1"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=as.numeric(RNA_matrix[Gene_name2,]))


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"),]
merge_matrix<-data.frame(merge_matrix)


select_matrix<-merge_matrix
project_id_list<-unique(select_matrix$project_id)
project_id_list<-project_id_list[which(!(project_id_list %in% c("TCGA-CHOL","TCGA-DLBC")))]
select_table<-select_matrix[select_matrix$project_id %in% project_id_list, ]
RColorBrewer::brewer.pal(3,"Set3")[1]
#绘图
TCGA_cor_p2 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#8DD3C7",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor()+
  #ylim(c(0,5))+
  ylab(paste0("Expression level (",Gene_name2,")"))+
  xlab(paste0("Expression level (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of C3 and C3AR1 in TCGA pan-cancer dataset (N = 10,327; RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free")
print(TCGA_cor_p2)



##C3AR1 SPP1
Gene_name1="C3AR1"
Gene_name2="SPP1"
select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=as.numeric(RNA_matrix[Gene_name1,]),
                         gene2=as.numeric(RNA_matrix[Gene_name2,]))


merge_matrix<-merge(phenotype_matrix,select_table,by="sample")

merge_matrix<-merge_matrix[which(merge_matrix$program=="TCGA"),]
merge_matrix<-merge_matrix[which(merge_matrix$sample_type=="Primary Tumor"|merge_matrix$sample_type=="Solid Tissue Normal"),]
merge_matrix<-data.frame(merge_matrix)


select_matrix<-merge_matrix
project_id_list<-unique(select_matrix$project_id)
project_id_list<-project_id_list[which(!(project_id_list %in% c("TCGA-CHOL","TCGA-DLBC")))]
select_table<-select_matrix[select_matrix$project_id %in% project_id_list, ]


#绘图
TCGA_cor_p3 <- ggplot(select_table, aes(x=gene1, y=gene2, group = 1))+
  geom_point(data = select_table,aes(x=gene1, y=gene2),size=0.5,color="#BEBADA",alpha=0.8)+
  geom_smooth(method="lm",size=0.5,se=F,color="black",linetype="dashed")+
  stat_cor()+
  #ylim(c(0,5))+
  ylab(paste0("Expression level (",Gene_name2,")"))+
  xlab(paste0("Expression level (",Gene_name1,")"))+
  #expand_limits(y = c(0,5),x=c(2,8))+
  ggtitle("Correlation of C3AR1 and SPP1 in TCGA pan-cancer dataset (N = 10,327; RNA-seq)")+
  theme_classic()+
  theme(plot.title = element_text(size = 12,hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text = element_blank(),
        plot.subtitle =element_text(size = 10,hjust = 0.0) ,
        strip.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  facet_wrap(~project_id,scales = "free")
print(TCGA_cor_p3)



total_p<-ggarrange(TCGA_cor_p2,TCGA_cor_p3,ncol = 1,nrow = 2)


pdf("14.Figure_plot/Figure5_S8_correlation.pdf",width = 12,height = 16)
print(total_p)
dev.off()