library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(ggpubr)
library(easyGgplot2)

immune_mt<-read.table("/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/Immune_cells.txt",sep="\t",header = T)
#immune_mt<-immune_mt[,which(!(grepl("s367578",colnames(immune_mt))) )]
#immune_mt_melt<-immune_mt_melt[which(immune_mt_melt$variable!="s367578"),]


es.dif<-immune_mt
rownames(es.dif)<-es.dif$Cell.Type
es.dif<-es.dif[which(rownames(es.dif)!="T cells gamma delta"&rownames(es.dif)!="T cells CD4 naive"),-1]
es.dif1<-es.dif[,which(grepl("N",colnames(es.dif)))]
es.dif2<-es.dif[,which(grepl("T",colnames(es.dif)))]
es.dif3<-es.dif[,which(grepl("P",colnames(es.dif)))]

if (T){
  library(pacman)
  es.dif1_scale <- scale(t(es.dif1))
  es.dif1.dist <- dist(es.dif1_scale, method = "manhattan")
  es.dif1.hclust <- hclust(es.dif1.dist, method = "ward.D2")
  es.dif1<-es.dif1[,es.dif1.hclust$order]
  
  es.dif2_scale <- scale(t(es.dif2))
  es.dif2.dist <- dist(es.dif2_scale, method = "manhattan")
  es.dif2.hclust <- hclust(es.dif2.dist, method = "ward.D2")
  es.dif2<-es.dif2[,es.dif2.hclust$order]
  
  
  es.dif3_scale <- scale(t(es.dif3))
  es.dif3.dist <- dist(es.dif3_scale, method = "manhattan")
  es.dif3.hclust <- hclust(es.dif2.dist, method = "ward.D2")
  es.dif3<-es.dif3[,es.dif3.hclust$order]
  
  es.dif<-cbind(es.dif1,es.dif2)
  es.dif<-cbind(es.dif,es.dif3)
}

es.dif<-es.dif[unique(c("Neutrophils","Macrophages M2","Macrophages M0",
                        "Monocytes","Dendritic cells activated",
                        "T cells regulatory (Tregs)","Dendritic cells resting","Plasma cells","Mast cells resting","T cells CD8",
                        "Mast cells activated",
                        rownames(es.dif))),]

head(es.dif)

coldata1<-data.frame(Type=rep("GN",ncol(es.dif)))
coldata1$Type[which(grepl("P",colnames(es.dif)))]<-"GCPM"
coldata1$Type[which(grepl("T",colnames(es.dif)))]<-"GC"

rownames(coldata1)<-colnames(es.dif)

Sample_info<-data.frame(Sample_ID=colnames(es.dif),Type=rep("N",ncol(es.dif)))

pheatmap_p1<-pheatmap::pheatmap(es.dif,
                   color=colorRampPalette(c('navy','navy','white','#DC0000FF','#DC0000FF'), bias=1)(50), border_color="grey",
                   annotation_col = coldata1,
                   #annotation_colors = ann_col,
                   show_rownames = T,
                   show_colnames = F,
                   scale="row",
                   cluster_rows = T,
                   cluster_cols = F,
                   treeheight_row = 15,
                   main = "Immune cell infiltration (Cibersort)",
                  annotation_names_col = F,
                  gaps_col = c(12,24)
)
library(cowplot)
pheatmap_p1<-plot_grid(pheatmap_p1$gtable)



library(reshape2)
immune_mt_melt<-melt(immune_mt,value.name = "Per",id.vars="Cell.Type")
immune_mt_melt<-immune_mt_melt[which(immune_mt_melt$Cell.Type!="T cells gamma delta"&immune_mt_melt$Cell.Type!="T cells CD4 naive"),]


immune_mt_melt$Type<-"GN"
immune_mt_melt$Type[which(grepl("P",immune_mt_melt$variable))]<-"GCPM"
immune_mt_melt$Type[which(grepl("T",immune_mt_melt$variable))]<-"GC"

Immune_cell_list<-unique(sort(immune_mt_melt$Cell.Type))

color_list<-RColorBrewer::brewer.pal(7,"Set2")[c(1:3)]
#names(color_list)<-c("Normal", "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c("N", "T" ,"P" )
immune_mt_melt$Type<-factor(immune_mt_melt$Type,levels=c("GN", "GC" ,"GCPM" ))
immune_mt_melt$variable<-gsub("N|P|T","",immune_mt_melt$variable)

immune_mt_melt$Per<-immune_mt_melt$Per*100
cell_list<-c("Macrophages M0","Macrophages M1","Dendritic cells resting","Macrophages M2","Mast cells resting","B cells memory","Plasma cells","T cells CD8","Monocytes")
select_table<-immune_mt_melt[immune_mt_melt$Cell.Type %in% cell_list,]
#my_comparisons <- list( c("Normal", "Intestinal"), c("Normal", "Diffuse"), c("Normal", "Mixed"),c("Normal", "Metastatic") )
box_p1<-ggpaired(select_table, x="Type", y="Per", fill="Type",id = "variable",
            add="jitter",line.color = "grey", line.size = 0.5,
            palette=color_list,
            xlab="Region", 
            ylab="Relative Abundance(%)", title = "Immune cell infiltration (Cibersort)",
            legend.title=" ",show.legend = F) + 
  theme_classic()+
  scale_fill_brewer("Set3")+
  stat_compare_means( label = "p.signif",method = "wilcox.test",paired = T,
                      comparisons=list(c("GC", "GN" ),c("GCPM", "GC" ),c("GCPM", "GN" )),
                      
                      tip.length = 0.01) +#配对t检验
  theme(legend.position = 'none',
        axis.text = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  facet_wrap(~Cell.Type,scales = "free")
print(box_p1)






RNA_data<-read.table("/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/Exp_TPM.mat",sep="\t",row.names=1,header=T)
#RNA_data<-RNA_data[,which(!(grepl("s379717|s371799",colnames(RNA_data))) )]

library(GSVA)
mymatrix <- RNA_data
mymatrix <- as.matrix(mymatrix)
mymatrix[1:5,1:5]
##load the gene set 
mysymbol <- read.table(file = "/home/zhengyq/data/single_cell/27.scMouse_zyy/8.RNA_Seq/xcell_marker.txt", sep = "\t", header = T, stringsAsFactors = F)
#mysymbol<-mysymbol[,c(2,1)]
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

es.dif <- gsva(mymatrix, gs, method = "ssgsea", ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)



es.dif1<-es.dif[,which(grepl("N",colnames(es.dif)))]
es.dif2<-es.dif[,which(grepl("T",colnames(es.dif)))]
es.dif3<-es.dif[,which(grepl("P",colnames(es.dif)))]

if (T){
  library(pacman)
  es.dif1_scale <- scale(t(es.dif1))
  es.dif1.dist <- dist(es.dif1_scale, method = "manhattan")
  es.dif1.hclust <- hclust(es.dif1.dist, method = "ward.D2")
  es.dif1<-es.dif1[,es.dif1.hclust$order]
  
  es.dif2_scale <- scale(t(es.dif2))
  es.dif2.dist <- dist(es.dif2_scale, method = "manhattan")
  es.dif2.hclust <- hclust(es.dif2.dist, method = "ward.D2")
  es.dif2<-es.dif2[,es.dif2.hclust$order]
  
  
  es.dif3_scale <- scale(t(es.dif3))
  es.dif3.dist <- dist(es.dif3_scale, method = "manhattan")
  es.dif3.hclust <- hclust(es.dif2.dist, method = "ward.D2")
  es.dif3<-es.dif3[,es.dif3.hclust$order]
  
  es.dif<-cbind(es.dif1,es.dif2)
  es.dif<-cbind(es.dif,es.dif3)
}




head(es.dif)

coldata1<-data.frame(Type=rep("GN",ncol(es.dif)))
coldata1$Type[which(grepl("P",colnames(es.dif)))]<-"GCPM"
coldata1$Type[which(grepl("T",colnames(es.dif)))]<-"GC"

rownames(coldata1)<-colnames(es.dif)

Sample_info<-data.frame(Sample_ID=colnames(es.dif),Type=rep("N",ncol(es.dif)))

pheatmap_p2<-pheatmap::pheatmap(es.dif,
                                color=colorRampPalette(c('navy','navy','white','#DC0000FF','#DC0000FF'), bias=1)(50), border_color="grey",
                                annotation_col = coldata1,
                                #annotation_colors = ann_col,
                                show_rownames = T,
                                show_colnames = F,
                                scale="row",
                                cluster_rows = T,
                                cluster_cols = F,
                                treeheight_row = 15,
                                main = "Immune cell infiltration (XCell)",
                                gaps_col = c(12,24)
)
library(cowplot)
pheatmap_p2<-plot_grid(pheatmap_p2$gtable)



library(reshape2)
immune_mt<-data.frame(Cell.Type=rownames(es.dif),es.dif)
immune_mt_melt<-melt(immune_mt,value.name = "Per",id.vars="Cell.Type")
immune_mt_melt<-immune_mt_melt[which(immune_mt_melt$Cell.Type!="T cells gamma delta"&immune_mt_melt$Cell.Type!="T cells CD4 naive"),]


immune_mt_melt$Type<-"GN"
immune_mt_melt$Type[which(grepl("P",immune_mt_melt$variable))]<-"GCPM"
immune_mt_melt$Type[which(grepl("T",immune_mt_melt$variable))]<-"GC"

Immune_cell_list<-unique(sort(immune_mt_melt$Cell.Type))

color_list<-RColorBrewer::brewer.pal(7,"Set2")[c(1:3)]
#names(color_list)<-c("Normal", "Intestinal","Diffuse", "Mixed", "Metastatic")
names(color_list)<-c("N", "T" ,"P" )
immune_mt_melt$Type<-factor(immune_mt_melt$Type,levels=c("GN", "GC" ,"GCPM" ))
immune_mt_melt$variable<-gsub("N|P|T","",immune_mt_melt$variable)


cell_list<-c("Adipocytes","Endothelial cells","DC","Macrophages","Macrophages M1","Macrophages M2","Mast cells","Memory B-cells","Plasma cells","MSC","Monocytes","Preadipocytes")

select_table<-immune_mt_melt[immune_mt$Cell.Type %in% cell_list,]


box_p2<-ggpaired(select_table, x="Type", y="Per", fill="Type",id = "variable",
                 add="jitter",line.color = "grey", line.size = 0.5,
                 palette=color_list,
                 xlab="Region", 
                 ylab="ssGSVA score", title = "Immune cell infiltration (XCell)",
                 legend.title=" ",show.legend = F) + 
  theme_classic()+
  scale_fill_brewer("Set3")+
  stat_compare_means( label = "p.signif",method = "wilcox.test",paired = T,
                      comparisons=list(c("GC", "GN" ),c("GCPM", "GC" ),c("GCPM", "GN" )),
                      
                      tip.length = 0.01) +#配对t检验
  #ylim(c(0.5,1.3))+
  theme(legend.position = 'none',
        axis.text = element_text(size = 10), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = 0.5,hjust=0.5),
        strip.background = element_blank())+
  facet_wrap(~Cell.Type,scales = "free",ncol = 3)
print(box_p2)






total_p1<-ggarrange(pheatmap_p1,box_p1,nrow=1,ncol=2,
                    widths  = c(1.3,1))

total_p2<-ggarrange(pheatmap_p2,box_p2,nrow=1,ncol=2,
                    widths  = c(1.3,1))

total_p<-ggarrange(total_p1,total_p2,ncol = 1,nrow = 2,heights = c(6,9))
pdf("14.Figure_plot/Figure_1_S2_RNA_seq.pdf",width = 12,height = 15)
print(total_p)
dev.off()