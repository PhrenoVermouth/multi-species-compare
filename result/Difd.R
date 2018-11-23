setwd("/Volumes/Mac_Workplace/Mac_workspace/leptin/multi-species-compare/result")
library(readxl)
library(ggplot2)
library(reshape2)
#library(ggsci)
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text( size = 10, face = "bold"),axis.title =element_text(size=10,face = "bold") ,axis.text=element_text(size=10))

gene <- read_excel("all_species_mat.xlsx")
gene_box <- gene[c(4:8),c(2:4)]
gene_box$rat <- gene$X__1[4:8]

gene_box_tran<-melt(
  gene_box,                                      #待转换的数据集名称
  id.vars=c("rat"),  #要保留的主字段
  variable.name="Human",         #转换后的分类字段名称（维度）
  value.name="Counts"            #转换后的度量值名称
)
gene_box_tran$rat <- factor(gene_box_tran$rat,levels = c("diff_w4","diff_w8","diff_w16","diff_w32","diff_w48"))
gene_box_tran$Counts <- as.numeric(gene_box_tran$Counts)

ggplot(gene_box_tran,aes(x=Human,y=Counts,fill=rat)) +
  geom_bar(stat='identity',position=position_dodge()) + 
  gran_theme +
  scale_fill_manual(values=c("LightGray","Cyan","Blue","magenta","red"))+
  theme(legend.position="none")
ggsave("diff_genes_overlap.pdf",dpi=300,width = 3,height = 3)


pathway <- read_excel("all_species_mat_kegg.xlsx")
pathway_box <- pathway[c(4:8),c(2:4)]
pathway_box$rat <- pathway$X__1[4:8]

pathway_box_tran<-melt(
  pathway_box,                                      #待转换的数据集名称
  id.vars=c("rat"),  #要保留的主字段
  variable.name="Human",         #转换后的分类字段名称（维度）
  value.name="Counts"            #转换后的度量值名称
)
pathway_box_tran$rat <- factor(pathway_box_tran$rat,levels = c("w4","w8","w16","w32","w48"))
pathway_box_tran$Counts <- as.numeric(pathway_box_tran$Counts)

ggplot(pathway_box_tran,aes(x=Human,y=Counts,fill=rat)) +
  geom_bar(stat='identity',position=position_dodge()) + 
  gran_theme +
  scale_fill_manual(values=c("LightGray","Cyan","Blue","magenta","red"))+
  theme(legend.position="none")
ggsave("diff_pathways_overlap.pdf",dpi=300,width = 3,height = 3)


