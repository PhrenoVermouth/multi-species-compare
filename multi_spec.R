setwd("/Volumes/Work/Mac_WorkSpace/leptin/multi_spec_comp")
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(dplyr)
library(randomForest)
library(org.Hs.eg.db)
library(org.Rn.eg.db)
library(ggsignif)
library(ggsci)
library(WGCNA)
library(reshape2)
library(cowplot)

gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text( size = 14, face = "bold"),axis.title =element_text(size=24,face = "bold") ,axis.text=element_text(size=16))

########################################
#################### expr matrix
########################################

hum <- read.csv("Human_normat_allgene_2.csv",row.names = 1)
hum <- hum[,c(grep("^c",colnames(hum)),grep("^h",colnames(hum)),grep("^s",colnames(hum)),grep("^n",colnames(hum)))]
mus <- read.csv("Mosue_normat.txt",sep=" ",row.names = 1)
rat <- read.csv("rat_all_normat_log.txt",row.names = 1)
maca <- read.csv("macaca_normat_log.txt",row.names = 1)
homolog <- read.csv('Hum_Mus_Rat_Macaca.csv',stringsAsFactors = F,header = T,na.strings = "") #10212
#dim(homolog)
homolog <- na.omit(homolog)
lsg <- as.vector(read.csv("liver.txt",sep="\t",header = F)$V1)
lsg <- bitr(lsg, fromType="REFSEQ" , toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
lsg <- unique(lsg$SYMBOL)
inflam <- as.character(read.csv("inflam.tab")$V1)

nash_up <- homolog[which(homolog$Hum %in% as.character(read.csv("Hum_nash_005_up.csv")$x)),]$Rat
nash_down <- homolog[which(homolog$Hum %in% as.character(read.csv("Hum_nash_005_down.csv")$x)),]$Rat
steatosis_up <- homolog[which(homolog$Hum %in% as.character(read.csv("Hum_steatosis_005_up.csv")$x)),]$Rat
steatosis_down <- homolog[which(homolog$Hum %in% as.character(read.csv("Hum_steatosis_005_down.csv")$x)),]$Rat


# nash <- as.character(read.csv("Hum_nash_005.csv")$x)
# steatosis <- as.character(read.csv("Hum_steatosis_005.csv")$x)
# nash  <- homolog[which(homolog$Hum %in% nash),]$Rat
# steatosis  <- homolog[which(homolog$Hum %in% steatosis),]$Rat



diff_hum <- as.vector(read.csv("Hum_NAFLD.csv",header = F)$V1)
diff_rat_up <- as.vector(read.csv("rat_diff_up.csv",header = F)$V1)
diff_rat_down <- as.vector(read.csv("rat_diff_down.csv",header = F)$V1)
diff_macaca <- as.vector(read.csv("macaca_diff.csv",header = F)$V1)
#diff_rat_w16 <- as.vector(read.csv("rat_diff_w16.csv",header = F)$V1)
diff_rat_up <- homolog[homolog[,4] %in% diff_rat_up,]$Hum
diff_rat_down <- homolog[homolog[,4] %in% diff_rat_down,]$Hum
diff_rat <- unique(c(diff_rat_down,diff_rat_up))
#diff_rat_w16 <- homolog[homolog[,4] %in% diff_rat_w16,]$Hum
diff_macaca <- homolog[homolog[,2] %in% diff_macaca,]$Hum




com_matrix <- cbind(hum[homolog$Hum,],mus[homolog$Mouse,],rat[homolog$Rat,],maca[homolog$Macaca,])
com_matrix <- na.omit(com_matrix) #9909
com_matrix <- apply(com_matrix,2, function(x) (x-mean(x))/sd(x))

com_matrix_inflam <- com_matrix[which(rownames(com_matrix) %in% inflam),]
#From now on variable 'com_matrix' are liver-specific genes by default unless re-declaraition.
com_matrix <- com_matrix[which(rownames(com_matrix) %in% lsg),] #162 #164 #196 #179
com_matrix_rat <- com_matrix[which(rownames(com_matrix) %in% diff_rat),] #15 #15 #17 #16
com_matrix_hum <- com_matrix[which(rownames(com_matrix) %in% diff_hum),] #22 #19 #24 #22
com_matrix_macaca <- com_matrix[which(rownames(com_matrix) %in% diff_macaca),] #14 #14 #18 #17

#co_diff <-  intersect(diff_macaca,diff_rat)
#com_matrix_rat16 <- com_matrix[which(rownames(com_matrix) %in% diff_rat_w16),]

########################################
#################### Pheatmap 验证物种差异性
########################################

##########inflam
p1 <- cor(com_matrix_inflam)[c(74:116),c(1:73)]

p1 <- p1[,-21] 
p1 <- p1[,-3] #Remove control.11 healthyobese15

###########lsg
p1 <- cor(com_matrix)[c(74:116),c(1:73)]

###########lsg+rat
p1 <- cor(com_matrix_rat)[c(74:116),c(1:73)]
p1 <- p1[,-55] 
p1 <- p1[,-3]#remove con_rep11 ste_rep9


##############This is an order manually picked from excel.
row_order = read.csv("row_order.csv",header = F,stringsAsFactors = F)
row_order <- row_order$V1

#colnames(a)
p1 <- p1[row_order,]
pheatmap(p1,
color = c(colorRampPalette(c("navy","white"))(30),colorRampPalette(c("white","firebrick3"))(20)),
cluster_cols = F,cluster_rows = F, border_color = NA)

p1 <- data.frame(p1)
p1$species <- rep(c("Mouse","Rat","Macaca"),c(8,29,6))

mydata<-melt(
  p1,                                      #待转换的数据集名称
  id.vars=c("species"),  #要保留的主字段
  variable.name="hum_class",         #转换后的分类字段名称（维度）
  value.name="cor"            #转换后的度量值名称
)
mydata$species <- factor(mydata$species,levels=c("Macaca","Rat","Mouse"))

anno_1 <- wilcox.test(mydata[mydata$species  == "Macaca", "cor"],
                      mydata[mydata$species  == "Rat", "cor"])$p.value
anno_2 <- wilcox.test(mydata[mydata$species  == "Rat", "cor"],
                      mydata[mydata$species  == "Mouse", "cor"])$p.value
ggplot(mydata,aes(x=species,y=cor,fill=species)) +
  geom_boxplot(outlier.color = "white")+
   geom_signif(annotations = c(formatC(anno_1, digits=3),formatC(anno_2, digits=3)),
        y_position = c(0.93,0.9), xmin=c(1,2), xmax=c(2,3),tip_length = c(0,0,0,0)) +
  gran_theme   + scale_fill_npg() 
ggsave('species_boxplot.pdf',dpi = 300)



########################################
############# Pheatmap 验证炎症与脂肪基因
########################################

############ 1.Both up-regu in rat & hum

diff_rat_up <- as.vector(read.csv("rat_diff_up.csv",header = F)$V1)
both_up <- intersect(diff_rat_up,c(nash_up,steatosis_up))
p <- na.omit(rat[both_up,][,c(colnames(rat)[grep("WT.$",colnames(rat))],
                           colnames(rat)[grep("KO.$",colnames(rat))])])
p <- data.frame(t(apply(p,1,function(x) (x-mean(x))/sd(x))))


############# 2. Both down-regu in rat & hum

diff_rat_down <- as.vector(read.csv("rat_diff_down.csv",header = F)$V1)
both_down <- intersect(diff_rat_down,c(nash_down,steatosis_down))
p <- na.omit(rat[both_down,][,c(colnames(rat)[grep("WT.$",colnames(rat))],
            colnames(rat)[grep("KO.$",colnames(rat))])])
p <- data.frame(t(apply(p,1,function(x) (x-mean(x))/sd(x))))


#############3. 1 plus 2
p <- na.omit(rat[c(both_up,both_down),][,c(colnames(rat)[grep("WT.$",colnames(rat))],
                                colnames(rat)[grep("KO.$",colnames(rat))])])
p <- data.frame(t(apply(p,1,function(x) (x-mean(x))/sd(x))))


both_heatmap <- pheatmap(p[,c(1:29)],cluster_cols = F,cluster_rows = F,border_color = NA, 
          color = c(colorRampPalette(c("navy","white","firebrick3"))(100)),
         breaks = seq(-4,4,length.out = 100),legend = F)
ggsave("bothup_down.pdf",both_heatmap,width = 4,height = 8)

######################K-means贴标签 必然分两类
cluster_p <- kmeans(p,2)
p$label <- as.character(cluster_p$cluster)
p <- p[order(p$label),]
p$gene <- rownames(p)
mydata<-melt(
  p,                                      #待转换的数据集名称
  id.vars=c("label","gene"),  #要保留的主字段
  variable.name="name",         #转换后的分类字段名称（维度）
  value.name="expr"            #转换后的度量值名称
)
mydata$time <- factor(str_split_fixed(mydata$name,"_",2)[,1],levels = c("W4","W8","W16","W32","W48"))
mydata$treat <- str_sub(str_split_fixed(mydata$name,"_",2)[,2],1,2)

both_boxplot <- ggplot(mydata,aes(x=time,y=expr,fill=treat)) +
  geom_boxplot()+
  gran_theme   + scale_fill_npg() + facet_grid(label~ .)

plot_grid(both_heatmap, both_boxplot,  align = "h")
ggsave("heatmap_boxplot.pdf",dpi=300,width = 12, height = 7.5)


#########################################################
############# 看先验的geneset能否对应物种相关性的时间关系
#########################################################

inflam_geneset <- as.character(read.csv("inflam_geneset_LP.csv",col.names = "V1")$V1)
fat_geneset <- as.character(read.csv("fat_geneset_LP.csv",col.names = "V1")$V1)
############Re-define com_matrix without liver-specific genes 
com_matrix <- cbind(hum[homolog$Hum,],mus[homolog$Mouse,],rat[homolog$Rat,],maca[homolog$Macaca,])
com_matrix <- na.omit(com_matrix) 
com_matrix <- apply(com_matrix,2, function(x) (x-mean(x))/sd(x))
############Heatmap of correlation
com_matrix_inf_lp <- com_matrix[which(rownames(com_matrix) %in% toupper(inflam_geneset)),]
p1 <- data.frame(cor(com_matrix)[c(74:116),c(1:73)])
p1 <- p1[,-55] 
p1 <- p1[,-3]
pheatmap(p1[row_order,],cluster_cols = F,cluster_rows = F, border_color = NA)

com_matrix_fat_lp <- com_matrix[which(rownames(com_matrix) %in% toupper(fat_geneset)),]
p1 <- cor(com_matrix_fat_lp)[c(74:116),c(1:73)]
p1 <- p1[,-55] 
p1 <- p1[,-3]

pheatmap(p1[row_order,],cluster_cols = F,cluster_rows = F, border_color = NA)

com_matrix_all_lp <- com_matrix[which(rownames(com_matrix) %in% toupper(c(inflam_geneset,fat_geneset))),]
p1 <- cor(com_matrix_all_lp)[c(74:116),c(1:73)]
p1 <- p1[,-55] 
p1 <- p1[,-3]
pheatmap(p1[row_order,],cluster_cols = F,cluster_rows = F, border_color = NA,color = c(colorRampPalette(c("navy","white","firebrick3"))(100)))

######3p1 <- apply(p1,2,function(x) (x-mean(x))/sd(x))
p1$sample <- rownames(p1)
mydata<-melt(
  p1,                                      #待转换的数据集名称
  id.vars="sample",  #要保留的主字段
  variable.name="hum_name",         #转换后的分类字段名称（维度）
  value.name="cor"            #转换后的度量值名称
)
mydata <- mydata[grep("^W\\d",mydata$sample),]
mydata$time <- factor(str_split_fixed(mydata$sample,"_",2)[,1],levels = c("W4","W8","W16","W32","W48"))
mydata$treat <- str_sub(str_split_fixed(mydata$sample,"_",2)[,2],1,-2)
mydata$human_class <- str_split_fixed(mydata$hum_name,"_",2)[,1]
ggplot(mydata,aes(x=human_class,y=cor,fill=treat)) +
  geom_boxplot()+
  gran_theme   + scale_fill_npg() + facet_grid(time~ .)


ids <- bitr(forest_top_gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
ggo <- groupGO(gene     = ids$ENTREZID,
OrgDb    = org.Hs.eg.db,
ont      = "BP",
level    = 3,
readable = TRUE)

barplot(ggo, drop=TRUE, showCategory=12) +
gran_theme+ guides(fill=FALSE)
ggsave("forest_go.png",dpi=300)
ggo <- groupGO(gene     = ids$ENTREZID,
OrgDb    = org.Hs.eg.db,
ont      = "MF",
level    = 3,
readable = TRUE)
barplot(ggo, drop=TRUE, showCategory=12) +
gran_theme+ guides(fill=FALSE)
ggsave("forest_go_MF.png",dpi=300)
ggo <- groupGO(gene     = ids$ENTREZID,
OrgDb    = org.Hs.eg.db,
ont      = "CC",
level    = 3,
readable = TRUE)
barplot(ggo, drop=TRUE, showCategory=12) +
gran_theme+ guides(fill=FALSE) #Remove the legend
ggo <- groupGO(gene     = ids$ENTREZID,
OrgDb    = org.Hs.eg.db,
ont      = "BP",
level    = 4,
readable = TRUE)

