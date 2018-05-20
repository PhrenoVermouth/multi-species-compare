setwd("/Volumes/Work/Mac_WorkSpace/leptin/multi_spec_comp")
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(dplyr)
library(randomForest)
library(org.Hs.eg.db)
library(ggsignif)
library(ggsci)
library(WGCNA)

gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text( size = 14, face = "bold"),axis.title =element_text(size=14,face = "bold") ,axis.text=element_text(size=16))

########################################
#################### expr matrix
########################################

hum <- read.csv("Human_normat_allgene.txt",sep=' ')
mus <- read.csv("Mosue_normat.txt",sep=" ",row.names = 1)
rat <- read.csv("rat_all_normat_log.txt",row.names = 1)
maca <- read.csv("macaca_normat_log.txt",row.names = 1)
homolog <- read.csv('Hum_Mus_Rat_Macaca.csv',stringsAsFactors = F,header = T,na.strings = "") #10212
#dim(homolog)
homolog <- na.omit(homolog)
lsg <- as.vector(read.csv("liver.txt",sep="\t",header = F)$V1)
lsg <- bitr(lsg, fromType="REFSEQ" , toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
lsg <- unique(lsg$SYMBOL)

#####################add inflamme here!!!!
inflam <- as.character(read.csv("inflam.tab")$V1)

diff_hum <- as.vector(read.csv("Hum_NAFLD.csv",header = F)$V1)
diff_rat <- as.vector(read.csv("rat_diff.csv",header = F)$V1)
diff_macaca <- as.vector(read.csv("macaca_diff.csv",header = F)$V1)
#diff_rat_w16 <- as.vector(read.csv("rat_diff_w16.csv",header = F)$V1)
diff_rat <- homolog[homolog[,4] %in% diff_rat,]$Hum
#diff_rat_w16 <- homolog[homolog[,4] %in% diff_rat_w16,]$Hum
diff_macaca <- homolog[homolog[,2] %in% diff_macaca,]$Hum
com_matrix <- cbind(hum[homolog$Hum,],mus[homolog$Mouse,],rat[homolog$Rat,],maca[homolog$Macaca,])
com_matrix <- na.omit(com_matrix) #9909
com_matrix <- apply(com_matrix,2, function(x) (x-mean(x))/sd(x))

com_matrix_inflam <- com_matrix[which(rownames(com_matrix) %in% inflam),]
#From now variable 'com_matrix' are liver-specific genes by default unless re-declaraition.
com_matrix <- com_matrix[which(rownames(com_matrix) %in% lsg),] #162 #164 #196 #179
com_matrix_rat <- com_matrix[which(rownames(com_matrix) %in% diff_rat),] #15 #15 #17 #16
com_matrix_hum <- com_matrix[which(rownames(com_matrix) %in% diff_hum),] #22 #19 #24 #22
com_matrix_macaca <- com_matrix[which(rownames(com_matrix) %in% diff_macaca),] #14 #14 #18 #17

#co_diff <-  intersect(diff_macaca,diff_rat)
#com_matrix_rat16 <- com_matrix[which(rownames(com_matrix) %in% diff_rat_w16),]

########################################
#################### Pheatmap
########################################

##########1st Try
p1 <- cor(com_matrix_rat)[c(72:115),c(1:71)]
pheatmap(p1[,c(str_subset(names(sort(p1["W32_KO1",])),"advanced"),str_subset(names(sort(p1["W32_KO1",])),"mild"))],
         color = c(colorRampPalette(c("navy","white"))(30),colorRampPalette(c("white","firebrick3"))(20)),
         cluster_cols = F,cluster_rows = F, border_color = NA)
###发现mildrep10,W32_WT1 异常 去除之 row27 col34
p1 <- p1[-27,-34]

##############This is an order manually picked from excel.
row_order = read.csv("row_order_minus32WT1.csv",header = F,stringsAsFactors = F)
row_order <- row_order$V1
#colnames(a)
p1 <- p1[row_order,]
pheatmap(p1[,c(str_subset(names(sort(p1["W32_KO1",])),"advanced"),str_subset(names(sort(p1["W32_KO1",])),"mild"))],
color = c(colorRampPalette(c("navy","white"))(30),colorRampPalette(c("white","firebrick3"))(20)),
cluster_cols = F,cluster_rows = F, border_color = NA)

#p1 <- cor(com_matrix_hum)[c(72:115),c(1:71)]
#pheatmap(p1[,c(str_subset(names(sort(p1["W32_KO1",])),"advanced"),str_subset(names(sort(p1["W32_KO1",])),"mild"))],
#         color = c(colorRampPalette(c("navy","white"))(30),colorRampPalette(c("white","firebrick3"))(20)),
#         cluster_cols = F,cluster_rows = F, border_color = NA)

##############GO of 179 ls-homo genes, as requested by LP
ids <- bitr(rownames(com_matrix), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
ggo <- groupGO(gene     = ids$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

barplot(ggo, drop=TRUE, showCategory=12) +
  gran_theme+ guides(fill=FALSE)
ggsave("179genes_go_bp.png",dpi=300)
ggo <- groupGO(gene     = ids$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "MF",
               level    = 3,
               readable = TRUE)
barplot(ggo, drop=TRUE, showCategory=12) +
  gran_theme+ guides(fill=FALSE)
ggsave("179genes_go_MF.png",dpi=300)


########################################
#################### Violin plot
########################################

#########Data cleaning.
library(reshape2)
mydata <- data.frame(p1,name=rownames(p1))
mydata$spec <- rep(c("Mouse","Rat","Macaca"),c(8,29,6))
mydata$treat <- rep(c("KO","WT","KO","WT","KO","WT","KO","WT"),c(4,4,9,9,6,5,3,3))
mydata$time <- rep(c("none","early","late","none"),c(8,18,11,6))
mydata<-melt(
mydata,                                      #待转换的数据集名称
id.vars=c("name","spec","treat","time"),  #要保留的主字段
variable.name="hum",         #转换后的分类字段名称（维度）
value.name="cor"            #转换后的度量值名称
)
mydata$hum <- str_split_fixed(mydata$hum,"_",2)[,1]
mydata$comp <- factor(paste0(mydata$spec,mydata$treat))



######证明一 三种物种的相关性比例差距
anno_1 <- wilcox.test(mydata[mydata$comp  == "MouseKO", "cor"],
                      mydata[mydata$comp  == "RatKO", "cor"])$p.value
anno_2 <- wilcox.test(mydata[mydata$comp  == "MouseWT", "cor"],
                      mydata[mydata$comp  == "RatWT", "cor"])$p.value
ggplot(mydata,aes(x=spec,y=cor,fill=treat)) +
  geom_violin()+
  geom_boxplot(position=position_dodge(0.9),width=0.2)+
  gran_theme   + scale_fill_npg() +
  geom_signif(annotations = c(formatC(anno_1, digits=3),formatC(anno_2, digits=3)),
              y_position = c(0.9,0.95), xmin=c(1.775,2.225), xmax=c(2.775,3.225),tip_length = c(0.7,0.1,0.7,0.2))
ggsave('species_violin.png',dpi = 300)

######证明二 早晚期Rat与treat方式
anno_1 <- wilcox.test(mydata[mydata$time  == "early" & mydata$comp == "RatWT", "cor"],
                      mydata[mydata$time  == "late" & mydata$comp == "RatWT", "cor"])$p.value
anno_2 <- wilcox.test(mydata[mydata$time  == "early" & mydata$comp == "RatKO", "cor"],
                      mydata[mydata$time  == "late" & mydata$comp == "RatKO", "cor"])$p.value

ggplot(mydata[which(mydata$spec == "Rat"),],aes(x= time,y= cor,fill=treat)) +
  geom_violin() +
  geom_boxplot(width=0.2,position=position_dodge(0.9)) +
  gran_theme + scale_fill_npg() +geom_signif(annotations = c(formatC(anno_2, digits=3),formatC("  NS")),
                                             y_position = c(0.9,0.95), xmin=c(0.775,1.225), xmax=c(1.775,2.225),tip_length = c(0.3,0.2,0.6,0.5))
 
ggsave('time_treat_violin.png',dpi = 300)

######证明三 早晚期Rat和hum脂肪肝类型的关系
anno_1 <- wilcox.test(mydata[mydata$time  == "early" & mydata$spec == "Rat" & mydata$hum == "advanced", "cor"],
                        mydata[mydata$time  == "early" & mydata$spec == "Rat" & mydata$hum == "mild", "cor"])$p.value
anno_2 <- wilcox.test(mydata[mydata$time  == "early" & mydata$spec == "Rat" & mydata$hum == "advanced", "cor"],
                      mydata[mydata$time  == "late" & mydata$spec == "Rat" & mydata$hum == "mild", "cor"])$p.value

ggplot(mydata[which(mydata$spec == "Rat"),],aes(x= time,y= cor,fill=hum)) +
  geom_violin() +
  geom_boxplot(width=0.2,position=position_dodge(0.9)) +
  gran_theme + scale_fill_npg() +
  geom_signif(annotations = c(formatC(anno_1, digits=3),formatC(anno_2, digits=3)),
              y_position = c(0.9,0.95), xmin=c(0.775,0.775), xmax=c(1.225,1.775),tip_length = c(0,0,0,0))
ggsave('time_hum_violin.png',dpi = 300)


########################################
#################### Random forest
########################################

############Data cleaning.

#com_matrix <- cbind(hum[homolog$Hum,],mus[homolog$Mouse,],rat[homolog$Rat,],maca[homolog$Macaca,])
#com_matrix <- na.omit(com_matrix) #9909
#com_matrix <- apply(com_matrix,2, function(x) (x-mean(x))/sd(x))
#hepa-specific
#com_matrix <- com_matrix[which(rownames(com_matrix) %in% lsg),] #162 #164 #196 #179

p = cor(com_matrix)[c(72:115),c(1:71)]
mydata <- data.frame(p,name=rownames(p))[-27,-34]
mydata$spec <- rep(c("Mouse","Rat","Macaca"),c(8,29,6)) 
mydata$treat <- rep(c("KO","WT","KO","WT","KO","WT","KO","WT"),c(4,4,9,9,6,5,3,3))
mydata$time <- rep(c("none","early","late","none"),c(8,18,11,6))
dim(mydata)
mean_row <- apply(mydata[,c(1:70)],1,mean)
quantile(mean_row,c(.33,.66))
mydata$class <- rep(0,43)
for(i in 1:length(mean_row)){
if(mean_row[i] <= 0.7494){mydata$class[i] = "fail"} 
else if(mean_row[i] >= 0.7582){mydata$class[i] = "okay"} 
else{mydata$class[i] = "medium"}
}
data_forest <- rownames(mydata[which(mydata$class != 'medium'),])
################## Remove : col34:mildrep10 col98 W32WT1
com_matrix <-  com_matrix[,c(-34,-98)]
###############
data_forest <- as.data.frame(t(com_matrix)[data_forest,])
#dim(data_forest)
data_forest$class <- rep(0,29) #22
data_forest[rownames(mydata[mydata$class=="okay",]),]$class = "okay"
data_forest[rownames(mydata[mydata$class=="fail",]),]$class = "fail"
data_forest$class <- factor(data_forest$class)


forest_if <- randomForest(formula = class ~ .,ntree=1000, data = data_forest, importance = TRUE, proximity = TRUE)
plot(forest_if)
#round(importance(forest_if), 2)
varImpPlot(forest_if)

#################GO
forest_gene <- as.data.frame(round(importance(forest_if), 2))
summary(forest_gene$MeanDecreaseAccuracy)
forest_top_gene <- rownames(forest_gene[which(forest_gene$MeanDecreaseAccuracy >= 2.46),]) #top 25%
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

########################################
#################### WGCNA
########################################

#### 1)read FPKM data & condition data ###############################
gene_FPKM <- rat[homolog[which(homolog$Hum %in% rownames(com_matrix)),]$Rat,]
colnames(gene_FPKM) <- gsub("Lep","KO",colnames(gene_FPKM))

condition_table <- data.frame(row.names = colnames(gene_FPKM),
                              time = str_split_fixed(colnames(gene_FPKM),"_",2)[,1],
                              treat = str_sub(gsub("Lep","KO",str_split_fixed(colnames(gene_FPKM),"_",2)[,2]),1,2))

## Data Processing：filter FPKM<1 in all sample 
#gene_FPKM_log=log2(gene_FPKM+1)
#gene_FPKM_log = gene_FPKM_log[apply(gene_FPKM_log,1,max)>=1,]
#gene_FPKM_log = gene_FPKM_log[apply(gene_FPKM_log,1,max)-apply(gene_FPKM_log,1,min)>=1,]

#GC_FPKM_log=log2(GC_gene_FPKM+1)
#result=apply(GC_FPKM_log,1,function(x) all(x<1))
#pos=which(result=="TRUE")
#GC_filter=GC_FPKM_log[-pos,] 
 

## filter for OC/GC
#gene_FPKM=gene_FPKM[apply(OC_gene_FPKM,1,max)>=1,];
#gene_FPKM=gene_FPKM[apply(OC_gene_FPKM,1,max)-apply(OC_gene_FPKM,1,min)>=1,]


##mad（c（1,2,3,4））=median(sum(abs(i-2.5)))*1.4826=1*1.4826
WGCNA_matrix=t(gene_FPKM[order(apply(gene_FPKM, 1, mad),decreasing = T),])

datTrait=data.frame(gsm = rownames(WGCNA_matrix),
                    time = condition_table$time,
                    treat = condition_table$treat)


####  cluster the samples
sampleTree=hclust(dist(WGCNA_matrix),method = "average")
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
p1<-plot(sampleTree,main="Sample clustring to detect outliers",sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)


#### 2)set BETA value(soft-thresholding power) ##############################################
#choose a set of soft-thresholding powers
powers=c(c(1:10), seq(from = 12, to=20, by=2))
#call the network topology analysis function
sft=pickSoftThreshold(WGCNA_matrix,powerVector = powers,verbose = 5)
#plot the result
par(mfrow=c(1,2))
cex1=0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main=paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
#this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="green")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

best_beta=sft$powerEstimate

### 3)set net (OC_best_beta=3,GC_best_beta=1)
net=blockwiseModules(WGCNA_matrix,
                     power = sft$powerEstimate,
                     maxBlockSize = 6000,
                     TOMType = "unsigned",minModuleSize = 10,
                     reassignThreshold = 0, mergeCutHeight = 0.25,
                     numericLabels = TRUE, pamRespectsDendro=FALSE,
                     saveTOMs = TRUE,
                     saveTOMFileBase = "AS-green-FPKM-TOM",
                     verbose = 3)
#count model number 
table(net$colors)


### 4)model visualization
#convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

## count numbers of genes and samples
nGenes=ncol(WGCNA_matrix)
nSamples=nrow(WGCNA_matrix)

## plot cluster_tree for samples
datExpr_tree<-hclust(dist(WGCNA_matrix),method = "average")
par(mar=c(0,5,2,0))
plot(datExpr_tree,main="Sample clustering",sub="",xlab="",cex.lab=2,
     cex.axis=1,cex.main=1,cex.lab=1)

sample_colors <- numbers2colors(as.numeric(factor(datTrait$time)), #different traits
                                colors = c(heat.colors(10)[c(1,3,5,7)],topo.colors(10)[c(5,4,3,2,1)]),signed = FALSE)


par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")


### 5)model and condition relationship
design=model.matrix(~0+datTrait$gsm)
colnames(design)=levels(datTrait$gsm)
moduleColors <- labels2colors(net$colors)

# Recalculate MEs with color labels
MEs0=moduleEigengenes(WGCNA_matrix,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### 6) analysis of inseresting genes
## calculate model and gene relevant matrix
modNames=substring(names(MEs),3)
geneModuleMembership=as.data.frame(cor(WGCNA_matrix,MEs,use="p"))
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))
names(geneModuleMembership)=paste("MM",modNames,sep = "")
names(MMPvalue)=paste("p.MM",modNames,sep = "")

## calculate condition and gene relevant matrix
Luminal = as.data.frame(design[,4])
names(Luminal) = "Luminal"
geneConditionSignificance=as.data.frame(cor(WGCNA_matrix,Luminal,use="p"))
GSPvalue=as.data.frame(corPvalueStudent(as.matrix(geneConditionSignificance),nSamples))
names(geneConditionSignificance)=paste("GS.",names(Luminal),sep = "")
names(GSPvalue)=paste("P.GS.",names(Luminal),sep = "")

## connect the two matrix
module="red"
column=match(module,modNames)
moduleGenes=moduleColors==module
sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneConditionSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Luminal",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

### 7) net visualization

### 8) get genenames of module (correlation bewteen model and sample > 0.5)

relation_matrix=moduleTraitCor
colnames(relation_matrix)=rownames(WGCNA_matrix)

#module_pos=which(relation_matrix>0.9,arr.ind = TRUE)
modules=apply(relation_matrix,2,function(x) rownames(relation_matrix)[which.max(x)])

num=length(rownames(module_pos))
gene_names=colnames(WGCNA_matrix)
for (i in 1:num){
  module=rownames(moduleTraitCor)[module_pos[i,1]]
  module=substr(module,3,nchar(module))
  outfile=colnames(moduleTraitCor)[module_pos[i,2]]
  mod_gene=gene_names[moduleColors==module]
  write.table(mod_gene,file = paste0("F:/lwj-works/2018-4-8_GAO_RNAseq/lwj_GAO_tempresult_week3(4.23-4.29)/pool_sample/GC/",outfile,"_specific_gene.txt"),sep = "\t",quote = FALSE,append = TRUE,row.names = FALSE,col.names = FALSE)
}


### 9) export model
