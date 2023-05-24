###########update 2023/5/24
####Louvain clustering
library(Seurat)
library(ggplot2)
library(patchwork)
library(fpc)
library (cluster)
library (vegan)
library(pheatmap)
library(corrplot) 
setwd("/path/test")
data = read.table(file = "1ibaq_rm0log_genesymbol_T_50_knn_scale_immune.txt",header = T,row.names = 1,sep = "\t",check.names = F)
n <- ncol(data)
data2 = as.sparse(data)
data2 = CreateSeuratObject(data2)
#data2 <- NormalizeData(data2)
data2 <- FindVariableFeatures(data2,selection.method = "disp")
data2 <- ScaleData(data2)
data2 <- RunPCA(data2, verbose = FALSE)
data3 <- FindNeighbors(data2, dims = 1:15)
data4 <- FindClusters(data3, resolution =1.15, verbose = FALSE)
data4@active.ident
tmp2 = data4@active.ident

tmp3 <- as.data.frame(tmp2)
colnames(tmp3) <- "group"

tmp3$group <- as.numeric(tmp3$group)
data5 <- as.data.frame(t(data))
group <- tmp3

stats <- cluster.stats(dist(data5), as.vector(group[,1]))
Calinski_Harabasz_index <- stats$avg.silwidth

rt=data      #input data
rt=t(rt)      #turn data
rt <- merge(group,rt,by = "row.names", all = T)
rownames(rt) <- rt[,1]
rt <- rt[,-1]
rt <- rt[order(rt$group),]
rt <- rt[,-1]
rt=t(rt)  

M=cor(rt, method = "spearman")   #relationship
M <- as.data.frame(M)
a <- min(M)
b <- max(M)

bk <- c(seq(a,0-(0-a)/40,by=(0-a)/40),seq(0,0-a,by=(0-a)/40))
n=ncol(M)
sameSample=colnames(M)

immData=cbind(M, NewSubtype=group[sameSample,])
immData <- as.data.frame(immData)

data=immData
data=data[order(data$NewSubtype),]
Type=data[,((n+1):ncol(data))]
Type=as.data.frame(Type)
rownames(Type)=rownames(data)
colnames(Type)="NewSubtype"
data=t(data[,1:n])
data=data[colnames(data),]
bioCol=c("#094B8D","#FF9900","#FF0000","#6E568C","#808080","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()

ICIcol=bioCol[1:length(levels(factor(Type$NewSubtype)))]
names(ICIcol)=levels(factor(Type$NewSubtype))
ann_colors[["NewSubtype"]]=ICIcol

#next plot

pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color=c(colorRampPalette(colors=c("#094B8D","white"))(length(bk)/2),
                 colorRampPalette(color=c("white","red"))(length(bk)/2)),
         cluster_cols =F,
         cluster_rows =F,
         show_colnames=F,
         show_rownames=F,
         border=FALSE,
         gaps_col = cumsum(table(Type$NewSubtype)),
         gaps_row =cumsum(table(Type$NewSubtype)),
         legend_breaks=seq(-0.4,0.4,0.2),breaks=bk)