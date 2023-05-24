########update 2023/5/24
########Analysis and visualization for singlecell RNA sequencing
library(Seurat)
library(dplyr)
library(stringr)
library(MySeuratWrappers)
library(patchwork)
library(RColorBrewer)
library(ggplot2)
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(data.table)
library(ggsci)
setwd("/path/test/") #set path

harmonyrna <- readRDS("merge_T.rds") 
test.list <- SplitObject(harmonyrna, split.by ="orig.ident")
test.list <- lapply(X = test.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = test.list,nfeatures = 1000)
test.anchors <- FindIntegrationAnchors(object.list = test.list, 
                                       anchor.features = features)
Integra.combined <- IntegrateData(anchorset = test.anchors)
DefaultAssay(Integra.combined) <- "integrated"
Integra.combined <- ScaleData(Integra.combined, verbose = FALSE)
Integra.combined <- RunPCA(Integra.combined, npcs = 50, verbose = FALSE)
Integra.combined <- RunUMAP(Integra.combined, reduction = "pca", dims = 1:50)
Integra.combined <- RunTSNE(Integra.combined, reduction = "pca", dims = 1:50)
Integra.combined <- FindNeighbors(Integra.combined, reduction = "pca", dims = 1:50)
Integra.combined <- FindClusters(Integra.combined, resolution = 1)
DimPlot(Integra.combined, reduction = "umap", label = TRUE, repel = TRUE)
ggsave(filename = "umap1cluster.pdf",width = 12,height =10,limitsize = FALSE)
DimPlot(Integra.combined, reduction = "tsne", label = TRUE, repel = TRUE)
ggsave(filename = "tsne1cluster.pdf",width = 12,height =10,limitsize = FALSE)
DefaultAssay(Integra.combined) <- "RNA"
features <- c(
  "CD79A",
  "CD3D",	"CD3G",
  "KLRD1",	
  "CPA3",	"TPSAB1",
  "CD14",	"CD68",
  "VWF",	"FLT1",
  "COL1A1",	"COL1A2",	
  "KRT19","KRT18")
DotPlot(Integra.combined, features = features) + RotatedAxis()
ggsave(filename = "dotPlotcell.pdf",width = 6,height =5,limitsize = FALSE)

VlnPlot(Integra.combined, features = features,pt.size = 0,stacked = T)
ggsave(filename = "VlnPlotcell.pdf",width = 10,height =10,limitsize = FALSE)

celltype <- read.table("cell annotion.txt",sep = '\t',header = T,check.names = F)
colnames(celltype) <- c("ClusterID","celltype")
Integra.combined@meta.data$celltype ="NA"
for(i in 1:nrow(celltype)){
  Integra.combined@meta.data[which(Integra.combined@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}
View(Integra.combined@meta.data)  
DimPlot(Integra.combined , reduction = "tsne", pt.size = 0.01,label = T,group.by = "celltype")
ggsave(filename = "tsneclustercell1.pdf",width = 6,height = 5,limitsize = FALSE)
DimPlot(Integra.combined , reduction = "umap", pt.size = 0.01,label = T,group.by = "celltype")
ggsave(filename = "umapclustercell1.pdf",width = 6,height = 5,limitsize = FALSE)
EP <- subset(Integra.combined,celltype=="Epithelial cell")
DefaultAssay(EP) <- "RNA"
View(EP@meta.data)
TCGA=read.table("Malignant_signature.txt",header = T,sep = "\t",stringsAsFactors = F) 
colnames(TCGA)=c("signature","gene")
for (i in unique(TCGA$signature)) {
  TCGA_small=TCGA%>%filter(signature==i)
  genes.for.scoring <- list(TCGA_small$gene)
  EP <- AddModuleScore(object = EP, features = genes.for.scoring, name = i) 
}
View(EP@meta.data)
EP@meta.data$CB=rownames(EP@meta.data)
metadf=EP@meta.data[,c("malig","no.malig","CB")] #Initial score
#kmeans
kmeans.result <- kmeans(metadf[,1:2],2)
kmeans_df <- data.frame(
  kmeans_class=kmeans.result$cluster
)
kmeans_df$CB=rownames(kmeans_df)
metadf=merge(metadf,kmeans_df,by="CB")
small1_metadf=metadf%>%filter(kmeans_class==1)
small2_metadf=metadf%>%filter(kmeans_class==2)
if( (mean(small1_metadf$malig1) < mean(small2_metadf$malig1)) & (mean(small1_metadf$no.malig1)>mean(small2_metadf$no.malig1)) ) {
  print("1:Normal cell; 2:Malignant cell")
  metadf$classification=ifelse(metadf$kmeans_class==1,"Normal cell","Malignant cell")
}else if (mean(small1_metadf$malig1) > mean(small2_metadf$malig1) & mean(small1_metadf$no.malig1) < mean(small2_metadf$no.malig1)){
  print("1:Malignant cell; 2:Normal cell")
  metadf$celltypeNT=ifelse(metadf$kmeans_class==1,"Malignant cell","Normal cell")
}else{
  print("error")
}
metadf=metadf[,c("CB","celltypeNT")]
EP@meta.data=EP@meta.data%>%inner_join(metadf,by = "CB")
rownames(EP@meta.data)=EP@meta.data$CB
View(EP@meta.data)
table(EP@meta.data$celltypeNT,EP@meta.data$celltypeNT)
table(Integra.combined@meta.data$celltype)
Integra.combined@meta.data$celltypeNT <- Integra.combined@meta.data$celltype
for (i in row.names(EP@meta.data)){
  Integra.combined@meta.data[i,"celltypeNT"] <- EP@meta.data[i,"celltypeNT"]
}
table(Integra.combined@meta.data$celltypeNT)
Integra.combined@meta.data$celltypeNT <- factor(Integra.combined@meta.data$celltypeNT,
levels = c("B cell","T cell","NK cell","Mast cell","Macrophage","Endothelial cell","CAF","Malignant cell","Normal cell"))

###Samples proportion
trait="celltypeNT"                 
score=Integra.combined@meta.data
sameSample=row.names(score)
rt=score[sameSample,]
bioCol=c("#106CA9",	"#B33B2A","#21804C","#D8822B",
         "#282424","#DDCD7E","#79508C","#B1A5C8","#FA5E5C",
         "#E00F0A","#3D6197","#C03D64","#576150"
)
bioCol=bioCol[1:length(unique(rt[,trait]))]
rt1=rt[,c(trait, "orig.ident")]
colnames(rt1)=c("trait", "orig.ident")
df=as.data.frame(table(rt1))
df=ddply(df, .(orig.ident), transform, percent = Freq/sum(Freq) * 100)
df=ddply(df, .(orig.ident), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
p=ggplot(df, aes(x = factor(orig.ident), y = percent, fill = trait)) +ggtitle("")+
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+theme(axis.text.x=element_text(size=12,  color = "black",angle = 90, vjust =0.5))+
  xlab("")+ ylab("Percent weight")+  guides(fill=guide_legend(title=trait))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5) 

pdf(file="proportion.pdf")
print(p)
dev.off()

###MKI67 proportion
Integra.combined@meta.data$MKI67 <- Integra.combined@assays$RNA@data["MKI67",]
Integra.combined@meta.data[Integra.combined@meta.data$MKI67>0,"MKI67"] <- "MKI67+"
Integra.combined@meta.data[Integra.combined@meta.data$MKI67==0,"MKI67"] <- "MKI67-"
Integra.combined@meta.data$celltype_MKI67 <- paste0(Integra.combined@meta.data$MKI67,Integra.combined@meta.data$celltypeNT)
table(Integra.combined@meta.data$celltype_MKI67)
tab <- table(Integra.combined@meta.data$orig.ident,Integra.combined@meta.data$celltype_MKI67)
tab <- as.matrix(tab)
for(i in 1:nrow(tab)) {
  tab[i,] <- tab[i,]/sum(tab[i,])
}
tab <- as.data.frame(tab)
colnames(tab) <- c("Sample","Cell","Fre")
tab2 <- tab[-(grep("MKI67-",tab$Cell)),]
subtype <- cbind(orig.ident2=Integra.combined@meta.data$orig.ident,subtype=Integra.combined@meta.data$subtype)
head(subtype)
subtype <- as.data.frame(subtype)
subtype <- subtype[!duplicated(subtype$orig.ident),]
for(i in 1:nrow(subtype)){
  tab2[which(tab2$Sample == subtype$orig.ident[i]),'Subtype'] <- subtype$subtype[i]
}
tab3 <- aggregate(Fre ~ Cell + Subtype, data = tab2, FUN = mean)
tab3$Percent <- round(tab3$Fre*100,1)
tab3$label=paste0(tab3$Percent, "%")

bioCol=c("#106CA9",	"#B33B2A","#21804C","#D8822B",
         "#DDCD7E","#79508C","#B1A5C8","#FA5E5C",
         "#E00F0A","#3D6197","#6D6466","#282424"
)
p=ggplot(tab3, aes(x = factor(Subtype), y = Percent, fill = Cell)) +#ggtitle("")+
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+theme(axis.title.x=element_text(angle=90,size=20,vjust=2))+
  xlab("Subtype")+ ylab("Percent weight")+  guides(fill=guide_legend(title="Cell"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  
  theme_bw() + theme(  panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       legend.position = "top",
                       panel.background = element_blank())
pdf(file="MKI67proportion.pdf")
print(p)
dev.off()

###CAR-T targets of CAF

Tumor <- readRDS("Merge_T.rds")
Normal <- readRDS("Merge_N.rds")
Tumor@meta.data$celltype <- gsub("CAF","CAFT",Tumor@meta.data$celltype)
Normal@meta.data$celltype <- gsub("CAF","CAFN",Normal@meta.data$celltype)
CAFT <- subset(Tumor,subset = (celltype=="CAFT"))
NT <- merge(CAFT,Normal)

gene <- read.table("gene_logFC1up.txt",sep = '\t')
mebran <- NT[gene$V1,]
future::plan("multiprocess", workers = 10)
options(future.globals.maxSize= 16000 * 1024^2)
CAFTmarker <- FindMarkers(mebran,ident.1="CAFT",logfc.threshold = 0, min.pct = 0.6,only.pos = T)
CAFTmarker.final <- CAFTmarker %>% filter(p_val_adj <0.05 &pct.2 < 0.1) 
nrow(CAFTmarker)
write.table(CAFTmarker.final,'CAFTmarker.final.txt',sep = '\t',col.names = T,row.names = T,quote = F)

##AND
data <- mebran@assays$RNA@counts 
data <- as.matrix(data)
data <- ifelse(data>0,1,data)
a0 <- c()
for(u in 1:nrow(data)){
  for(d in 1:nrow(data)){
    a <- paste0(rownames(data)[u],"_",rownames(data)[d])
    a0 <- c(a,a0)
  }}
gene_pairs <- a0

new_gene_matrix <- matrix(0, nrow = length(gene_pairs), ncol = ncol(data), dimnames = list(gene_pairs, colnames(data)))
for (i in 1:length(gene_pairs)) {
  gene_pair <- strsplit(gene_pairs[i], "_")[[1]]
  new_gene_matrix[i, ] <- data[gene_pair[1], ] + data[gene_pair[2], ]
}


new_gene_matrix2 <- ifelse(new_gene_matrix<2,0,new_gene_matrix)
pair <- CreateSeuratObject(new_gene_matrix2)
pair@meta.data$celltype <- mebran@meta.data$celltype
table(pair@meta.data$celltype)
DefaultAssay(pair) <- "RNA"
Idents(pair) <- pair@meta.data$celltype
CAFTANDmarker <- FindMarkers(pair,ident.1="CAFT",logfc.threshold = 0, min.pct = 0.6,only.pos = T)
CAFTANDmarker.final <-  CAFTANDmarker %>% filter(p_val_adj <0.05 &pct.2 < 0.1) 
write.table(CAFTANDmarker.final,'CAFTANDmarker.final.txt',sep = '\t',col.names = T,row.names = T,quote = F)

##OR

new_gene_matrix3 <- ifelse(new_gene_matrix<1,0,new_gene_matrix)
pair <- CreateSeuratObject(new_gene_matrix3)
pair@meta.data$celltype <- mebran@meta.data$celltype
DefaultAssay(pair) <- "RNA"
Idents(pair) <- pair@meta.data$celltype
CAFTORmarker <- FindMarkers(pair,ident.1="CAFT",logfc.threshold = 0, min.pct = 0.6,only.pos = T)
CAFTORmarker.final <-  CAFTORmarker %>% filter(p_val_adj <0.05 &pct.2 < 0.1) 
View(CAFTORmarker.final)

write.table(CAFTORmarker.final,'CAFTORmarker.final.txt',sep = '\t',col.names = T,row.names = T,quote = F)

##NOT
gene_up <- read.table("gene_logFC1up.txt",sep = '\t')
mebran_up <- NT[gene_up$V1,]
data_up <- mebran_up@assays$RNA@counts 
data_up <- as.matrix(data_up)
data_up <- ifelse(data_up>0,1,data_up)

gene_down <- read.table("gene_logFC1down.txt",sep = '\t')
mebran_down <- NT[gene_down$V1,]
data_down <- mebran_down@assays$RNA@counts 
data_down <- as.matrix(data_down)
data_down <- ifelse(data_down>0,1,data_down)
a0 <- c()
for(u in 1:nrow(data_up)){
  for(d in 1:nrow(data_down)){
    a <- paste0(rownames(data_up)[u],"_",rownames(data_down)[d])
    a0 <- c(a,a0)
  }}
gene_pairs <- a0


new_gene_matrix <- matrix(0, nrow = length(gene_pairs), ncol = ncol(data_up), dimnames = list(gene_pairs, colnames(data_up)))
for (i in 1:length(gene_pairs)) {
  gene_pair <- strsplit(gene_pairs[i], "_")[[1]]
  new_gene_matrix[i, ] <- data_up[gene_pair[1], ] - data_down[gene_pair[2], ]
}


new_gene_matrix <- ifelse(new_gene_matrix<1,0,new_gene_matrix)
pair <- CreateSeuratObject(new_gene_matrix)
pair@meta.data$celltype <- mebran_up@meta.data$celltype
DefaultAssay(pair) <- "RNA"
Idents(pair) <- pair@meta.data$celltype
CAFTNOTmarker <- FindMarkers(pair,ident.1="CAFT",logfc.threshold = 0, min.pct = 0.6,only.pos = T)
CAFTNOTmarker.final <-  CAFTNOTmarker %>% filter(p_val_adj <0.05 &pct.2 < 0.1) 
nrow(CAFTNOTmarker.final)

write.table(CAFTNOTmarker.final,'CAFTNOTmarker.final.txt',sep = '\t',col.names = T,row.names = T,quote = F)

###cell-cell communication

data.input <- GetAssayData(Integra.combined, assay = "RNA", slot = "data")
identity <- subset(Integra.combined@meta.data, select = "celltypeingtegra_cafNNMT")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "celltypeingtegra_cafNNMT")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor", "Cell-Cell Contact" ))
cellchat@DB <- CellChatDB.use 
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 10)
options(future.globals.maxSize= 16000 * 1024^2)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
pathways.show <- c("CXCL" ,"VEGF","MIF","LAMININ")        #Specific pathway
pairLR  <- extractEnrichedLR(cellchat, signaling = pathways.show , geneLR.return = FALSE)
pairLR2 <- as.data.frame(pairLR[c(5,7:13,47),])
colnames(pairLR2)[1]<- "interaction_name"
pairLR2
netVisual_bubble(cellchat, sources.use = c(6,7), targets.use =c(1:5,8) , remove.isolate = FALSE,pairLR.use =pairLR2)