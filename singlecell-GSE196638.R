rm(list = ls())

library(remotes)
library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(SeuratData)
library(viridis)
library(harmony)
library(SingleR)
library(celldex)

###Set the working path
setwd("E:/Rflie/GSE196638")
####reading data
pbmc.data1 <- Read10X(data.dir = "NL1") 
pbmc.data2 <- Read10X(data.dir = "NL2") 
pbmc.data3 <- Read10X(data.dir = "NL3")
pbmc.data4 <- Read10X(data.dir = "TXP24") 
pbmc.data5 <- Read10X(data.dir = "TXP25") 
pbmc.data6 <- Read10X(data.dir = "TXP26") 

###create a seurat object
Normal1 <- CreateSeuratObject(pbmc.data1, project = "Normal", min.cells = 3, min.features = 200)
Normal2 <- CreateSeuratObject(pbmc.data2, project = "Normal", min.cells = 3, min.features = 200)
Normal3 <- CreateSeuratObject(pbmc.data3, project = "Normal", min.cells = 3, min.features = 200)
Emphysema1 <- CreateSeuratObject(pbmc.data4, project = "Emphysema", min.cells = 3, min.features = 200)
Emphysema2 <- CreateSeuratObject(pbmc.data5, project = "Emphysema", min.cells = 3, min.features = 200)
Emphysema3 <- CreateSeuratObject(pbmc.data6, project = "Emphysema", min.cells = 3, min.features = 200)

####线粒体基因
Normal1[['percent_mito']] <- PercentageFeatureSet(Normal1, pattern = "^MT-")
Normal2[['percent_mito']] <- PercentageFeatureSet(Normal2, pattern = "^MT-")
Normal3[['percent_mito']] <- PercentageFeatureSet(Normal3, pattern = "^MT-")
Emphysema1[['percent_mito']] <- PercentageFeatureSet(Emphysema1, pattern = "^MT-")
Emphysema2[['percent_mito']] <- PercentageFeatureSet(Emphysema2, pattern = "^MT-")
Emphysema3[['percent_mito']] <- PercentageFeatureSet(Emphysema3, pattern = "^MT-")

{#####
  VlnPlot(Normal1, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Normal2, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Normal3, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Emphysema1, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Emphysema2, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Emphysema3, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  
  ####nFeature_RNA,nCount_RNA
  ngenemax1 <- unname(quantile(Normal1$nFeature_RNA,0.90))
  umimax1 <- unname(quantile(Normal1$nCount_RNA,0.90))
  ngenemax2 <- unname(quantile(Normal2$nFeature_RNA,0.90))
  umimax2 <- unname(quantile(Normal2$nCount_RNA,0.90))
  ngenemax3 <- unname(quantile(Normal3$nFeature_RNA,0.90))
  umimax3 <- unname(quantile(Normal3$nCount_RNA,0.90))
  ngenemax4 <- unname(quantile(Emphysema1$nFeature_RNA,0.90))
  umimax4 <- unname(quantile(Emphysema1$nCount_RNA,0.90))
  ngenemax5 <- unname(quantile(Emphysema2$nFeature_RNA,0.90))
  umimax5 <- unname(quantile(Emphysema2$nCount_RNA,0.90))
  ngenemax6 <- unname(quantile(Emphysema3$nFeature_RNA,0.90))
  umimax6 <- unname(quantile(Emphysema3$nCount_RNA,0.90))
  
  ###质控筛选
  Normal1 <- subset(x = Normal1, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax1 & nCount_RNA > 0 & nCount_RNA  < umimax1 & percent_mito < 10)
  Normal2 <- subset(x = Normal2, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax2 & nCount_RNA > 0 & nCount_RNA  < umimax2 & percent_mito < 10)
  Normal3 <- subset(x = Normal3, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax3 & nCount_RNA > 0 & nCount_RNA  < umimax3 & percent_mito < 10)
  Emphysema1 <- subset(x = Emphysema1, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax4 & nCount_RNA > 0 & nCount_RNA  < umimax4 & percent_mito < 10)
  Emphysema2 <- subset(x = Emphysema2, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax5 & nCount_RNA > 0 & nCount_RNA  < umimax5 & percent_mito < 10)
  Emphysema3 <- subset(x = Emphysema3, subset = nFeature_RNA > 200 & nFeature_RNA < ngenemax6 & nCount_RNA > 0 & nCount_RNA  < umimax6 & percent_mito < 10)
  
  #####
  VlnPlot(Normal1, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Normal2, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Normal3, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Emphysema1, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Emphysema2, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
  VlnPlot(Emphysema3, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
}


# MergeData
Human <- merge(x = Normal1, y = c(Normal2,Normal3,Emphysema1,Emphysema2,Emphysema3))
table(Human@meta.data$orig.ident)

all.genes=rownames(Human)
Human <- ScaleData(Human, features = all.genes)
###Data were normalized, hypervariable genes were searched, homogenized, and PCA dimensionality reduction was performed

Human <- NormalizeData(Human, normalization.method = "LogNormalize", scale.factor = 10000) 
Human <- FindVariableFeatures(Human, selection.method = "vst", nfeatures = 2000)
Human <- RunPCA(Human, features = VariableFeatures(object = Human),verbose = FALSE)###PCA分析

#####Use harmony to remove batch effects
Human<-RunHarmony(Human,group.by.vars = "orig.ident")
#####harmony法是根据PCA后的结果进行批次效应的去除，单线程跑，省内存
#####group.by.vars = "orig.ident"，是表示两个样本间去批次效应

#####使用harmony后的结果挑选合适的PC数聚类
ElbowPlot(Human, ndims = 50, reduction = "harmony")####根据散点图和CD文献选择拐点PC数16
Human <- FindNeighbors(Human, reduction= "harmony",dim=1:27) %>% FindClusters(resolution = 0.4)##dim=1:16拐点数##0.5分成12簇  1.2分成16
####一定要指定harmony(reduction="harmony"),要不然harmony白做
table(Human@meta.data$seurat_clusters)##查看每一类有多少个细胞，可以看到一共有17个簇


####UMAP的绘制
Human <- RunUMAP(Human,reduction= "harmony", dim=1:20)
plot1= DimPlot(Human,reduction = "umap", label=T)
plot2= DimPlot(Human,reduction = "umap", group.by = "orig.ident")

plot1 + plot2

# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
Human.markers <- FindAllMarkers(Human, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)


##write.csv(Human.markers ,"markers2.csv")

celltype <- c("T cell","Monocyte","NK T cell","Fibroblast","CD8+ T cell","Macrophage","NK cell",
              "Endothelial cell","B cell","Mast cell","Epithelial cell","Endothelial cell",
              "Smooth muscle cell","Endothelial cell","Monocyte","Macrophage","Endothelial cell",
              "Plasma cell","Fibroblast","Epithelial cell","Plasmacytoid dendritic cell(pDC)",
              "Plasma cell","Epithelial cell")

names(celltype ) <- levels(Human)
Human <- RenameIdents(Human, celltype)
P1 <- DimPlot(Human, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
P2 <- DimPlot(Human,reduction = "umap", group.by = "orig.ident",pt.size = 0.5)

P1+P2
DimPlot(Human,split.by = "orig.ident",label = T)+ NoLegend()

table(Idents(Human), Human$orig.ident)#各组不同细胞群细胞数

###组间细胞差异的展示
table(Human$orig.ident)
prop.table(table(Idents(Human)))
table(Idents(Human), Human$orig.ident)#各组不同细胞群细胞数

Cellratio <- prop.table(table(Idents(Human), Human$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))

library(ggplot2)
library(RColorBrewer)
fills <- colorRampPalette(brewer.pal(5, "Set3"))(15)

ggplot(Cellratio) + 
  geom_bar(aes(x = Freq, y= Var2, fill = Var1),stat = "identity",width = 0.5,size = 0.1,colour = '#222222')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Sample')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.1, linetype="solid"))+
  scale_fill_manual(values = rev(fills),
                    name = "celltype",
                    guide = guide_legend(ncol = 1,
                                         keywidth = 1,
                                         keyheight = 1))
#给meta.data添加一列，并暂时赋予空值
Human@meta.data$celltype="NA"
celltype=data.frame(ClusterID=c('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'),celltype=celltype, stringsAsFactors = F)
#赋值
for(i in 1:nrow(celltype)){Human@meta.data[which(Human@meta.data$seurat_clusters==celltype$ClusterID[i]),'celltype']<-celltype$celltype[i]}

####提取巨噬细胞进一步分群
Epithelial<- Human[,Human@meta.data$celltype %in% c('Epithelial cell')]
#计算 分群
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)##求表达量
Epithelial<- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)###寻找高变基因的数量

#(3)数据缩放
all.genes <- rownames(Epithelial)
Epithelial <- ScaleData(Epithelial, features = all.genes)

#(4)线性降维
#执行PCA线性降维
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))###PCA分析

DimPlot(Epithelial, reduction = "pca",group.by = 'orig.ident')

##高阶PCA分析
Epithelial <-RunHarmony(Epithelial,group.by.vars = "orig.ident")


#####使用harmony后的结果挑选合适的PC数聚类
ElbowPlot(Epithelial, ndims = 50, reduction = "harmony")####根据散点图和CD文献选择拐点PC数16
Epithelial <- FindNeighbors(Epithelial, reduction= "harmony",dim=1:12) %>% FindClusters(resolution = 0.4)##dim=1:16拐点数##0.5分成12簇  1.2分成16
####一定要指定harmony(reduction="harmony"),要不然harmony白做
table(Epithelial@meta.data$seurat_clusters)##查看每一类有多少个细胞，可以看到一共有17个簇


####UMAP的绘制
Epithelial <- RunUMAP(Epithelial,reduction= "harmony", dim=1:20)
plot1= DimPlot(Epithelial,reduction = "umap", label=T)
plot2= DimPlot(Epithelial,reduction = "umap", group.by = "orig.ident")

plot1 + plot2

# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
Epithelial.markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
write.csv(Epithelial.markers ,"Epithelial.markers.csv")

#(9)利用先验知识定义细胞类型
#加上注释
new.cluster.ids1 <- c("Alveolar epithelial cell Type 2", 'SLC16A7+ cell','Club cell','Alveolar pneumocyte Type I','Epithelial cell',
                      'Lung cell','Basal cell',"Macrophage","T cell",
                      "Fibroblast","Ciliated cell")

names(new.cluster.ids1) <- levels(Epithelial)
Epithelial<- RenameIdents(Epithelial, new.cluster.ids1)
P1 <- DimPlot(Epithelial, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()
P2 <- DimPlot(Epithelial,reduction = "umap", group.by = "orig.ident",pt.size = 0.5)

P1+P2
DimPlot(Macrophage,split.by = "orig.ident")

###组间细胞差异的展示
table(Epithelial$orig.ident)
prop.table(table(Idents(Epithelial)))
table(Idents(Epithelial), Epithelial$orig.ident)#各组不同细胞群细胞数

Cellratio <- prop.table(table(Idents(Epithelial), Epithelial$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))

library(ggplot2)
library(RColorBrewer)
fills <- colorRampPalette(brewer.pal(5, "Set3"))(15)

ggplot(Cellratio) + 
  geom_bar(aes(x = Freq, y= Var2, fill = Var1),stat = "identity",width = 0.5,size = 0.1,colour = '#222222')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Sample')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA, color="black", size=0.1, linetype="solid"))+
  scale_fill_manual(values = rev(fills),
                    name = "celltype",
                    guide = guide_legend(ncol = 1,
                                         keywidth = 1,
                                         keyheight = 1))
