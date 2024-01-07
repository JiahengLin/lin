rm(list=ls())
pacman::p_load(tidyverse,microeco,magrittr)
#install.packages("readxl")
#install.packages("openxlsx")
rm(list = ls())
{library(readxl) 
  library(openxlsx)
  library(MicrobiotaProcess)
  library(phyloseq)
  library(ggtree)
  library(ggtreeExtra)
  library(tidyverse)}
library(microbiomeMarker)
###矩阵准备
data  <- read.csv("E:/metabolic/Smoker/LEfSeSmoker.csv",head=T,sep=",")
head(data)
colnames(data) <- data[1,]
data <- data[-1,]
otu <- data.frame(data[,-1], row.names=data[,1])
otu[c(1:6),c(1:6)]
otu <- as.matrix(otu)
class(otu) <- "numeric"
otu <- data.frame(otu)
names(otu)<-gsub("\\Healthy","Non-Smoker",names(otu))
##分组信息
group<- data.frame( row.names = colnames(otu),
                    group = c(rep("Non-Smoker",9),rep("Smoker",6)))  

##菌群名称
library(tidyr)
df<-data.frame(gsub("\\|",",",data[,1]))
names(df) <- "ID"
df<-separate(df, ID, into = c("d","p","c","o","f","g","s"),sep = ",")
rownames(df) <- data[,1]
###计算
head(otu)[1:6,1:6]; head(group)[1:2, ]; head(df)[,1:6]
dataset <- microtable$new(sample_table = group,
                          otu_table = otu, 
                          tax_table = df)
lefse <- trans_diff$new(dataset = dataset, 
                        method = "lefse", 
                        group = "group", 
                        alpha = 0.05, 
                        lefse_subgroup = NULL,
                        p_adjust_method = "none")
##画图
lefse$plot_diff_cladogram(use_taxa_num = 500, 
                          use_feature_num = 50, 
                          clade_label_level = 5, 
                          group_order = c("Non-Smoker", "Smoker"))


rm(list = ls())
library(vegan)
library(ggplot2)
data  <- read.table("E:/microbe/Smoker/bray_curtis_distance_matrix.csv",head=T,sep=",")
rownames(data) <- data[,1]
data <- data[,-1]
pcoa <- cmdscale (data,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p
###加组
group <- read.table("E:/microbe/Smoker/group.csv",head=T,sep=",")
df <- merge(pc12,group,by="samples")
color=c("#98F5FF","#EE3B3B")
p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=group))+
  theme_bw()+
  geom_point(size=4)+
  theme(panel.grid = element_blank())+
  ylim(-0.6,1.2)+
  #geom_vline(xintercept = 0,lty="dashed")+
  #geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=4)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PC1 ",pc[1],"%"),
       y=paste0("PC2 ",pc[2],"%"))+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#1597A5","#EE3B3B"))+
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20,angle=90),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20),
        panel.grid=element_blank())
p1

p1 + stat_ellipse(data=df,geom = "polygon",level=0.95,linetype = 1,size=1,aes(fill=group),alpha=0.1,show.legend = T)

library(ggplot2)
library(pheatmap)
library(scales)
library(plyr)
express<- read.table("E:/microbe/COPD/Group.COPD_vs_Healthy.diff.csv",head=T,sep=",")
express <- express[c(1:11),c(1, 6, 10:27)]
matrix <- express[,c(2:20)]
rownames(matrix) <- express[,1]
data<-matrix[,c(2:19)]
data <- data [,c(10:18,1:9)]
head(data)
####除去符号
names(data)<-gsub("\\_PH","",names(data))
head(data)
pheatmap(data)
annotation_row <- express[,2]
annotation_row <- data.frame(annotation_row)
rownames(annotation_row) <- express[,1]
names(annotation_row)[1] <- 'Class'
head(annotation_row)
group<- data.frame(samples = colnames(data),
                   group = c(rep("Healthy",9),rep("COPD",9)))
annotation_col<- data.frame(group[,2])
rownames(annotation_col) <- group[,1]
names(annotation_col)[1] <- 'Group'
head(annotation_col)
pheatmap(data, annotation_row = annotation_row)
ann_colors = list(Group = c(Healthy='#b0d45d', COPD = '#ffe788'), Class = c(down='#76afda',up='#f06152'))
pheatmap(data, annotation_row = annotation_row, annotation_col = annotation_col, 
         cluster_cols = F, cluster_rows = T, scale = "row", color = colorRampPalette(c("#ffc556", "white", "#fc9272" ))(40),
         annotation_colors = ann_colors,border_color = "white", cellwidth = 32, cellheight = 20)

##单个样本
genetop<- read.table("E:/lncRNAseq/Genus_abund_top30.csv",head=T,sep=",")
matrix<-genetop[,c(2:25)]
rownames(matrix)=genetop[,1]
matrix<-matrix[,c(1:15)]
group<- data.frame(samples = colnames(matrix),
                   group = c(rep("Healthy",9),rep("Smoker",6)),
                   group_number = c(rep("1", 9),rep("2",6)))
data3 <- matrix
library(reshape2)
library(ggplot2)
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))
data4 <- melt(data3, id='Taxonomy')
names(data4)[2] <- 'samples'
data5 <- merge(data4, group, by = 'samples')
p = ggplot( data5, aes( x = samples, weight = value, fill = Taxonomy))+
  geom_bar( position = "stack")
##分组样本
genetop<- read.table("E:/microbe/COPD/Phylum_abund_Group.csv",head=T,sep=",")
matrix<-genetop[,c(2:3)]
rownames(matrix)=genetop[,1]
matrix<-matrix[,c(2,1)]
data3 <- matrix
library(reshape2)
library(ggplot2)
data3$Taxonomy <- factor(rownames(data3), levels = rev(rownames(data3)))
data4 <- melt(data3, id='Taxonomy')
names(data4)[2] <- 'samples'
names(data4)[1] <- 'Phyla'
p = ggplot( data4, aes( x = samples, weight = value, fill = Phyla))+
  geom_bar( position = "stack")
color=c("#4a1486","#7CC767","#d73027","#7A142C","#167153","#eb998b",
        "#562e3c","#dcf2ff","#35212e","#562e3c","#a14462","#eb998b",
        "#fddbc8","#42465c","#356d67","#4c9568", "#7fb961","#b0d45d",
        "#ffe788","#b20000", "#f06152","#7d4444","#9e6c69","#cca69c",
        "#5066a1","#76afda","#abddff","#e8743c","#ffc556","#8264CC","#E0367A")
Smoker <- p +  scale_fill_manual( values = color)+
  labs(x=paste0(""),y=paste0("Relative Abundance(%)"))+#将x、y轴标题改为贡献度
  theme_bw()+theme(panel.border = element_blank(),panel.grid=element_blank(),
                   axis.line.y = element_line(colour = "black",size = 1))+#隐藏网格线
  scale_y_continuous(expand = c(0,0)) 
Smoker
ggsave(Smoker, filename = 'IPvsHealthy.pdf')

library(ggplot2)
alph<- read.table("E:/microbe/Smoker/alpha.csv",head=T,sep=",")
matrix <- alph[,c(2:7)]
rownames(matrix) <- alph[,1]
group <- factor(c(rep("Healthy",9),rep("Smoker",6)))
table(group)
alph$group <- group
p = ggplot(alph, aes(x=group, y=shannon)) + 
  geom_boxplot()

p+theme_classic() + scale             


