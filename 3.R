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
         
