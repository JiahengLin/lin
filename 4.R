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


