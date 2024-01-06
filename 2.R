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

