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
 