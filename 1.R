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
