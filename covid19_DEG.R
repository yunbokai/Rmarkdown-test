# covid19_DEG

library(tidyverse)
library(ggpubr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(hdf5r)
library(tidyverse)
library(magrittr)
library(EnhancedVolcano) #火山图的包 不想要不安也行
library(dplyr)
library(DESeq2)

rm(list = ls())
Exp <- as.data.frame(fread("GSE147507_RawReadCounts_Human.tsv"))
dim(Exp)
# [1] 21797    79
rownames(Exp) <- Exp$V1
Exp <- Exp[,-1]

# 所有列名,看一下数据样本的组成情况
colnames(Exp)

dim(Exp)
# [1] 21797    78

###看一下Series的组成情况
series <- data.frame(Series=c(rep(NA,78)),Sample=c(rep(NA,78)))
for (i in (1:length(colnames(Exp)))){
  series[i,1] <- unlist(strsplit(colnames(Exp)[i], "_"))[1]
  series[i,2]<- colnames(Exp)[i]
}
table(series$Series)#大的分组应该就是这几个了，拿出来分别分析
# Series1 Series15 Series16  Series2  Series3  Series4  Series5  Series6  Series7  Series8  Series9 
# 6        4        9        6        4        4        6        6        6        9       18 

###差异分析一组一组看吧，想分析哪些数据就自己把他们选出来
###因为每个series的数据组成都还有点不一样，没法写成函数了

#---------------------------------------------------------------------------------

###举例1：比如要第1个series的两组比较
colnames(Exp)
##观察数据，其实就是Series1_NHBE_Mock & Series1_NHBE_SARS-CoV-2各有三个生物学重复
myseries <- series%>%filter(Series == paste0("Series",1)) #choose your group here, for example:1

#condition信息，其实就是分组信息
myseries$Condition <- c(rep(NA,length(myseries$Series)))
for (i in(1:length(myseries$Series))) {
  myseries$Condition[i] <- unlist(strsplit(myseries$Sample[i], "_"))[length(unlist(strsplit(myseries$Sample[i], "_")))-1]
}

##DESeq2的conditon因子构建，代表分组，比如我这里就就是group Mock和 group SARS-CoV-2
condition <- factor(myseries$Condition)

##构建表达矩阵数据data，按分析需要，提取我们需要的表达矩阵
data <- Exp[,myseries$Sample]

##开始差异分析 DESeq2
dds <- DESeqDataSetFromMatrix(countData =data, 
                              DataFrame(condition), design= ~ condition )#构建DES数据类型
dds <- DESeq(dds)
save(dds,file = "Series1_Mock_vs_SARS_dds.rda")#储存dds结果
resultsNames(dds)#[1] "Intercept"                    "condition_SARS.CoV.2_vs_Mock"即结果是SARS.CoV.2_vs_Mock的形式展示的
res <- results(dds)
summary(res)
table(res$padj<0.05) #可取P值小于0.05的结果
res[order(res$padj), ] %>% head   #查看排在前几位的差异情况 

##随意筛选了一下结果，可以自己确定结果筛选或者处理dds结果
res <- res[order(res$padj),]# 按照 padj 值进行排序
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
##得到csv格式的差异表达分析结果
write.csv(resdata,file= "Series1_Mock_vs_SARS.csv",row.names = F)

##火山图简单看一下结果
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
#---------------------------------------------------------------------------------


###举例2：比如要第16个series的三组比较
##观察数据，其实就是Series16_A549-ACE2_Mock & Series16_A549-ACE2_SARS-CoV-2 &Series16_A549-ACE2_SARS-CoV-2_Rux各有三个生物学重复
myseries <- series%>%filter(Series == paste0("Series",16))#choose your group here: choose 16

##确定分组信息，把分组的Mock SARS-CoV-2 Rux信息提取出来
myseries$Condition <- c(rep(NA,length(myseries$Series)))
for (i in(1:length(myseries$Series))) {
  myseries$Condition[i] <- unlist(strsplit(myseries$Sample[i], "_"))[length(unlist(strsplit(myseries$Sample[i], "_")))-1]
}

##DESeq2的conditon因子构建，代表分组，比如我这里就就是group Mock和 group SARS-CoV-2和 group Rux
condition <- factor(myseries$Condition)

##构建表达矩阵数据data，把需要差异分析的表达矩阵提取出来
data <- Exp[,myseries$Sample]

##开始差异分析
dds <- DESeqDataSetFromMatrix(countData =data, 
                              DataFrame(condition), design= ~ condition )#构建DES数据类型
dds <- DESeq(dds)
save(dds,file = "Series16_Mock_vs_SARS_vs_Rux_dds.rda")#储存dds结果
resultsNames(dds)
# [1] "Intercept"                    "condition_Rux_vs_Mock"        "condition_SARS.CoV.2_vs_Mock" 可以看到其实是Rux组和SARS-CoV-2组分别和Mock组进行了比较

##三个分组和两个分组在提取结果的时候略有不同,要分别提取结果
res1 <- results(dds, name="condition_Rux_vs_Mock",tidy=TRUE)   
res2 <- results(dds, name="condition_SARS.CoV.2_vs_Mock",tidy=TRUE)    
summary(res1)   ##结果的简单统计
summary(res2)   ##结果的简单统计

#对res1分析
res <- res1
table(res$padj<0.05) #可取P值小于0.05的结果
res[order(res$padj), ] %>% head   #查看排在前几位的差异情况 

#随意筛选了一下结果，可以自己确定结果筛选或者处理dds结果
res <- res[order(res$padj),]# 按照 padj 值进行排序
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
#得到csv格式的差异表达分析结果
write.csv(resdata,file= "Series16_Rux_vs_Mock.csv",row.names = F)

#火山图简单看一下结果

EnhancedVolcano(res,
                lab =res$row,
                x = 'log2FoldChange',
                y = 'pvalue')

#对res2分析
res <- res2
table(res$padj<0.05) #可取P值小于0.05的结果
res[order(res$padj), ] %>% head   #查看排在前几位的差异情况 

#随意筛选了一下结果，可以自己确定结果筛选或者处理dds结果
res <- res[order(res$padj),]# 按照 padj 值进行排序
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
#得到csv格式的差异表达分析结果
write.csv(resdata,file= "Series16_SARS.CoV.2_vs_Mock.csv",row.names = F)

#火山图简单看一下结果
EnhancedVolcano(res,
                lab = res$row,
                x = 'log2FoldChange',
                y = 'pvalue')
