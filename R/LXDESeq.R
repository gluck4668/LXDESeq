

LXDESeq <- function(gene_set){

all_packages <- data.frame(installed.packages())

pack <- data.frame(c("BiocManager","openxlsx","dplyr"))

bioc_pack <- data.frame(c("DESeq2",'factoextra'))

pack$type <- pack[,1] %in% all_packages$Package

for (i in 1:nrow(pack)){if (!requireNamespace(pack[i,1], quietly=TRUE))
  install.packages(pack[i,1],update = F,ask = F)}
rm(i)

for (i in 1:nrow(bioc_pack)){if (!requireNamespace(bioc_pack[i,1], quietly=TRUE))
  BiocManager::install (bioc_pack[i,1],update = F,ask = F) }

rm(i)

packages <- c(pack[,1],bioc_pack[,1])

for(i in packages){
  library(i, character.only = T)}

rm(i)

#--------------------------------------

if(dir.exists("analysis results")==F)
  dir.create("analysis results")

gene_set_df <- read.xlsx(gene_set)

gene_set_df <- na.omit(gene_set_df)

gene_df <- gene_set_df[,-1]

gene_df <- round(gene_df,0) # DESeq2 要求数据全部为整数，不能有小数

rownames(gene_df) <- gene_set_df[,1]

col_names <- colnames(gene_df)

condition <- factor(gsub("\\d+$", "", col_names))

group_df <- data.frame(row.names=col_names,condition)

dds <- DESeqDataSetFromMatrix(countData = gene_df, colData = group_df, design = ~ condition)


keep <- rowSums(counts(dds) >= 10) >= 3  #过滤低表达基因，至少有3个样品都满足10个以上的reads数  
dds <- dds[keep, ]


# 利用层次聚类的方法评估一下数据质量（vst标准化---hclust)

vsd <- vst(dds,blind = TRUE)    #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
sampleDists <- dist(t(assay(vsd))) 
res1 <- hcut(sampleDists, k = 2, stand = FALSE,hc_method ="average" ) 
fviz_dend(res1, 
          rect_fill = T,
          cex = 1,  # 字体大小
          color_labels_by_k=T, # 字体颜色
          horiz=T # 平行放置
          )


# 利用PCA主成分分析的方法评估一下数据质量。

rld <- vst(dds, blind=FALSE)     #vst()函数效果和rlog（）一样，且速度更快。
plotPCA(rld, intgroup="condition",ntop=500)


# 利用DESeq（）函数标准化dds矩阵；
dds1 <- DESeq(dds)    # 将数据标准化，必要步骤！！！
resultsNames(dds1)    # 查看结果的名称。
dds1$condition        #默认后者的处理组比前面的对照组。
res <- results(dds1)  # 必要步骤！！！
summary(res)          #看一下结果的概要信息，p值默认小于0.1。

# 提取差异分析结果。
table(res$padj < 0.05)        #padj 即矫正后的P值。看看有多少差异基因满足所设的P值要求。TRUE的数值为满足要求的基因个数。
res <- res[order(res$padj),]  #按照padj 进行升序排列
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds1, normalized=TRUE)),by="row.names",sort=FALSE)

resdata <- resdata[,c(-2,-4,-5)]

colnames(resdata)[1] <- "gene_symbol"

write.xlsx (resdata, "analysis results/DESeq2 analysis result.xlsx")   

#提取显著性差异的基因，并保存结果。
#大家通常是提取padj<0.05，Log2 (fold change)>1或者<-1的差异表达基因
#diff_gene_deseq2 <-subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))   #或者直接写为diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

#DEGs <- dplyr::filter(resdata, pvalue < 0.05 & abs(log2FoldChange) > 1)

DEGs <- dplyr::filter(resdata, pvalue < 0.05)

DEGs <- na.omit(DEGs)

up_DEGs <- dplyr::filter(DEGs, log2FoldChange >0)
down_DEGs <- dplyr::filter(DEGs, log2FoldChange <0)  
  
write.xlsx(DEGs,"analysis results/Differently expressed genes (DEGs).xlsx")  #生成的差异分析结果文件名为diff_gene_X666_CK.csv，可改动。

write.xlsx(up_DEGs,"analysis results/Differently expressed genes (DEGs)_up.xlsx")

write.xlsx(down_DEGs,"analysis results/Differently expressed genes (DEGs)_down.xlsx")

#print("---------------------------------------------------------------------")

print("Differently expressed genes were successfully analyzed, and the resuts can be found in the folder of <analysis results>")

}
