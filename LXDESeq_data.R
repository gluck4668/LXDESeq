
setwd("D:/Desktop/Råå¼å/LXDESeq")
library(openxlsx)

gene_set_example <- read.xlsx("GSE205390_gene_expression_Set.xlsx")

usethis::use_data(gene_set_example,overwrite = T)

rm(list=ls())

data(gene_set_example)


