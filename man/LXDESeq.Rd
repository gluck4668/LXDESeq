\name{LXDESeq::GEO gene set analysis}
\alias{LXDESeq}
\title{LXDESeq::GEO gene set analysis}
\usage{
LXDESeq(gene_set)
}
\description{
LXDESeq::GEO gene set analysis.
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXDESeq")

library(LXDESeq)

??LXDESeq

#-------------------
data(gene_set_example)

#------------------

rm(list=ls())

setwd("D:/Desktop/R_example/LXDESeq_example")

gene_set = "GSE205390_gene_expression_Set.xlsx"


LXDESeq(gene_set)

}
