install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)
library(ggplot2)
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# function to retrieve rankz.mrna data for a given gene
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

# create list of genes that are of interest and create normalized phenotype data
genelist <- c("Heatr6", "Ddx52", "Hnf1b", "Ccl6", "Synrg")

Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)

# find correlation of the genes of interest and Ins_tAUC/Ins_0min
for(i in 1:5){
  x <- cor(gene.expression.rankz(genelist[i]), Ins_tAUC_norm)
  y <- cor(gene.expression.rankz(genelist[i]), Ins_0min_norm)
  print("********************")
  print(genelist[i])
  print(paste("Ins_tAUC correlation:", x, sep = " "))
  print(paste("Ins_0min correlation:", y, sep = " "))
}

# None of the genes, when added into the mixed model as covariate, drops the peak on chr11 by a significant amount.

#####################################################
# Looking at the genes of interest in the BTBR data #
#####################################################
rm(list = ls())

# load BTBR data
load("/Users/s-ishimt/Desktop/Jax_SSP_17/BTBR/BTBR_Data/BTBR.clean.data.Rdata")

# create list of genes that are of interest and create normalized phenotype data
genelist <- c("Heatr6", "Ddx52", "Hnf1b", "Ccl6", "Synrg")

INS_4wk <- phenotypes.rz$INS.4wk
INS_6wk <- phenotypes.rz$INS.6wk
INS_8wk <- phenotypes.rz$INS.8wk
INS_10wk <- phenotypes.rz$INS.10wk


for(i in 1:5){
  islet_cor <- cor(islet.rz[,which(annot$gene_symbol == genelist[i])], INS_10wk, use = "complete.obs")
  adipose_cor <- cor(adipose.rz[,which(annot$gene_symbol == genelist[i])], INS_10wk, use = "complete.obs")
  gastroc_cor <- cor(gastroc.rz[,which(annot$gene_symbol == genelist[i])], INS_10wk, use = "complete.obs")
  hypo_cor <- cor(hypo.rz[,which(annot$gene_symbol == genelist[i])], INS_10wk, use = "complete.obs")
  kidney_cor <- cor(kidney.rz[,which(annot$gene_symbol == genelist[i])], INS_10wk, use = "complete.obs")
  liver_cor <- cor(liver.rz[,which(annot$gene_symbol == genelist[i])], INS_10wk, use = "complete.obs")
  print("************************")
  print(genelist[i])
  print(paste("islet tissue correlation:", islet_cor, sep = " "))
  print(paste("adipose tissue correlation:", adipose_cor, sep = " "))
  print(paste("gastroc tissue correlation:", gastroc_cor, sep = " "))
  print(paste("hypo tissue correlation:", hypo_cor, sep = " "))
  print(paste("kidney tissue correlation:", kidney_cor, sep = " "))
  print(paste("liver tissue correlation:", liver_cor, sep = " "))
}


for(i in 1:5){
  islet_cor <- cor(islet.rz[,which(annot$gene_symbol == genelist[i])], INS_4wk, use = "complete.obs")
  adipose_cor <- cor(adipose.rz[,which(annot$gene_symbol == genelist[i])], INS_4wk, use = "complete.obs")
  gastroc_cor <- cor(gastroc.rz[,which(annot$gene_symbol == genelist[i])], INS_4wk, use = "complete.obs")
  hypo_cor <- cor(hypo.rz[,which(annot$gene_symbol == genelist[i])], INS_4wk, use = "complete.obs")
  kidney_cor <- cor(kidney.rz[,which(annot$gene_symbol == genelist[i])], INS_4wk, use = "complete.obs")
  liver_cor <- cor(liver.rz[,which(annot$gene_symbol == genelist[i])], INS_4wk, use = "complete.obs")
  print("************************")
  print(genelist[i])
  print(paste("islet tissue correlation:", islet_cor, sep = " "))
  print(paste("adipose tissue correlation:", adipose_cor, sep = " "))
  print(paste("gastroc tissue correlation:", gastroc_cor, sep = " "))
  print(paste("hypo tissue correlation:", hypo_cor, sep = " "))
  print(paste("kidney tissue correlation:", kidney_cor, sep = " "))
  print(paste("liver tissue correlation:", liver_cor, sep = " "))
}