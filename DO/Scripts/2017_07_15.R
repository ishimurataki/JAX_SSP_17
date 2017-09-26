install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

G167_G33_lm <- lm(G167_ins_sec_norm ~ G33_ins_sec_norm, na.action = na.exclude)
G167_G33_res <- as.data.frame(G167_G33_lm$residuals)

which(c(1:378) %in%  rownames(G167_G33_res) == FALSE)

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

newrow <- NA

G167_G33_res <- insertRow(G167_G33_res, newrow, 66)
G167_G33_res <- insertRow(G167_G33_res, newrow, 67)
G167_G33_res <- insertRow(G167_G33_res, newrow, 68)

G167_G33_res <- G167_G33_res[,1]

for(i in 1:21771){
  if(abs(cor(rankz.mrna[,i], G167_G33_res, use = "complete.obs")) > 0.4){
    print("********************")
    print(cor(rankz.mrna[,i], G167_G33_res, use = "complete.obs"))
    print(annot.mrna$symbol[i])
  }
}






