#########################################################################
# Script for creating dataframe of p-values when 11 clinical phenotypes #
#               are regressed for all gene expression data              #
#########################################################################

# load clinical phenotype data
pheno_data <- read.csv("usable_clin_pheno_data.csv", as.is = TRUE)

# load library
library(dplyr)

# create matrix of the 11 clinical phenotypes
clin_pheno11 <- pheno_data[,c(1,2,11,12,31,54,32,126,151,152,55,56,57)]

# check if the phenotypes require normalization
quartz()
par(mfrow=c(3,4))
for(i in 1:11){
  hist(clin_pheno11[,i+2], breaks = 20, main = colnames(clin_pheno11)[i+2], xlab =  colnames(clin_pheno11)[i+2])
}

quartz()
par(mfrow=c(3,4))
for(i in 1:11){
  hist(log(clin_pheno11[,i+2]), breaks = 20, main = paste("log transformed", colnames(clin_pheno11)[i+2], sep = " "), xlab =  colnames(clin_pheno11)[i+2])
}

# log transform the 11 clinical phenotypes to normalize the data
for (i in 1:11){
  clin_pheno11[,i+2] <- log(clin_pheno11[,i+2])
}

# replace all infinite terms with NA
clin_pheno11$Gcg_secreted[which(is.infinite(clin_pheno11$Gcg_secreted))] <- NA

# function to compute adjusted pvalues adjusted for one covariate
pvalue_1cov <- function(x,f,c){
  fit1 <- lm(x ~ c+f)
  fit0 <- lm(x ~ c)
  anova(fit0, fit1)[2,6]
}

# create dataframe
df1 <- vector('numeric')
pvalue_1cov_allgene <- function(clin_pheno11_number){
  for(i in 1:21771){
    print(i)
    df1 <- c(df1, pvalue_1cov(clin_pheno11[,clin_pheno11_number], rankz.mrna[,i], pheno_data$sex))
  }
  return(p.adjust(df1, method="BH"))
}

df2 <- data.frame(matrix(ncol = 0, nrow = 21771), stringsAsFactors = FALSE)
for(i in 1:11){
  df2 <- cbind(df2, pvalue_1cov_allgene(i+2))
  colnames(df2)[i] <- colnames(clin_pheno11)[i+2]
}

df2 <- cbind(annot.mrna$id, annot.mrna$symbol, df2)

# exportation of newly created csv file
write.csv(x=df2, file="Pvaluetable_11clinpheno.csv", row.names = FALSE)

