install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)
library(ggplot2)
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# find all genes under the chr 11 peak from 80 - 88 mbp
chr11intervalgenes <- annot.mrna[which(annot.mrna$chr == "11"), ]
chr11intervalgenes <- chr11intervalgenes[which(chr11intervalgenes$middle_point > 80*10^6 &
                                                        chr11intervalgenes$middle_point < 88*10^6), ]

# creating dataframe with selected chr11 genes:
chr11intervalgenesdf <- rankz.mrna[,which(colnames(rankz.mrna) %in% chr11intervalgenes$id)]
colnames(chr11intervalgenesdf) <- chr11intervalgenes$symbol

# finding correlation between every selected chr11 gene and Ins_tAUC and Ins_0min
Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)

Ins_tAUCcor <- vector("numeric")
for(i in 1:133){
  cor <- cor(Ins_tAUC_norm, chr11intervalgenesdf[,i], use = "complete.obs")
  Ins_tAUCcor <- c(Ins_tAUCcor, cor)
}

Ins_0mincor <- vector("numeric")
for(i in 1:133){
  cor <- cor(Ins_0min_norm, chr11intervalgenesdf[,i], use = "complete.obs")
  Ins_0mincor <- c(Ins_0mincor, cor)
}

chr11expression_insulin_cordf <- cbind(Ins_tAUCcor, Ins_0mincor)
rownames(chr11expression_insulin_cordf) <- chr11intervalgenes$symbol

#####################################################
# Isolated islet secretion phenotypes p-value table #
#####################################################

Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)
G33_ins_sec_norm <- log(pheno_data$G33_ins_secrete)
G83_ins_sec_norm <- log(pheno_data$G83_ins_secrete)
G167_ins_sec_norm <- log(pheno_data$G167_ins_secrete)
G83_FC_ins_secrete <- log(pheno_data$G83_FC_ins_secrete)
G167_FC_ins_secrete <- log(pheno_data$G167_FC_ins_secrete)
G33_fract_ins_secrete <- log(pheno_data$G33_fract_ins_secrete)
G83_fract_ins_secrete <- log(pheno_data$G83_fract_ins_secrete)
G167_fract_ins_secrete <- log(pheno_data$G167_fract_ins_secrete)

# creating dataframe with the selected clinical phenotypes
ins_test_clinphenos <- cbind(pheno_data[,1:2], 
                             Ins_tAUC_norm, Ins_0min_norm, G33_ins_sec_norm, 
                             G83_ins_sec_norm, G167_ins_sec_norm, G83_FC_ins_secrete, 
                             G167_FC_ins_secrete, G33_fract_ins_secrete, 
                             G83_fract_ins_secrete, G167_fract_ins_secrete)

# function to compute adjusted pvalues adjusted for one covariate
pvalue_1cov <- function(x,f,c){
  fit1 <- lm(x ~ c+f)
  fit0 <- lm(x ~ c)
  anova(fit0, fit1)[2,6]
}

# create dataframe
df1 <- vector('numeric')
pvalue_1cov_allgene <- function(ins_test_clinphenos_number){
  for(i in 1:21771){
    print(i)
    df1 <- c(df1, pvalue_1cov(ins_test_clinphenos[,ins_test_clinphenos_number], rankz.mrna[,i], pheno_data$sex))
  }
  return(p.adjust(df1, method="BH"))
}

df2 <- data.frame(matrix(ncol = 0, nrow = 21771), stringsAsFactors = FALSE)
for(i in 3:12){
  df2 <- cbind(df2, pvalue_1cov_allgene(i))
  colnames(df2)[i-2] <- colnames(ins_test_clinphenos)[i]
}

dftest <- cbind(annot.mrna, df2)

# exportation of newly created csv file
write.csv(x=dftest, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_19/Pvaluetable_Ins_test.csv", row.names = FALSE)

#######################################################################################################
# Snps associated with both high LOD scores and changes in nucleic sequences in chr11 interval genes  #
#######################################################################################################
# read in SNP association datafile
Ins_tAUC_topSNPs <- read.csv(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_18/top_SNPs_Ins_tAUC_chr11.csv")
Ins_0min_topSNPs <- read.csv(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_18/top_SNPs_Ins_0min_chr11.csv")

# create matching SNP data for tAUC and 0min, concatenate into one dataframe
Ins_tAUC_topSNPs <- Ins_tAUC_topSNPs[which(Ins_tAUC_topSNPs$snp %in% Ins_0min_topSNPs$snp), -c(1,6,9,10,13), drop = TRUE]
colnames(Ins_tAUC_topSNPs)[9] <- "Ins_tAUC_lod"
Ins_0min_topSNPs<- Ins_0min_topSNPs[which(Ins_0min_topSNPs$snp %in% Ins_tAUC_topSNPs$snp), -c(1,6,9,10,13), drop = TRUE]
Ins_0min_lod <- Ins_0min_topSNPs[match(Ins_tAUC_topSNPs$snp_id, Ins_0min_topSNPs$snp_id),9,drop = TRUE]
tAUC_0min_SNPs_lod <- cbind(Ins_tAUC_topSNPs, Ins_0min_lod)

# importating sanger snp data for genes between 80-86 mbp
sangersnpdata <- read.csv(file = "https://www.sanger.ac.uk/sanger/Mouse_SnpViewer_Export/SNPs/rel-1303?sv=31;sn=512;st=156241930;loc=11:80000000-86000000;format=csv")

# find SNPS that are in both datasets
sangersnpdata <- sangersnpdata[which(sangersnpdata$dbSNP %in% tAUC_0min_SNPs_lod$snp_id),]
tAUC_0min_SNPs_lod <- tAUC_0min_SNPs_lod[which(tAUC_0min_SNPs_lod$snp_id %in% sangersnpdata$dbSNP),]
sangersnpdata <- sangersnpdata[match(tAUC_0min_SNPs_lod$snp_id, sangersnpdata$dbSNP), -c(1,2,4)]
associatedSNPS <- cbind(tAUC_0min_SNPs_lod, sangersnpdata)

# exportation of newly created csv file
write.csv(x=associatedSNPS, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_19/associatedSNPS.csv", row.names = FALSE)
