#########################################
# Batch Mediation Analysis and Plotting #
#########################################

install.packages("devtools")
library(devtools)
library(qtl2)
library(qtl2convert)
library(ggplot2)
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/mediation.scan.R")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/plot.mediation.R")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/gmb.coordinates.R")

load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/mouse.chrlen.rda")
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)
rownames(pheno_data) <- pheno_data[,1]

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

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)
addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# creating elements for mediation analysis of chr 11 peak on Ins_tAUC
Ins_tAUC_norm <- log(pheno_data[,18])
qtl11 <- genoprobs[,,dimnames(genoprobs)[[3]] == "11_83591414"]

# creating list of all the elements
Ins_tAUC_mediationDF <- list(target = Ins_tAUC_norm, mediator = rankz.mrna, 
                             annotation = annot.mrna, covar = addcovar, qtl.geno = qtl11)

# run mediation analysis
med <- mediation.scan(target = Ins_tAUC_mediationDF$target,
                      mediator = Ins_tAUC_mediationDF$mediator,
                      annotation = Ins_tAUC_mediationDF$annotation,
                      covar = Ins_tAUC_mediationDF$covar,
                      qtl.geno = Ins_tAUC_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med)                       

#####################################################
# Mediation Analysis for chr 14 peak on G33_ins_sec #
#####################################################
G33_ins_sec_norm <- log(pheno_data[,55])
qtl14 <- genoprobs[,,dimnames(genoprobs)[[3]] == "14_91585596"]

# creating list of all the elements
G33_ins_sec_mediationDF <- list(target = G33_ins_sec_norm, mediator = rankz.mrna, 
                             annotation = annot.mrna, covar = addcovar, qtl.geno = qtl14)

med <- mediation.scan(target = G33_ins_sec_mediationDF$target,
                      mediator = G33_ins_sec_mediationDF$mediator,
                      annotation = G33_ins_sec_mediationDF$annotation,
                      covar = G33_ins_sec_mediationDF$covar,
                      qtl.geno = G33_ins_sec_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med)   

#####################################################
# Mediation Analysis for chr 14 peak on G33_ins_sec #
#####################################################
G167_ins_sec_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,7,drop = FALSE],
                          kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
rownames(G167_ins_sec_qtl)[which.max(G167_ins_sec_qtl[,1])]
  
G167_ins_sec_norm <- log(pheno_data[,57])
qtl3 <- genoprobs[,,dimnames(genoprobs)[[3]] == "3_133549858"]

# creating list of all the elements
G167_ins_sec_mediationDF <- list(target = G167_ins_sec_norm, mediator = rankz.mrna, 
                                annotation = annot.mrna, covar = addcovar, qtl.geno = qtl3)

med <- mediation.scan(target = G167_ins_sec_mediationDF$target,
                      mediator = G167_ins_sec_mediationDF$mediator,
                      annotation = G167_ins_sec_mediationDF$annotation,
                      covar = G167_ins_sec_mediationDF$covar,
                      qtl.geno = G167_ins_sec_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med) 

# Cyp2u1 may be a gene of interest
Cyp2u1 <- gene.expression.rankz("Cyp2u1")
for(i in 3:12){
  print("*****************")
  print(cor(Cyp2u1, ins_test_clinphenos[,i], use = "complete.obs"))
  print(colnames(ins_test_clinphenos)[i])
}

# Trmt10a may be a gene of interest
Trmt10a <- gene.expression.rankz("Trmt10a")
for(i in 3:12){
  print("*****************")
  print(cor(Trmt10a, ins_test_clinphenos[,i], use = "complete.obs"))
  print(colnames(ins_test_clinphenos)[i])
}

G167_fract_ins_sec_norm <- log(pheno_data[,70])
qtl3 <- genoprobs[,,dimnames(genoprobs)[[3]] == "3_135570855"]

# creating list of all the elements
G167_fract_ins_sec_norm <- list(target = G167_fract_ins_sec_norm, mediator = rankz.mrna, 
                                 annotation = annot.mrna, covar = addcovar, qtl.geno = qtl3)

med <- mediation.scan(target = G167_fract_ins_sec_norm$target,
                      mediator = G167_fract_ins_sec_norm$mediator,
                      annotation = G167_fract_ins_sec_norm$annotation,
                      covar = G167_fract_ins_sec_norm$covar,
                      qtl.geno = G167_fract_ins_sec_norm$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med)  
