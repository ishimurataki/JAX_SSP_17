install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)
library(ggplot2)
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# Establishing wanted clinical phenotypes:
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

# creating qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

#################################################################################################
# qtl scan of Ins_tAUC_norm and Ins_tAUC_norm regressed with different insulin secretion traits #
#################################################################################################
addcovar_G33_ins_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,5], data=annot.samples)[,-1]
addcovar_G83_ins_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,6], data=annot.samples)[,-1]
addcovar_G167_ins_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,7], data=annot.samples)[,-1]
addcovar_G83_FC_ins_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,8], data=annot.samples)[,-1]
addcovar_G167_FC_ins_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,9], data=annot.samples)[,-1]
addcovar_G33_fract_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,10], data=annot.samples)[,-1]
addcovar_G83_fract_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,11], data=annot.samples)[,-1]
addcovar_G167_fract_sec <- model.matrix(~Sex + Generation + diet_days + ins_test_clinphenos[,12], data=annot.samples)[,-1]

Ins_tAUC_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
Ins_tAUC_reg_G33_ins_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G33_ins_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G83_ins_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G83_ins_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G167_ins_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G167_ins_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G83_FC_ins_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G83_FC_ins_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G167_FC_ins_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G167_FC_ins_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G33_fract_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G33_fract_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G83_fract_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G83_fract_sec, cores = 4, reml=TRUE)
Ins_tAUC_reg_G167_fract_sec <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar_G167_fract_sec, cores = 4, reml=TRUE)


pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/tAUC_regressed.pdf", width = 20, height = 24)
par(mfrow = c(5,2))
plot_scan1(x = Ins_tAUC_qtl, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_qtl)[1])
plot_scan1(x = Ins_tAUC_reg_G33_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G33_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G83_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G83_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G167_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G167_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G83_FC_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G83_FC_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G167_FC_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G167_FC_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G33_fract_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G33_fract_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G83_fract_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G83_fract_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G167_fract_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G167_fract_ins_sec")
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/tAUC_regressedV2.pdf", width = 20, height = 12)
plot_scan1(x = Ins_tAUC_qtl, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_qtl)[1])
plot_scan1(x = Ins_tAUC_reg_G33_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G33_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G83_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G83_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G167_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G167_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G83_FC_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G83_FC_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G167_FC_ins_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G167_FC_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G33_fract_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G33_fract_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G83_fract_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G83_fract_ins_sec")
plot_scan1(x = Ins_tAUC_reg_G167_fract_sec, map = map, lodcolumn = 1, main = "Ins_tAUC|G167_fract_ins_sec")
dev.off()

quartz()
par(mfrow = c(4,2))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G33_ins_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G83_ins_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G167_ins_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G83_FC_ins_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G167_FC_ins_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G33_fract_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G83_fract_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))
plot_scan1(x = Ins_tAUC_qtl - Ins_tAUC_reg_G167_fract_sec, map = map, lodcolumn = 1, main = colnames(Ins_tAUC_reg_G33_ins_sec)[1], ylim = c(-3,3))

##########################################################
# qtl scan of Glu_tAUC_norm regressed with Ins_tAUC_norm #
##########################################################
Glu_tAUC_norm <- log(pheno_data[,16,drop = FALSE])
hist(Glu_tAUC_norm[,1])

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]
addcovar_Ins_tAUC <- model.matrix(~ Sex + Generation + diet_days + Ins_tAUC_norm, data = annot.samples)[,-1]

Glu_tAUC_qtl <- scan1(genoprobs = probs, pheno = Glu_tAUC_norm[,1,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
Glu_tAUC_reg_Ins_tAUC_qtl <- scan1(genoprobs = probs, pheno = Glu_tAUC_norm[,1,drop = FALSE],
                      kinship = K, addcovar = addcovar_Ins_tAUC, cores = 4, reml=TRUE)
quartz()
par(mfrow = c(2,1))
plot_scan1(x = Glu_tAUC_qtl, map = map, lodcolumn = 1, main = "Glu_tAUC")
plot_scan1(x = Glu_tAUC_reg_Ins_tAUC_qtl, map = map, lodcolumn = 1, main = "Glu_tAUC|ns_tAUC")

##########################################################
# more exploration of chr11 on qtl map of Ins_tAUC_norm #
##########################################################
Ins_tAUC_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
Ins_0min_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,4, drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

rownames(Ins_tAUC_qtl)[which(Ins_tAUC_qtl[,1] == max(Ins_tAUC_qtl[,1], na.rm = TRUE))]
rownames(Ins_0min_qtl)[which(Ins_0min_qtl[,1] == max(Ins_0min_qtl[,1], na.rm = TRUE))]

Ins_0min_qtl[which(rownames(Ins_0min_qtl) == "11_83591414"),1]

qtl11 <- genoprobs[,,dimnames(genoprobs)[[3]] == "11_83591414"]

anova(lm(Ins_tAUC_norm ~ pheno_data$sex + qtl11))
anova(lm(Ins_0min_norm ~ pheno_data$sex + qtl11))
anova(lm(Ins_tAUC_norm ~ pheno_data$sex + Ins_0min_norm + qtl11))
anova(lm(Ins_0min_norm ~ pheno_data$sex + Ins_tAUC_norm + qtl11))

for(i in 1:21771){
  if(abs(cor(rankz.mrna[,i], Ins_tAUC_norm, use = "complete.obs")) > 0.4){
    print("***************")
    print(annot.mrna$symbol[i])
    print(cor(rankz.mrna[,i], Ins_tAUC_norm, use = "complete.obs"))
  }
}

for(i in 1:21771){
  x <- anova(lm(rankz.mrna[,i] ~ pheno_data$sex + qtl11))
  if(x$`Pr(>F)`[2] < 0.000005 & 
     abs(cor(rankz.mrna[,i], Ins_tAUC_norm, use = "complete.obs")) > 0.3){
    print("********************")
    print(annot.mrna$symbol[i])
    print(x$`Pr(>F)`[2])
    print(cor(rankz.mrna[,i], Ins_tAUC_norm, use = "complete.obs"))
  }
}

##################################################
# Glucose stimulated insulin secretion slope qtl #
##################################################
G33_ins_sec_norm <- log(pheno_data$G33_ins_secrete)
G83_ins_sec_norm <- log(pheno_data$G83_ins_secrete)
G167_ins_sec_norm <- log(pheno_data$G167_ins_secrete)

Insulin_sec_df <- cbind(pheno_data[,1:2, drop = FALSE], G33_ins_sec_norm, G83_ins_sec_norm, G167_ins_sec_norm)
Insulin_sec_df <- rbind(c(NA,NA,3.3,8.3,16.7), Insulin_sec_df)

slopelist <- vector("numeric")
for(i in 2:379){
  lm <- lm(as.numeric(Insulin_sec_df[i,3:5]) ~ as.numeric(Insulin_sec_df[1,3:5]))
  slope <- as.numeric(lm$coefficients[2])
  slopelist <- c(slopelist, slope)
}
slopelist <- cbind(pheno_data[,1:2, drop = FALSE], slopelist)

slopes_qtl <- scan1(genoprobs = probs, pheno = slopelist[,3,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/Ins_secretion_slope_qtl.pdf", width = 15, height = 6)
plot_scan1(x = slopes_qtl, map = map, lodcolumn = 1, main = "Ins_secretion_slope")
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/Ins_secretion_slope_chr3_effectplot.pdf", width = 20, height = 12)
chr <- 3
qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = slopelist[,3,drop = FALSE],
                      kinship = K[[chr]], addcovar = addcovar, cores = 4)
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = slopes_qtl[,1, drop = FALSE], main = "Ins_secretion_slope")
dev.off()

#####################################
# More exploration of the chr3 peak #
#####################################
G167_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,7,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G83_FC_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,8,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G167_FC_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,9,drop = FALSE],
                         kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G167_frac_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,12,drop = FALSE],
                           kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

# Find confidence interval for chr 3 peak
G167_ins_qtlchr3 <- G167_ins_qtl[grep("3_", rownames(G167_ins_qtl)), , drop = FALSE]
G167_ins_qtlchr3 <- G167_ins_qtlchr3[-grep("13_", rownames(G167_ins_qtlchr3)), , drop = FALSE]
G83_FC_ins_qtlchr3 <- G83_FC_ins_qtl[grep("3_", rownames(G83_FC_ins_qtl)), , drop = FALSE]
G83_FC_ins_qtlchr3 <- G83_FC_ins_qtlchr3[-grep("13_", rownames(G83_FC_ins_qtlchr3)), , drop = FALSE]
G167_FC_ins_qtlchr3 <- G167_FC_ins_qtl[grep("3_", rownames(G167_FC_ins_qtl)), , drop = FALSE]
G167_FC_ins_qtlchr3 <- G167_FC_ins_qtlchr3[-grep("13_", rownames(G167_FC_ins_qtlchr3)), , drop = FALSE]
G167_frac_ins_qtlchr3 <- G167_frac_ins_qtl[grep("3_", rownames(G167_frac_ins_qtl)), , drop = FALSE]
G167_frac_ins_qtlchr3 <- G167_frac_ins_qtlchr3[-grep("13_", rownames(G167_frac_ins_qtlchr3)), , drop = FALSE]

rownames(G167_ins_qtlchr3)[which(G167_ins_qtlchr3[,1] == max(G167_ins_qtlchr3[,1]))]    #"3_133549858"
rownames(G83_FC_ins_qtlchr3)[which(G83_FC_ins_qtlchr3[,1] == max(G83_FC_ins_qtlchr3[,1]))]    #"3_135647724"
rownames(G167_FC_ins_qtlchr3)[which(G167_FC_ins_qtlchr3[,1] == max(G167_FC_ins_qtlchr3[,1]))]    #"3_137348141"
rownames(G167_frac_ins_qtlchr3)[which(G167_frac_ins_qtlchr3[,1] == max(G167_frac_ins_qtlchr3[,1]))]    #"3_140779374

G167_ins_qtlchr3[which(rownames(G167_ins_qtlchr3)== "3_135570855"),1]
G83_FC_ins_qtlchr3[which(rownames(G83_FC_ins_qtlchr3)== "3_135570855"),1]
G167_FC_ins_qtlchr3[which(rownames(G167_FC_ins_qtlchr3)== "3_135570855"),1]
G167_frac_ins_qtlchr3[which(rownames(G167_frac_ins_qtlchr3)== "3_135570855"),1] # best peak to use is 3_135570855

# get genoprobs for the locus on chr3
qtl3 <- genoprobs[,,dimnames(genoprobs)[[3]] == "3_135570855"]

# find genes that are highly correlated with both the genotype of the locus at chr 3 and the clinical phenotypes of interest
for(i in 1:21771){
  x <- anova(lm(rankz.mrna[,i] ~ pheno_data$sex + qtl3))
  if(x$`Pr(>F)`[2] < 0.05 & 
     abs(cor(rankz.mrna[,i], G167_ins_sec_norm, use = "complete.obs")) > 0.3){
    print("********************")
    print(annot.mrna$symbol[i])
    print(x$`Pr(>F)`[2])
    print(cor(rankz.mrna[,i], G167_ins_sec_norm, use = "complete.obs"))
  }
}

#####################################
# More exploration of the chr9 peak #
#####################################
G167_ins_qtlchr9 <- G167_ins_qtl[grep("9_", rownames(G167_ins_qtl)), , drop = FALSE]
G167_ins_qtlchr9 <- G167_ins_qtlchr9[-grep("19_", rownames(G167_ins_qtlchr9)), , drop = FALSE]
G167_FC_ins_qtlchr9 <- G167_FC_ins_qtl[grep("9_", rownames(G167_FC_ins_qtl)), , drop = FALSE]
G167_FC_ins_qtlchr9 <- G167_FC_ins_qtlchr9[-grep("19_", rownames(G167_FC_ins_qtlchr9)), , drop = FALSE]
G167_frac_ins_qtlchr9 <- G167_frac_ins_qtl[grep("9_", rownames(G167_frac_ins_qtl)), , drop = FALSE]
G167_frac_ins_qtlchr9 <- G167_frac_ins_qtlchr9[-grep("19_", rownames(G167_frac_ins_qtlchr9)), , drop = FALSE]

rownames(G167_ins_qtlchr9)[which(G167_ins_qtlchr9[,1] == max(G167_ins_qtlchr9[,1]))]    # "9_66804359"
rownames(G167_FC_ins_qtlchr9)[which(G167_FC_ins_qtlchr9[,1] == max(G167_FC_ins_qtlchr9[,1]))]    # "9_60117197"
rownames(G167_frac_ins_qtlchr9)[which(G167_frac_ins_qtlchr9[,1] == max(G167_frac_ins_qtlchr9[,1]))] # "9_67041196"

G167_ins_qtlchr9[which(rownames(G167_ins_qtlchr9)== "9_66804359"),1]
G167_FC_ins_qtlchr9[which(rownames(G167_FC_ins_qtlchr9)== "9_66804359"),1]
G167_frac_ins_qtlchr9[which(rownames(G167_frac_ins_qtlchr9)== "9_66804359"),1]

View(cbind(G167_ins_qtlchr9, G167_FC_ins_qtlchr9, G167_frac_ins_qtlchr9))

9_66804359
qtl9 <- genoprobs[,,dimnames(genoprobs)[[3]] == "9_66804359"]

for(i in 1:21771){
  x <- anova(lm(rankz.mrna[,i] ~ pheno_data$sex + qtl9))
  if(x$`Pr(>F)`[2] < 0.05 & 
     abs(cor(rankz.mrna[,i], G167_ins_sec_norm, use = "complete.obs")) > 0.3){
    print("********************")
    print(annot.mrna$symbol[i])
    print(x$`Pr(>F)`[2])
    print(cor(rankz.mrna[,i], G167_ins_sec_norm, use = "complete.obs"))
  }
}

Aph1b <- gene.expression.rankz("Aph1b")
cor(Aph1b, G167_ins_sec_norm, use = "complete.obs")
anova(lm(Aph1b ~ pheno_data$sex + qtl9))
anova(lm(G167_ins_sec_norm ~ pheno_data$sex + Aph1b))
anova(lm(G167_ins_sec_norm ~ pheno_data$sex + Aph1b + qtl9))

addcovar_Aph1b <- model.matrix(~Sex + Generation + diet_days + Aph1b, data=annot.samples)[,-1]
G167_ins_reg_Aph1b_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,7,drop = FALSE],
                      kinship = K, addcovar = addcovar_Aph1b, cores = 4, reml = TRUE)
plot_scan1(x = G167_ins_reg_Aph1b_qtl, map = map, lodcolumn = 1, main = "G167_ins_sec|Aph1b")


######################################
# More exploration of the chr14 peak #
######################################
Klf5 <- gene.expression.rankz("Klf5")

Ins_sec_Klf5_df <- cbind(pheno_data[,1:2, drop = FALSE], G33_ins_sec_norm, G33_fract_ins_secrete, G83_fract_ins_secrete, Klf5)

Klf5_qtl <- scan1(genoprobs = probs, pheno = Ins_sec_Klf5_df[,6,drop = FALSE],
                   kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/Klf5_qtl.pdf", width = 15, height = 6)
plot_scan1(x = Klf5_qtl, map = map, lodcolumn = 1, main = "Klf5")
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/Klf5_chr14_effectplot.pdf", width = 20, height = 12)
chr <- 14
qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = Ins_sec_Klf5_df[,6,drop = FALSE],
                      kinship = K[[chr]], addcovar = addcovar, cores = 4)
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = Klf5_qtl[,1, drop = FALSE], main = "Klf5")
dev.off()

addcovar_Klf5 <- model.matrix(~Sex + Generation + diet_days + Klf5, data=annot.samples)[,-1]
G33_ins_sec_qtl <- scan1(genoprobs = probs, pheno = Ins_sec_Klf5_df[,3,drop = FALSE],
                                  kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G33_ins_sec_reg_Klf5_qtl <- scan1(genoprobs = probs, pheno = Ins_sec_Klf5_df[,3,drop = FALSE],
                         kinship = K, addcovar = addcovar_Klf5, cores = 4, reml = TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_25/G33_ins_secregressed.pdf", width = 15, height = 4)
plot_scan1(x = G33_ins_sec_qtl, map = map, lodcolumn = 1, main = "G33_ins_sec")
plot_scan1(x = G33_ins_sec_reg_Klf5_qtl, map = map, lodcolumn = 1, main = "G33_ins_sec|Klf5")
dev.off()

rownames(G33_ins_sec_qtl)[G33_ins_sec_qtl[,1] == max(G33_ins_sec_qtl[,1], na.rm = TRUE)]

qtl14 <- genoprobs[,,dimnames(genoprobs)[[3]] == "14_91585596"]

anova(lm(G33_ins_sec_norm ~ pheno_data$sex + qtl14))
anova(lm(G33_ins_sec_norm ~ pheno_data$sex + Klf5 + qtl14))

anova(lm(Klf5 ~ pheno_data$sex + G33_ins_sec_norm + qtl14))

quartz()
ggplot(mapping = aes(x = Ins_sec_Klf5_df$Klf5, y = Ins_sec_Klf5_df$G33_ins_sec_norm, color = pheno_data$sex)) + geom_point() + 
  xlab("Klf5") + ylab("G33_ins_sec_norm") + 
  ggtitle(paste("G33_ins_sec_norm", "Klf5", sep = " vs "), 
          subtitle = cor( Ins_sec_Klf5_df$Klf5, Ins_sec_Klf5_df$G33_ins_sec_norm, use = "complete.obs"))

quartz()
ggplot(mapping = aes(x = Ins_sec_Klf5_df$Klf5, y = G83_ins_sec_norm, color = pheno_data$sex)) + geom_point() + 
  xlab("Klf5") + ylab("G83_ins_sec_norm") + 
  ggtitle(paste("G83_ins_sec_norm", "Klf5", sep = " vs "), 
          subtitle = cor(Ins_sec_Klf5_df$Klf5, G83_ins_sec_norm, use = "complete.obs"))
