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
Glu_tAUC_norm <- log(pheno_data$Glu_tAUC)
Ins_per_islet <- log(pheno_data$Ins_per_islet)

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

# function to get gene expression given gene name from rankz.mrna data
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

# Qtl analysis of G33_fract_ins_secrete and G83_fract_ins_secrete
G33_frac_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,10,drop = FALSE],
                          kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G83_frac_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,11,drop = FALSE],
                          kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G33_frac_ins_qtlchr14 <- G33_frac_ins_qtl[grep("14_", rownames(G33_frac_ins_qtl)),]
names(G33_frac_ins_qtlchr14)[which.max(G33_frac_ins_qtlchr14)]
# set maximum LOD peak position to 14_91792449
qtl14 <- genoprobs[,,dimnames(genoprobs)[[3]] == "14_91792449"]

# Batch Mediation of qtl14 on G33_fract_ins_secrete and G83_fract_ins_secrete
library(devtools)
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/mediation.scan.R")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/plot.mediation.R")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/gmb.coordinates.R")
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/mouse.chrlen.rda")

G33_ins_frac_mediationDF <- list(target = G33_fract_ins_secrete, mediator = rankz.mrna, 
                                annotation = annot.mrna, covar = addcovar, qtl.geno = qtl14)

med <- mediation.scan(target = G33_ins_frac_mediationDF$target,
                      mediator = G33_ins_frac_mediationDF$mediator,
                      annotation = G33_ins_frac_mediationDF$annotation,
                      covar = G33_ins_frac_mediationDF$covar,
                      qtl.geno = G33_ins_frac_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med)   

G83_ins_frac_mediationDF <- list(target = G83_fract_ins_secrete, mediator = rankz.mrna, 
                                 annotation = annot.mrna, covar = addcovar, qtl.geno = qtl14)

med <- mediation.scan(target = G83_ins_frac_mediationDF$target,
                      mediator = G83_ins_frac_mediationDF$mediator,
                      annotation = G83_ins_frac_mediationDF$annotation,
                      covar = G83_ins_frac_mediationDF$covar,
                      qtl.geno = G83_ins_frac_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med)  

# analysis of Stc1 gene expression
Stc1 <- gene.expression.rankz("Stc1")
Stc1 <- cbind(pheno_data[,1:2], Stc1)
Stc1_qtl <- scan1(genoprobs = probs, pheno = Stc1[,3,drop = FALSE],
                          kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
plot_scan1(x = Stc1_qtl, map = map, lodcolumn = 1, main = colnames(Stc1_qtl)[1])

chr <- 14
qtl.coef1 <- scan1coef(genoprobs = probs[,chr], pheno =  Stc1[,3,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
quartz()
plot(x = qtl.coef1, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = Stc1_qtl, main = "Stc1 Expression")

qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = Stc1[,3,drop = FALSE],
                      kinship = K[[chr]], addcovar = addcovar, cores = 4)
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = Stc1_qtl, main = "Stc1 Expression")

peak14_Stc1 <- genoprobs[,,dimnames(genoprobs)[[3]] == "14_68832847"] 
addcovar_Stc1_chr14peak <- model.matrix(~Sex + peak14_Stc1 + Generation + diet_days, data=annot.samples)[,-1]

Stc1_peak14cov_qtl <- scan1(genoprobs = probs, pheno = Stc1[,3,drop = FALSE],
                  kinship = K, addcovar = addcovar_Stc1_chr14peak, cores = 4, reml = TRUE)
quartz()
plot_scan1(x = Stc1_qtl, map = map, lodcolumn = 1, main = )

# Stc1 expression as an added covariate for Insulin fractional excretion clinical phenotypes
Stc1 <- gene.expression.rankz("Stc1")
addcovar_Stc1 <- model.matrix(~Sex + Stc1 + Generation + diet_days, data=annot.samples)[,-1]

G33_frac_ins_Stc1cov_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,10,drop = FALSE],
                          kinship = K, addcovar = addcovar_Stc1, cores = 4, reml = TRUE)
G83_frac_ins_Stc1cov_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,11,drop = FALSE],
                          kinship = K, addcovar = addcovar_Stc1, cores = 4, reml = TRUE)

quartz()
plot_scan1(x = G33_frac_ins_qtl, map = map, lodcolumn = 1, main = "G33 Fractional Ins Secretion | Stc1 Expression", col = "#FC81D0")
plot_scan1(x = G33_frac_ins_Stc1cov_qtl, map = map, lodcolumn = 1, add = TRUE, col = "#26A4C1")
quartz()
plot_scan1(x = G83_frac_ins_qtl, map = map, lodcolumn = 1, main = "G83 Fractional Ins Secretion | Stc1 Expression", col = "#FC81D0")
plot_scan1(x = G83_frac_ins_Stc1cov_qtl, map = map, lodcolumn = 1, add = TRUE, col = "#26A4C1")

# Stc1 scatter plots
Stc1_df <- cbind(pheno_data[,1:2], Stc1, G33_ins_sec_norm, G83_ins_sec_norm, G33_fract_ins_secrete, G83_fract_ins_secrete, Ins_per_islet)
quartz()
ggplot(data = Stc1_df, mapping = aes(x = Stc1, y = Ins_per_islet, color = sex)) + geom_point(size = 1.7) + 
  labs(x = "Stc1 Expression", y = "Ins per Islet") +
  ggtitle("Insulin per Islet vs. Stc1 Expression", subtitle = paste("Pearson's Correlation:", 
                                                                        round(cor(Stc1, Ins_per_islet, use = "complete.obs"), digits = 3), sep = " "))
quartz()
ggplot(data = Stc1_df, mapping = aes(x = Stc1, y = G83_fract_ins_secrete, color = sex)) + geom_point(size = 1.7) + 
  labs(x = "Stc1 Expression", y = "Fractional Insulin Secretion (8.3 mM Glucose Stimulation)") +
  ggtitle("Fractional Insulin Secretion vs. Stc1 Expression", subtitle = paste("Pearson's Correlation:", 
                                                                    round(cor(Stc1, G83_fract_ins_secrete, use = "complete.obs"), digits = 3), sep = " "))
quartz()
ggplot(data = Stc1_df, mapping = aes(x = Stc1, y = G83_ins_sec_norm, color = Ins_per_islet)) + geom_point(size = 1.7) + 
  labs(x = "Stc1 Expression", y = "Insulin Secretion (8.3 mM Glucose Stimulation)") +
  scale_color_gradient(low = "#0DED0D", high = "#0D5EED") +
  ggtitle("Insulin Secretion vs. Stc1 Expression", subtitle = paste("Pearson's Correlation:", 
                                                               round(cor(Stc1, G83_ins_sec_norm, use = "complete.obs"), digits = 3), sep = " "))

# Investigation of chr 1 peak on fractional insulin secretion traits
G33_frac_ins_qtl_chr1 <- G33_frac_ins_qtl[grep("1_", rownames(G33_frac_ins_qtl)),, drop = FALSE]
G33_frac_ins_qtl_chr1 <- G33_frac_ins_qtl_chr1[-grep("11_", rownames(G33_frac_ins_qtl_chr1)),, drop = FALSE]
rownames(G33_frac_ins_qtl_chr1)[which.max(G33_frac_ins_qtl_chr1)]

G83_frac_ins_qtl_chr1 <- G83_frac_ins_qtl[grep("1_", rownames(G83_frac_ins_qtl)),, drop = FALSE]
G83_frac_ins_qtl_chr1 <- G83_frac_ins_qtl_chr1[-grep("11_", rownames(G83_frac_ins_qtl_chr1)),, drop = FALSE]
rownames(G83_frac_ins_qtl_chr1)[which.max(G83_frac_ins_qtl_chr1)]

qtl1 <- genoprobs[,,dimnames(genoprobs)[[3]] == "1_32260090"]

# Batch Mediation of qtl14 on G33_fract_ins_secrete and G83_fract_ins_secrete
library(devtools)
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/mediation.scan.R")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/plot.mediation.R")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/gmb.coordinates.R")
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/mouse.chrlen.rda")

G33_ins_frac_mediationDF <- list(target = G33_fract_ins_secrete, mediator = rankz.mrna, 
                                 annotation = annot.mrna, covar = addcovar, qtl.geno = qtl14)

med <- mediation.scan(target = G33_ins_frac_mediationDF$target,
                      mediator = G33_ins_frac_mediationDF$mediator,
                      annotation = G33_ins_frac_mediationDF$annotation,
                      covar = G33_ins_frac_mediationDF$covar,
                      qtl.geno = G33_ins_frac_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Intermediate_Scripts/plot.mediation.chrcolors.R")
quartz()
plot.mediation_letmesleep(med)
quartz()
plot(med)   

G83_ins_frac_mediationDF <- list(target = G83_fract_ins_secrete, mediator = rankz.mrna, 
                                 annotation = annot.mrna, covar = addcovar, qtl.geno = qtl1)

med <- mediation.scan(target = G83_ins_frac_mediationDF$target,
                      mediator = G83_ins_frac_mediationDF$mediator,
                      annotation = G83_ins_frac_mediationDF$annotation,
                      covar = G83_ins_frac_mediationDF$covar,
                      qtl.geno = G83_ins_frac_mediationDF$qtl.geno)

# Plot mediation results and identify the mediator                      
quartz()
plot(med)

addcovar_G33 <- model.matrix(~Sex + G33_ins_sec_norm + Generation + diet_days, data=annot.samples)[,-1]
Ins_tAUC_covqtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                      kinship = K, addcovar = addcovar_G33, cores = 4, reml = TRUE)

quartz()
plot_scan1(x = Ins_tAUC_covqtl, map = map, lodcolumn = 1, col = "#FC81D0", main = "OGTT - tAUC | Ins. Sec. 3.3mM")
plot_scan1(x = Ins_tAUC_qtl, map = map, lodcolumn = 1,  add = TRUE, col = "#26A4C1")
