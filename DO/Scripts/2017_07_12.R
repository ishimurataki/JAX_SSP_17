load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)
rownames(pheno_data) <- pheno_data[,1]

# creating improved dataframe:

Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Glu_tAUC_norm <- log(pheno_data$Glu_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)
Glu_0min_norm <- log(pheno_data$Glu_0min)
oGTT_weight_norm <- log(pheno_data$oGTT_weight)
G33_ins_sec_norm <- log(pheno_data$G33_ins_secrete)
G83_ins_sec_norm <- log(pheno_data$G83_ins_secrete)
G167_ins_sec_norm <- log(pheno_data$G167_ins_secrete)
weight_sac_norm <- log(pheno_data$weight_sac)

# creating dataframe with the selected clinical phenotypes
taki_phenos <- cbind(pheno_data[,1:2], 
                     Ins_tAUC_norm, Glu_tAUC_norm, Ins_0min_norm, Glu_0min_norm, oGTT_weight_norm,
                     G33_ins_sec_norm, G83_ins_sec_norm, G167_ins_sec_norm, weight_sac_norm)
# creating qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# establishing genes of interest for chromosome 11 peak on Ins_tAUC and Ins_0min
E230016K23Rik <- gene.expression.rankz("E230016K23Rik")
A14Rik <- gene.expression.rankz("5530401A14Rik")
Slfn3 <- gene.expression.rankz("Slfn3")
Slfn2 <- gene.expression.rankz("Slfn2")

# genes of interest on chr 11 dataframe 
chr11genesdf <- cbind(pheno_data[,1:2], E230016K23Rik, A14Rik, Slfn3, Slfn2)

# qtl scans and allele effects for the chr 11 genes:
chr11genesqtl <- scan1(genoprobs = probs, pheno = chr11genesdf[,3:6,drop = FALSE],
                       kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/chr11genes_qtl.pdf"), width = 18, height = 10)
par(mfrow = c(2,2))
for(i in 1:4){
  plot_scan1(x = chr11genesqtl, map = map, lodcolumn = i, main = colnames(chr11genesqtl)[i])
}
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/chr11genes_AlleleEffects.pdf", width = 20, height = 12)
for(i in 1:4){
  print(i)
  chr <- 11
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = chr11genesdf[,i + 2,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = chr11genesqtl[,i, drop = FALSE], main = colnames(chr11genesqtl)[i])
}
dev.off()


# creating covariate dataframe including the selected genes on chromosome 11:
chr11covariates <- cbind(annot.samples, E230016K23Rik, A14Rik, Slfn3, Slfn2)

# Adding different gene expressions as covariates:
addcovar_baseline <- model.matrix(~Sex + Generation + diet_days, data=chr11covariates)[,-1]
addcovar_E230016K23Rik <- model.matrix(~Sex + Generation + diet_days + E230016K23Rik, data=chr11covariates)[,-1]
addcovar_A14Rik <- model.matrix(~Sex + Generation + diet_days + A14Rik, data=chr11covariates)[,-1]
addcovar_Slfn3 <- model.matrix(~Sex + Generation + diet_days + Slfn3, data=chr11covariates)[,-1]
addcovar_Slfn2 <- model.matrix(~Sex + Generation + diet_days + Slfn2, data=chr11covariates)[,-1]

# performing qtl scans with different covariates

A14Rikdf <- cbind(pheno_data[,1:2], A14Rik)
qtl1.1 <- scan1(genoprobs = probs, pheno = A14Rikdf[,3,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
plot_scan1(x = qtl1.1, map = map, lodcolumn = 1, chr = 11, main = colnames(qtl1.1)[1])


qtl11 <- genoprobs[,,which(dimnames(genoprobs)[[3]] == "11_81574923")]
anova(lm(Ins_tAUC_norm ~ E230016K23Rik + qtl11))

qtlhelp <-  scan1(genoprobs = probs, pheno = taki_phenos[,3:5,drop = FALSE],
                  kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

rownames(qtlhelp)[which(qtlhelp[,1] == max(qtlhelp[,1]))]
rownames(qtlhelp)[which(qtlhelp[,3] == max(qtlhelp[,3]))]


