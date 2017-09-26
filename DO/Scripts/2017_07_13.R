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

# qtl scans of glucose stimulated insulin response phenotypes regressed with each other

addcovarG33 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G33_ins_sec_norm, data=annot.samples)[,-1]
addcovarG83 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G83_ins_sec_norm, data=annot.samples)[,-1]

G83_regressG33qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,9,drop = FALSE],
                       kinship = K, addcovar = addcovarG33, cores = 4, reml = TRUE)
G167_regressG33qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                           kinship = K, addcovar = addcovarG33, cores = 4, reml = TRUE)
G167_regressG83qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                       kinship = K, addcovar = addcovarG83, cores = 4, reml = TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/Insulin_Response_Regressed.pdf"), width = 18, height = 20)
par(mfrow = c(3,1))
plot_scan1(x = G83_regressG33qtl, map = map, lodcolumn = 1, main = "G83_ins_sec|G33_ins_sec")
plot_scan1(x = G167_regressG33qtl, map = map, lodcolumn = 1, main = "G167_ins_sec|G33_ins_sec")
plot_scan1(x = G167_regressG83qtl, map = map, lodcolumn = 1, main = "G167_ins_sec|G83_ins_sec")
dev.off()

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/Insulin_Response_Regressed2.pdf"), width = 18, height = 10)
plot_scan1(x = G83_regressG33qtl, map = map, lodcolumn = 1, main = "G83_ins_sec|G33_ins_sec")
plot_scan1(x = G167_regressG33qtl, map = map, lodcolumn = 1, main = "G167_ins_sec|G33_ins_sec")
plot_scan1(x = G167_regressG83qtl, map = map, lodcolumn = 1, main = "G167_ins_sec|G83_ins_sec")
dev.off()

# Chr 3,9, and 14 effect plots for 3 insulin response phenotypes
pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/chr3genes_AlleleEffects.pdf", width = 20, height = 12)
chr <- 3
qtl.blup1 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,9,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovarG33, cores = 4)
qtl.blup2 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG33, cores = 4)
qtl.blup3 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG83, cores = 4)
plot(x = qtl.blup1, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = G83_regressG33qtl, main = "G83_ins_sec|G33_ins_sec")
plot(x = qtl.blup2, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G167_regressG33qtl, main = "G167_ins_sec|G33_ins_sec")
plot(x = qtl.blup3, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G167_regressG83qtl, main = "G167_ins_sec|G83_ins_sec")
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/chr9genes_AlleleEffects.pdf", width = 20, height = 12)
chr <- 9
qtl.blup1 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,9,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG33, cores = 4)
qtl.blup2 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG33, cores = 4)
qtl.blup3 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG83, cores = 4)
plot(x = qtl.blup1, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G83_regressG33qtl, main = "G83_ins_sec|G33_ins_sec")
plot(x = qtl.blup2, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G167_regressG33qtl, main = "G167_ins_sec|G33_ins_sec")
plot(x = qtl.blup3, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G167_regressG83qtl, main = "G167_ins_sec|G83_ins_sec")
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_12/chr14genes_AlleleEffects.pdf", width = 20, height = 12)
chr <- 14
qtl.blup1 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,9,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG33, cores = 4)
qtl.blup2 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG33, cores = 4)
qtl.blup3 <- scan1blup(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovarG83, cores = 4)
plot(x = qtl.blup1, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G83_regressG33qtl, main = "G83_ins_sec|G33_ins_sec")
plot(x = qtl.blup2, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G167_regressG33qtl, main = "G167_ins_sec|G33_ins_sec")
plot(x = qtl.blup3, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = G167_regressG83qtl, main = "G167_ins_sec|G83_ins_sec")
dev.off()

##################################################################################################################
# Chromosome 3 peak on G83_Ins_sec|G33_ins_sec, G167_ins_sec|G33_ins_sec, and G167_ins_sec|G83_ins_sec mediation #
##################################################################################################################
G83_regressG33qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,9,drop = FALSE],
                           kinship = K, addcovar = addcovarG33, cores = 4, reml = TRUE)
G167_regressG33qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                            kinship = K, addcovar = addcovarG33, cores = 4, reml = TRUE)
G167_regressG83qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                            kinship = K, addcovar = addcovarG83, cores = 4, reml = TRUE)





