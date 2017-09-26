load("DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# Establishing wanted clinical phenotypes:
Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Glu_tAUC_norm <- log(pheno_data$Glu_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)
Glu_0min_norm <- log(pheno_data$Glu_0min)
oGTT_weight_norm <- log(pheno_data$oGTT_weight)
G33_ins_sec_norm <- log(pheno_data$G33_ins_secrete)
G83_ins_sec_norm <- log(pheno_data$G83_ins_secrete)
G167_ins_sec_norm <- log(pheno_data$G167_ins_secrete)
weight_sac_norm <- log(pheno_data$weight_sac)
Ins_per_islet_norm <- log(pheno_data$Ins_per_islet)

# creating dataframe with the selected clinical phenotypes
taki_phenos <- cbind(pheno_data[,1:2], 
                     Ins_tAUC_norm, Glu_tAUC_norm, Ins_0min_norm, Glu_0min_norm, oGTT_weight_norm,
                     G33_ins_sec_norm, G83_ins_sec_norm, G167_ins_sec_norm, weight_sac_norm, Ins_per_islet_norm)

# creating qtl maps for all of the clinical phenotypes

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

qtl1.1<- scan1(genoprobs = probs, pheno = taki_phenos[,3,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.2<- scan1(genoprobs = probs, pheno = taki_phenos[,4,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.3<- scan1(genoprobs = probs, pheno = taki_phenos[,5,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.4<- scan1(genoprobs = probs, pheno = taki_phenos[,6,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.5<- scan1(genoprobs = probs, pheno = taki_phenos[,7,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.6<- scan1(genoprobs = probs, pheno = taki_phenos[,8,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.7<- scan1(genoprobs = probs, pheno = taki_phenos[,9,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.8<- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.9<- scan1(genoprobs = probs, pheno = taki_phenos[,11,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.10<- scan1(genoprobs = probs, pheno = taki_phenos[,12,drop = FALSE],
               kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

chr <- 11
qtl.coef1 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,3,drop = FALSE],
                     kinship = K[[chr]], addcovar = addcovar)
qtl.coef2 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,4,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef3 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,5,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef4 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,6,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef5 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,7,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef6 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,8,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef7 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,9,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef8 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef9 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,11,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef10 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,12,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)


quartz()
plot(x = qtl.coef1, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.1, main = colnames(qtl1.1)[1])
quartz()
plot(x = qtl.coef2, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.2, main = colnames(qtl1.2)[1])
quartz()
plot(x = qtl.coef3, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.3, main = colnames(qtl1.3)[1])
quartz()
plot(x = qtl.coef4, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.4, main = colnames(qtl1.4)[1])
quartz()
plot(x = qtl.coef5, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.5, main = colnames(qtl1.5)[1])
quartz()
plot(x = qtl.coef6, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.6, main = colnames(qtl1.6)[1])
quartz()
plot(x = qtl.coef7, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.7, main = colnames(qtl1.7)[1])
quartz()
plot(x = qtl.coef8, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.8, main = colnames(qtl1.8)[1])
quartz()
plot(x = qtl.coef9, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.9, main = colnames(qtl1.9)[1])
quartz()
plot(x = qtl.coef10, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.10, main = colnames(qtl1.10)[1])

chr <- 9
qtl.coef6 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,8,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef7 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,9,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
qtl.coef8 <- scan1coef(genoprobs = probs[,chr], pheno =  taki_phenos[,10,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)

quartz()
plot(x = qtl.coef6, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.6, main = colnames(qtl1.6)[1])
quartz()
plot(x = qtl.coef7, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.7, main = colnames(qtl1.7)[1])
quartz()
plot(x = qtl.coef8, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1.8, main = colnames(qtl1.8)[1])


