load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
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

# creating scatterplots with every pair of clinical phenotype
pheno_allpheno_cor <- function(colnumber){
  for(i in 3:12){
    if(i > colnumber){
      x <- ggplot(mapping = aes(x = taki_phenos[,colnumber], y = taki_phenos[,i], color = taki_phenos$sex)) + geom_point() + 
        xlab(colnames(taki_phenos)[colnumber]) + ylab(colnames(taki_phenos)[i]) + 
        ggtitle(paste(colnames(taki_phenos)[i], colnames(taki_phenos)[colnumber], sep = " vs "), 
                subtitle = cor(taki_phenos[,colnumber], taki_phenos[,i], use = "complete.obs"))
      pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Initial_phenocor_plots/", 
                       paste(colnames(taki_phenos)[i], colnames(taki_phenos)[colnumber], sep = "_vs_"), ".pdf", sep = ""))
      plot(x)
      dev.off()
    }}}

for(i in 3:12){
  pheno_allpheno_cor(i)
}

# creating qtl maps for all of the clinical phenotypes

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

qtl1.1 <- scan1(genoprobs = probs, pheno = taki_phenos[,3:12,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Initial_phenocor_qtls/allphenos_qtl.pdf"), width = 18, height = 20)
par(mfrow = c(5,2))
for(i in 1:10){
  plot_scan1(x = qtl1.1, map = map, lodcolumn = i, main = colnames(qtl1.1)[i])
}
dev.off()

for(i in 1:10){
  pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Initial_phenocor_qtls/", 
                   paste(colnames(qtl1.1)[i], "_qtl", ".pdf", sep = ""), sep = ""), width = 15, height = 6)
  plot_scan1(x = qtl1.1, map = map, lodcolumn = i, main = colnames(qtl1.1)[i])
  dev.off()
}

# creating fasting insulin levels phenotype dataframe
Ins_6wk_norm <- log(pheno_data$Ins_6wk)
Ins_10wk_norm <- log(pheno_data$Ins_10wk)
Ins_14wk_norm <- log(pheno_data$Ins_14wk)
Ins_sac_norm <- log(pheno_data$Ins_sac)

insulin_data <- cbind(pheno_data[,1:2], 
                      Ins_6wk_norm, Ins_10wk_norm, Ins_14wk_norm, Ins_sac_norm)

# create qtl maps of the insulin phenotypes
qtl1.2 <- scan1(genoprobs = probs, pheno = insulin_data[,3:6,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Insulin_phenotypes_qtls/all_insulin_phenos_qtl.pdf"),width = 15, height = 15)
par(mfrow = c(4,1))
for(i in 1:4){
  plot_scan1(x = qtl1.2, map = map, lodcolumn = i, main = colnames(qtl1.2)[i])
}
dev.off()

for(i in 1:4){
  pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Insulin_phenotypes_qtls/", 
                   paste(colnames(qtl1.2)[i], "_qtl", ".pdf", sep = ""), sep = ""), width = 15, height = 6)
  plot_scan1(x = qtl1.2, map = map, lodcolumn = i, main = colnames(qtl1.2)[i])
  dev.off()
}


# playin around
Ins1 <- gene.expression.rankz("Ins1")
Ins2 <- gene.expression.rankz("Ins2")
Irs1 <- gene.expression.rankz("Irs1")
Irs2 <- gene.expression.rankz("Irs2")

Ins_Irs_dataframe <- cbind(pheno_data[,1:2], 
                           Ins1, Ins2, Irs1, Irs2)

qtl1.3 <- scan1(genoprobs = probs, pheno = Ins_Irs_dataframe[,3:6,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Ins_Irs_qtls/all_Ins_Irs_qtl.pdf"),width = 15, height = 15)
par(mfrow = c(4,1))
for(i in 1:4){
  plot_scan1(x = qtl1.3, map = map, lodcolumn = i, main = colnames(qtl1.3)[i])
}
dev.off()

for(i in 1:4){
  pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/Ins_Irs_qtls/", 
                   paste(colnames(qtl1.3)[i], "_qtl", ".pdf", sep = ""), sep = ""), width = 15, height = 6)
  plot_scan1(x = qtl1.3, map = map, lodcolumn = i, main = colnames(qtl1.3)[i])
  dev.off()
}

Gck <- gene.expression.rankz("Gck")

Glu_sac_norm <- log(pheno_data$Glu_sac)
Glu_6wk_norm <- log(pheno_data$Glu_6wk)
Glu_10wk_norm <- log(pheno_data$Glu_10wk)
Glu_14wk_norm <- log(pheno_data$Glu_14wk)

insulin_tAUC_data <- cbind(pheno_data[,1:2], 
                           Ins_0min_norm, Ins_tAUC_norm, Ins_sac_norm)

qtl2.1 <- scan1(genoprobs = probs, pheno = insulin_tAUC_data[,3:5,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

summary(qtl2.1)

which(qtl2.1[,1] == max(qtl2.1[,1]))
which(qtl2.1[,2] == max(qtl2.1[,2]))
which(qtl2.1[,3] == max(qtl2.1[,3]))

quartz()
par(mfrow = c(3,1))
for(i in 1:3){
  plot_scan1(x = qtl2.1, map = map, lodcolumn = i, main = colnames(qtl2.1)[i])
}


for(i in 1:3){
  quartz()
  chr <- 11
  qtl.coef = scan1coef(genoprobs = probs[,chr], pheno =  insulin_tAUC_data[,i+2,drop = FALSE],
                       kinship = K[[chr]], addcovar = addcovar)
  plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = qtl2.1, main = colnames(qtl2.1)[i])
}


for(i in 1:3){
  quartz()
  chr <- 11
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = insulin_tAUC_data[,i+2,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = qtl2.1, main = colnames(qtl2.1)[i])
}
