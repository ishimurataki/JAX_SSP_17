# This intent of this script is purely to explore the mediational tests that can be executed by
# qtl scans with different additive covariates. In this exploration, gene expression of Gcc1 will be regressed on Il6st 
# expression and the magnitude of the LOD drop will be calculated.

# loading the data
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# function to get gene expression given gene name from rankz.mrna data
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

# Creating the dataframe with Gcc1 and Il6st expression:
Gcc1 <- gene.expression.rankz("Gcc1")
Il6st <- gene.expression.rankz("Il6st")
Gcc1_Il6st_df <- cbind(pheno_data[,1:2], Gcc1, Il6st)

# Generating covariates dataframe
covariatesdf <- cbind(annot.samples, Il6st)

# creating qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=covariatesdf)[,-1]
addcovar_Il6st <- model.matrix(~Sex + Generation + diet_days + Il6st, data=covariatesdf)[,-1]

# Qtl mapping of Gcc1 expression with and without Il6st expression as a covariate
Gcc1_scan <- scan1(genoprobs = probs, pheno = Gcc1_Il6st_df[,3,drop = FALSE],
                  kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

Gcc1_Il6stcov_scan <- scan1(genoprobs = probs, pheno = Gcc1_Il6st_df[,3,drop = FALSE],
                   kinship = K, addcovar = addcovar_Il6st, cores = 4, reml=TRUE)

# plotting both qtl maps and visualizing the difference when Il6st expression is added as a covariate
quartz()
par(mfrow = c(2,1))
plot_scan1(x = Gcc1_scan, map = map, lodcolumn = 1, chr = 13, main = colnames(Gcc1_scan)[1])
plot_scan1(x = Gcc1_Il6stcov_scan, map = map, lodcolumn = 1, chr = 13, main = colnames(Gcc1_scan)[1])

peak <- rownames(Gcc1_scan)[which(Gcc1_scan[,1] == max(Gcc1_scan[grep("11_", rownames(Gcc1_scan)),1]))]
Gcc1_scan[which(rownames(Gcc1_scan) == peak), 1] - Gcc1_Il6stcov_scan[which(rownames(Gcc1_Il6stcov_scan) == peak), 1]
