library(qtl2)
library(qtl2convert)

load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)
rownames(pheno_data) <- pheno_data[,1]

adiponectin_clean <- pheno_data[
  which(!is.na(pheno_data$adiponectin) == TRUE),
  75, drop = FALSE
]

annot.samples$diet_days <- pheno_data$diet_days
annot.samples_clean <- annot.samples[
  which(rownames(annot.samples) %in% rownames(adiponectin_clean)), , drop = FALSE]

genoprobs_clean <- genoprobs[which(rownames(genoprobs) %in% rownames(adiponectin_clean)), , ]
probs_clean <- probs_doqtl_to_qtl2(probs = genoprobs_clean, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")

K_clean = calc_kinship(probs = probs_clean, type = "loco", cores = 4)


addcovar_clean <- model.matrix(~Sex + Generation + diet_days, data= annot.samples_clean)[,-1]

qtl_clean <- scan1(genoprobs = probs_clean, pheno = adiponectin_clean[,1, drop = FALSE],
             kinship = K_clean, addcovar = addcovar_clean, cores = 4, reml=TRUE)

quartz()
plot_scan1(x = qtl_clean, map = map, lodcolumn = 1, main = colnames(qtl_clean)[1])

chr <- 16
qtl.coef = scan1coef(genoprobs = probs_clean[,chr], pheno = adiponectin_clean[,1,drop = FALSE],
                     kinship = K_clean[[chr]], addcovar = addcovar_clean)

quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl_clean, main = colnames(qtl_clean)[1])

qtl.blup <- scan1blup(genoprobs = probs_clean[,chr], pheno = adiponectin_clean[,1,drop = FALSE],
                      kinship = K_clean[[chr]], addcovar = addcovar_clean, cores = 4)

quartz()
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl_clean, main = colnames(qtl_clean)[1])



