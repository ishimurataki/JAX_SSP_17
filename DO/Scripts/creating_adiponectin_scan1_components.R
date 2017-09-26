# creating adiponectin phenotype
pheno <- pheno_data[, 75, drop = FALSE]
save(pheno, file="Adiponectin_Pheno.Rda")

# creating adiponectin probs
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
save(probs,file="Adiponectin_Probs.Rda")

# creating map
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
save(map, file = "Adiponectin_Map.Rda")

# creating adiponectin covariates
rownames(pheno_data) <- rownames(expr.mrna)
annot.samples$diet_days <- pheno_data$diet_days
addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]
save(addcovar, file = "Adiponectin_Covar.Rda")

# creating adiponectin kinship
kinship = calc_kinship(probs = probs, type = "loco", cores = 4)
save(kinship, file = "Adiponectin_Kinship.Rda")

save(pheno, probs, addcovar, kinship, file="Adiponectin_scan1_components.Rda")

rm(list=ls())
library(qtl2)
library(qtl2convert)

setwd("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Adiponectin")
load("Adiponectin_scan1_components.Rda")

qtl <- scan1(genoprobs = probs, pheno = pheno[,1, drop = FALSE],
             kinship = kinship, addcovar = addcovar, cores = 4, reml=TRUE)



