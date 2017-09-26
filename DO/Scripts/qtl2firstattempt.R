##############################
# DO Data QTL2 June 22, 2017 #
##############################

# remove all preexisiting objects and variables in the environmnet
rm(list=ls())

library(tidyverse)
library(plotly)
library(qtl2)   # Loads qtl2geno, qtl2scan & qtl2plot.
library(qtl2convert)
library(RSQLite)
library(dplyr)

# set working directory to the physical location of the data
setwd("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data")

#importation of the DO data
load("DO378_islet.RData")
pheno_data <- read.csv("usable_clin_pheno_data.csv", as.is = TRUE)

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")

rownames(pheno_data) <- rownames(expr.mrna)
annot.samples$diet_days <- pheno_data$diet_days

###################
# Kinship matrix. #
###################
K = calc_kinship(probs = probs, type = "loco", cores = 4)

######################
# Set up covariates. #
######################
addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

#########
# Scan1 #
#########
qtl <- scan1(genoprobs = probs, pheno = log(pheno_data[,12,drop = FALSE]),
            kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

str(qtl)

################
# Plotting QTL #
################

# Linkage mapping plot.
quartz()
plot_scan1(x = qtl, map = map, lodcolumn = 1, main = colnames(qtl)[1])

# Calculate linear model coefficients on Chr 14.
chr <- 7
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno = log(pheno_data[,12,drop = FALSE]),
                     kinship = K[[chr]], addcovar = addcovar)

# Plot coefficients with LOD cruve.
quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl, main = colnames(qtl)[1])

# Calculate BLUP coefficients on Chr 14.

qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = log(pheno_data[,12,drop = FALSE]),
                     kinship = K[[chr]], addcovar = addcovar, cores = 4)

quartz()
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl, main = colnames(qtl)[1])

