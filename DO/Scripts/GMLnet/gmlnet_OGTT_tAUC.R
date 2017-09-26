library(devtools)
library(qtl2)
library(qtl2convert)
library(ggplot2)

load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# Convert probs to qtl2 format.
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)
addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# Fix the SNP positions to be in bp.
snps$bp = snps$bp * 1e6

# establish phenotypes:
Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)

phenos <- cbind(pheno_data[,1:2, drop = FALSE], Ins_tAUC_norm, Ins_0min_norm)

# association mapping with only the snps that cause missense mutations
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Missense_Association_function.R")

# Perform association mapping on proximal Chr 11 (10 - 25 Mb).
chr = 11
start = 82
end = 85
markers = read.delim("/Users/s-ishimt/Desktop/marker_grid_64K.txt")
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/association_mapping_function.R")

assoc.tAUC <- assoc_mapping_missense(probs = probs, pheno = phenos, idx = 3, 
                                  addcovar = addcovar, k = K, markers = markers, 
                                  chr = chr, start = start, end = end, ncores = 4)

assoc.0min <- assoc_mapping_missense(probs = probs, pheno = phenos, idx = 4, 
                                     addcovar = addcovar, k = K, markers = markers, 
                                     chr = chr, start = start, end = end, ncores = 4)

assoc.tAUC <- assoc_mapping(probs = probs, pheno = phenos, idx = 3, 
                                     addcovar = addcovar, k = K, markers = markers, 
                                     chr = chr, start = start, end = end, ncores = 4)

assoc.0min <- assoc_mapping(probs = probs, pheno = phenos, idx = 4, 
                                     addcovar = addcovar, k = K, markers = markers, 
                                     chr = chr, start = start, end = end, ncores = 4)


quartz()
layout(matrix(1:2, 2, 1))
plot_snpasso(scan1output = assoc.tAUC[[1]], snpinfo = assoc.tAUC[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.0min[[1]], snpinfo = assoc.0min[[2]], 
             drop.hilit = 1)

# perfomring glmnet on OGTT-tAUC and Ins and 0min.
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/GMLnet/gmlnet_function.R")

glm.tAUC = assoc_glmnet(pheno = phenos[,3, drop = F],
                         genoprobs = probs, covar = addcovar, map = map, 
                         assoc = assoc.tAUC[[1]], snpinfo = assoc.tAUC[[2]], lod.thr = 1)

glm.0min = assoc_glmnet(pheno = phenos[,4, drop = F],
                        genoprobs = probs, covar = addcovar, map = map, 
                        assoc = assoc.0min[[1]], snpinfo = assoc.0min[[2]], lod.thr = 1)

tAUC_scores <- glm.tAUC[[1]]
zeromin_scores <- glm.0min[[1]]

nonzero_tAUC <- tAUC_scores[tAUC_scores[,1] != 0,,drop = FALSE]
nonzero_0min <- zeromin_scores[zeromin_scores[,1] != 0,,drop = FALSE]

tAUC_rows <- rownames(nonzero_tAUC)
zeromin_rows <- rownames(nonzero_0min)

nonzero_tAUC <- as.data.frame(nonzero_tAUC)
nonzero_0min <- as.data.frame(nonzero_0min)
colnames(nonzero_tAUC) <- "score"
colnames(nonzero_0min) <- "score"

nonzero_tAUC <- abs(nonzero_tAUC)
nonzero_0min <- abs(nonzero_0min)

nonzero_tAUC <- nonzero_tAUC[order(nonzero_tAUC$score),, drop = FALSE]
nonzero_0min <- nonzero_0min[order(nonzero_0min$score),, drop = FALSE]

tAUC_snpinfo <- glm.tAUC[[2]][tAUC_rows,, drop = FALSE]
zeromin_snpinfo <- glm.0min[[2]][zeromin_rows,, drop = FALSE]

quartz()
layout(matrix(1:2, 2, 1))
plot_snpasso(scan1output = abs(glm.tAUC[[1]]), snpinfo = glm.tAUC[[2]], 
             cex = 1.25, ylab = "coef", main = "OGTT - tAUC")
plot_snpasso(scan1output = abs(glm.0min[[1]]), snpinfo = glm.0min[[2]],
             cex = 1.25, ylab = "coef", main = "Ins 0min")

