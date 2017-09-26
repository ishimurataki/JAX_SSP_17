# GMLnet
# Dan Gatti
# August 4, 2017

library(qtl2convert)

# Load in the data.
load("DG_env.RData")
ls()

# Convert probs to qtl2 format.
genoprobs = probs_doqtl_to_qtl2(probs = genoprobs, map = snps, pos_column = "pos")

# Set the rownames in phenotype.
pheno = as.matrix(final_variables)
rownames(pheno) = rownames(add_covar)

# Fix the SNP positions to be in bp.
snps$bp = snps$bp * 1e6

# Phenotypes to map.
# alpha
# foodAvg
library(qtl2)

# Perform linkage mapping to verify that we get the same peak.
qtl.alpha = scan1(genoprobs = genoprobs, pheno = pheno[,"alpha", drop = F],
                  kinship = kin, addcovar = add_covar, cores = 4)

qtl.food = scan1(genoprobs = genoprobs, pheno = pheno[,"foodAvg", drop = F],
                 kinship = kin, addcovar = add_covar, cores = 4)

quartz()
layout(matrix(1:2, 2, 1))
plot_scan1(x = qtl.alpha, map = map, main = "alpha")
plot_scan1(x = qtl.food, map = map, main = "foodAvg")

layout(matrix(1:2, 2, 1))
plot_scan1(x = qtl.alpha, map = map, chr = 11, main = "alpha")
plot_scan1(x = qtl.food,  map = map, chr = 11, main = "foodAvg")

# Perform BLUP coef estimation on Chr 11.
chr = 11
blup.alpha = scan1coef(genoprobs = genoprobs[,chr], pheno = pheno[,"alpha", drop = F],
                       kinship = kin[[chr]], addcovar = add_covar)

blup.food = scan1coef(genoprobs = genoprobs[,chr], pheno = pheno[,"foodAvg", drop = F],
                      kinship = kin[[chr]], addcovar = add_covar)

layout(matrix(1:2, 2, 1))
plot_coefCC(x = blup.alpha, map = map)
plot_coefCC(x = blup.food,  map = map)

# Perform association mapping on proximal Chr 11 (10 - 25 Mb).
chr = 11
start = 10
end = 25
assoc.alpha = assoc_mapping(probs = genoprobs, pheno = final_variables, idx = 1,
                            addcovar = add_covar, k = kin, markers = snps,
                            chr = chr, start = start, end = end, ncores = 4)

assoc.food = assoc_mapping(probs = genoprobs, pheno = final_variables, idx = 6,
                           addcovar = add_covar, k = kin, markers = snps,
                           chr = chr, start = start, end = end, ncores = 4)

assoc.BW16 = assoc_mapping(probs = genoprobs, pheno = final_variables, idx = 5,
                           addcovar = add_covar, k = kin, markers = snps,
                           chr = chr, start = start, end = end, ncores = 4)

quartz()
layout(matrix(1:2, 2, 1))
plot_snpasso(scan1output = assoc.alpha[[1]], snpinfo = assoc.alpha[[2]], 
             drop.hilit = 1)
plot_snpasso(scan1output = assoc.food[[1]], snpinfo = assoc.food[[2]], 
             drop.hilit = 1)

quartz()
plot_snpasso(scan1output = assoc.BW16[[1]], snpinfo = assoc.BW16[[2]],
             drop.hilit = 1)


glm.alpha = assoc_glmnet(pheno = pheno[,"alpha", drop = F],
                         genoprobs = genoprobs, covar = add_covar, map = map, 
                         assoc = assoc.alpha[[1]], snpinfo = assoc.alpha[[2]], lod.thr = 1)

glm.food = assoc_glmnet(pheno = pheno[,"foodAvg", drop = F],
                        genoprobs = genoprobs, covar = add_covar, map = map, 
                        assoc = assoc.alpha[[1]], snpinfo = assoc.alpha[[2]], lod.thr = 1)

glm.BW16 = assoc_glmnet(pheno = pheno[,"BW16", drop = F],
                        genoprobs = genoprobs, covar = add_covar, map = map, 
                        assoc = assoc.alpha[[1]], snpinfo = assoc.alpha[[2]], lod.thr = 1)

alpha_scores <- glm.alpha[[1]]
food_scores <- glm.food[[1]]

nonzero_alpha <- alpha_scores[alpha_scores[,1] != 0,,drop = FALSE]
nonzero_food <- food_scores[food_scores[,1] != 0,,drop = FALSE]


alpha_rows <- rownames(nonzero_alpha)
food_rows <- rownames(nonzero_food)

nonzero_alpha <- as.data.frame(nonzero_alpha)
nonzero_food <- as.data.frame(nonzero_food)
colnames(nonzero_alpha) <- "score"
colnames(nonzero_food) <- "score"

nonzero_alpha <- abs(nonzero_alpha)
nonzero_food <- abs(nonzero_food)

nonzero_alpha <- nonzero_alpha[order(nonzero_alpha$score),, drop = FALSE]
nonzero_food <- nonzero_food[order(nonzero_food$score),, drop = FALSE]

alpha_snpinfo <- glm.alpha[[2]][alpha_rows,, drop = FALSE]
food_snpinfo <- glm.food[[2]][food_rows,, drop = FALSE]


layout(matrix(1:3, 3, 1))
plot_snpasso(scan1output = abs(glm.alpha[[1]]), snpinfo = glm.alpha[[2]], 
             cex = 1.25, ylab = "coef", main = "alpha")
plot_snpasso(scan1output = abs(glm.food[[1]]), snpinfo = glm.food[[2]],
             cex = 1.25, ylab = "coef", main = "foodAvg")
plot_snpasso(scan1output = abs(glm.BW16[[1]]), snpinfo = glm.BW16[[2]],
             cex = 1.25, ylab = "coef", main = "BW16")

find_peaks(scan1_output = qtl.alpha, map = map, threshold = 5, drop = 1.5)
find_peaks(scan1_output = qtl.BW16, map = map, threshold = 5, drop = 1.5)
find_peaks(scan1_output = qtl.foodAvg, map = map, threshold = 5, drop = 1.5)

nonzero_score <- numeric()
idx <- 1
for(i in scores) {
  if(i > 0) {
    nonzero_score <- c(nonzero_score, scores[idx])
  }
  idx <- idx + 1
}
