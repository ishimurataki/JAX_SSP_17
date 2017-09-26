install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

library(qtl2)
library(qtl2convert)

#importation of the DO data
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")

rownames(pheno_data) <- rownames(expr.mrna)
annot.samples$diet_days <- pheno_data$diet_days

pheno_data$coat_color <- replace(pheno_data$coat_color, which(pheno_data$coat_color=="agouti"), 0)
pheno_data$coat_color <- replace(pheno_data$coat_color, which(pheno_data$coat_color=="chinchilla"), 0)
pheno_data$coat_color <- replace(pheno_data$coat_color, which(pheno_data$coat_color=="WSB"), 0)
pheno_data$coat_color <- replace(pheno_data$coat_color, which(pheno_data$coat_color=="white"), 1)
pheno_data$coat_color <- replace(pheno_data$coat_color, which(pheno_data$coat_color=="black"), -1)

pheno_data$coat_color <- as.numeric(pheno_data$coat_color)

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
qtl1 <- scan1(genoprobs = probs, pheno = pheno_data[,5,drop = FALSE],
             kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

str(qtl)

################
# Plotting QTL #
################

# Linkage mapping plot.
quartz()
plot_scan1(x = qtl1, map = map, lodcolumn = 1, main = colnames(qtl1)[1])

# Calculate linear model coefficients on Chr 14.
chr <- 7
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno = pheno_data[,5,drop = FALSE],
                     kinship = K[[chr]], addcovar = addcovar)

# Plot coefficients with LOD cruve.
quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl1, main = colnames(qtl1)[1])

# Calculate BLUP coefficients on Chr 14.

qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = pheno_data[,5,drop = FALSE],
                      kinship = K[[chr]], addcovar = addcovar, cores = 4)

quartz()
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl, main = colnames(qtl)[1])

for(i in 1:157){
  if(class(pheno_data[,i]) == "numeric"){
    print(i)
    paste("qtl of", colnames(pheno_data)[i], sep = " ") <- scan1(genoprobs = probs, pheno = pheno_data[,i,drop = FALSE],
                                                                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
    if(max(paste("qtl of", colnames(pheno_data)[i], sep = " ")[,1]) < 8){
      rm(paste("qtl of", colnames(pheno_data)[i]))
    }
  }
}
  
for(i in c(1:9, 11:157)){
  if(class(pheno_data[,i]) == "numeric" | class(pheno_data[,i]) == "integer"){
    qtl <- scan1(genoprobs = probs, pheno = pheno_data[,i,drop = FALSE],
                                                                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
    if(max(qtl[,1]) < 8){
      rm(qtl)
    }
    else{
      assign(paste("qtl_of", colnames(pheno_data)[i], sep = "_"), qtl)
      print("******************")
      print(paste(colnames(pheno_data)[i], "success"), sep = " ")
      print(paste(max(qtl[,1]), "on chromosome", substr(rownames(qtl)[which(qtl[,1] == max(qtl[,1]))], 1,2), sep = " "))
    }
  }
}
