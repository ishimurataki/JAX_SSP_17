#################################################################################
# Investigations of chr 3 peak on multiple clinical phenotype scans of interest #
#################################################################################

install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)
library(ggplot2)
load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)

rownames(pheno_data) <- pheno_data[,1]

# Establishing wanted clinical phenotypes:
Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)
Ins_0min_norm <- log(pheno_data$Ins_0min)
G33_ins_sec_norm <- log(pheno_data$G33_ins_secrete)
G83_ins_sec_norm <- log(pheno_data$G83_ins_secrete)
G167_ins_sec_norm <- log(pheno_data$G167_ins_secrete)
G83_FC_ins_secrete <- log(pheno_data$G83_FC_ins_secrete)
G167_FC_ins_secrete <- log(pheno_data$G167_FC_ins_secrete)
G33_fract_ins_secrete <- log(pheno_data$G33_fract_ins_secrete)
G83_fract_ins_secrete <- log(pheno_data$G83_fract_ins_secrete)
G167_fract_ins_secrete <- log(pheno_data$G167_fract_ins_secrete)

# creating dataframe with the selected clinical phenotypes
ins_test_clinphenos <- cbind(pheno_data[,1:2], 
                             Ins_tAUC_norm, Ins_0min_norm, G33_ins_sec_norm, 
                             G83_ins_sec_norm, G167_ins_sec_norm, G83_FC_ins_secrete, 
                             G167_FC_ins_secrete, G33_fract_ins_secrete, 
                             G83_fract_ins_secrete, G167_fract_ins_secrete)

# creating qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# qtl effect plots for chr 3
ins_phenos_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3:12,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_21/Iso_Insulin_chr3.pdf", width = 20, height = 12)
for(i in 1:10){
  print(i)
  chr <- 3
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = ins_test_clinphenos[,i+2,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = ins_phenos_qtl[,i, drop = FALSE], main = colnames(ins_phenos_qtl)[i])
}
dev.off()

#################################
# Chr 3 peak mediation analysis #
#################################

Ins_tAUC_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G167_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,7,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G83_FC_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,8,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G167_FC_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,9,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G167_frac_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,12,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

# Find confidence interval for chr 3 peak
Ins_tAUC_qtlchr3 <- Ins_tAUC_qtl[grep("3_", rownames(Ins_tAUC_qtl)), , drop = FALSE]
Ins_tAUC_qtlchr3 <- Ins_tAUC_qtlchr3[-grep("13_", rownames(Ins_tAUC_qtlchr3)), , drop = FALSE]
G167_ins_qtlchr3 <- G167_ins_qtl[grep("3_", rownames(G167_ins_qtl)), , drop = FALSE]
G167_ins_qtlchr3 <- G167_ins_qtlchr3[-grep("13_", rownames(G167_ins_qtlchr3)), , drop = FALSE]
G83_FC_ins_qtlchr3 <- G83_FC_ins_qtl[grep("3_", rownames(G83_FC_ins_qtl)), , drop = FALSE]
G83_FC_ins_qtlchr3 <- G83_FC_ins_qtlchr3[-grep("13_", rownames(G83_FC_ins_qtlchr3)), , drop = FALSE]
G167_FC_ins_qtlchr3 <- G167_FC_ins_qtl[grep("3_", rownames(G167_FC_ins_qtl)), , drop = FALSE]
G167_FC_ins_qtlchr3 <- G167_FC_ins_qtlchr3[-grep("13_", rownames(G167_FC_ins_qtlchr3)), , drop = FALSE]
G167_frac_ins_qtlchr3 <- G167_frac_ins_qtl[grep("3_", rownames(G167_frac_ins_qtl)), , drop = FALSE]
G167_frac_ins_qtlchr3 <- G167_frac_ins_qtlchr3[-grep("13_", rownames(G167_frac_ins_qtlchr3)), , drop = FALSE]

rownames(Ins_tAUC_qtlchr3)[which(Ins_tAUC_qtlchr3[,1] == max(Ins_tAUC_qtlchr3[,1]))]    #"3_135224922"
rownames(G167_ins_qtlchr3)[which(G167_ins_qtlchr3[,1] == max(G167_ins_qtlchr3[,1]))]    #"3_133549858"
rownames(G83_FC_ins_qtlchr3)[which(G83_FC_ins_qtlchr3[,1] == max(G83_FC_ins_qtlchr3[,1]))]    #"3_135647724"
rownames(G167_FC_ins_qtlchr3)[which(G167_FC_ins_qtlchr3[,1] == max(G167_FC_ins_qtlchr3[,1]))]    #"3_137348141"
rownames(G167_frac_ins_qtlchr3)[which(G167_frac_ins_qtlchr3[,1] == max(G167_frac_ins_qtlchr3[,1]))]    #"3_140779374

Ins_tAUC_qtlchr3[which(rownames(Ins_tAUC_qtlchr3)== "3_135570855"),1]
G167_ins_qtlchr3[which(rownames(G167_ins_qtlchr3)== "3_135570855"),1]
G83_FC_ins_qtlchr3[which(rownames(G83_FC_ins_qtlchr3)== "3_135570855"),1]
G167_FC_ins_qtlchr3[which(rownames(G167_FC_ins_qtlchr3)== "3_135570855"),1]
G167_frac_ins_qtlchr3[which(rownames(G167_frac_ins_qtlchr3)== "3_135570855"),1] # best peak to use is 3_135570855

# set interval interest to 133.5 - 138.5 mbp
# Generate list of genes between  133.5 - 138.5 mbp on chromosome 3 using annot.mrna.
peak3intervalgenes <- annot.mrna[which(annot.mrna$chr == "3"), ]
peak3intervalgenes <- peak3intervalgenes[which(peak3intervalgenes$start > 133.5*10^6), ]
peak3intervalgenes <- peak3intervalgenes$symbol[which(peak3intervalgenes$end < 138.5*10^6)]

# generate mRNAseq dataframe for the list of genes in peak3intervalgenes using rankz.mrna
peak3intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak3intervalgenes)]
colnames(peak3intervalgenesdf) <- peak3intervalgenes

# combine peak11intervalgenesdataframe with annot.samples
peak3intervalgenesdf <- cbind(annot.samples, peak3intervalgenesdf)

# mediate chromosome 3 peak with gene expression of the genes that lie in the region of interest
peak <- "3_135570855"
peak3LODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:36){
  print(i)
  addgenecovar <- model.matrix(~Sex + Generation + diet_days + peak3intervalgenesdf[,i + 8], 
                                      data=peak3intervalgenesdf)[,-1]
  Ins_tAUCcovar_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3,drop = FALSE],
                kinship = K, addcovar = addgenecovar, cores = 4, reml=TRUE)
  G167_inscovar_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,7,drop = FALSE],
                        kinship = K, addcovar = addgenecovar, cores = 4, reml = TRUE)
  G83_FCcovar_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,8,drop = FALSE],
                          kinship = K, addcovar = addgenecovar, cores = 4, reml = TRUE)
  G167_FCcovar_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,9,drop = FALSE],
                           kinship = K, addcovar = addgenecovar, cores = 4, reml = TRUE)
  G167_fraccovar_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,12,drop = FALSE],
                             kinship = K, addcovar = addgenecovar, cores = 4, reml = TRUE)
  Ins_tAUC_LODdrop <- Ins_tAUC_qtl[which(rownames(Ins_tAUC_qtl) == peak), 1] - Ins_tAUCcovar_qtl[which(rownames(Ins_tAUCcovar_qtl) == peak), 1]
  G167_ins_LODdrop <- G167_ins_qtl[which(rownames(G167_ins_qtl) == peak), 1] - G167_inscovar_qtl[which(rownames(G167_inscovar_qtl) == peak), 1]
  G83_FC_ins_LODdrop <- G83_FC_ins_qtl[which(rownames(G83_FC_ins_qtl) == peak), 1] - G83_FCcovar_ins_qtl[which(rownames(G83_FCcovar_ins_qtl) == peak), 1]
  G167_FC_ins_LODdrop <- G167_FC_ins_qtl[which(rownames(G167_FC_ins_qtl) == peak), 1] - G167_FCcovar_ins_qtl[which(rownames(G167_FCcovar_ins_qtl) == peak), 1]
  G167_frac_ins_LODdrop <- G167_frac_ins_qtl[which(rownames(G167_frac_ins_qtl) == peak), 1] - G167_fraccovar_ins_qtl[which(rownames(G167_fraccovar_ins_qtl) == peak), 1]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak3intervalgenesdf)[i + 8], Ins_tAUC_LODdrop, G83_FC_ins_LODdrop, G167_ins_LODdrop, G167_FC_ins_LODdrop, G167_frac_ins_LODdrop)
  peak3LODdropDF <- rbind(peak3LODdropDF, geneLODrow)
}

write.csv(x=peak3LODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_21/peak3LODdropDF.csv", row.names = FALSE)

