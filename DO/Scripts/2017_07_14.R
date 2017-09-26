##################################################################################################################
# Chromosome 3 peak on G83_Ins_sec|G33_ins_sec, G167_ins_sec|G33_ins_sec, and G167_ins_sec|G83_ins_sec mediation #
##################################################################################################################
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)
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

# creating dataframe with the selected clinical phenotypes
taki_phenos <- cbind(pheno_data[,1:2], 
                     Ins_tAUC_norm, Glu_tAUC_norm, Ins_0min_norm, Glu_0min_norm, oGTT_weight_norm,
                     G33_ins_sec_norm, G83_ins_sec_norm, G167_ins_sec_norm, weight_sac_norm)

# creating qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# Perform scans of G83_ins_sec and G167_ins_sec regressed with G33_ins_sec/G83_ins_sec
addcovarG33 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G33_ins_sec_norm, data=annot.samples)[,-1]
addcovarG83 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G83_ins_sec_norm, data=annot.samples)[,-1]

G83_regressG33qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,9,drop = FALSE],
                           kinship = K, addcovar = addcovarG33, cores = 4, reml = TRUE)
G167_regressG33qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                            kinship = K, addcovar = addcovarG33, cores = 4, reml = TRUE)
G167_regressG83qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                            kinship = K, addcovar = addcovarG83, cores = 4, reml = TRUE)

# Find confidence interval for chr 3 peak
G83_regressG33qtlchr3 <- G83_regressG33qtl[grep("3_", rownames(G83_regressG33qtl)), , drop = FALSE]
G83_regressG33qtlchr3 <- G83_regressG33qtlchr3[-grep("13_", rownames(G83_regressG33qtlchr3)), , drop = FALSE]
G167_regressG33qtlchr3 <- G167_regressG33qtl[grep("3_", rownames(G167_regressG33qtl)), , drop = FALSE]
G167_regressG33qtlchr3 <- G167_regressG33qtlchr3[-grep("13_", rownames(G167_regressG33qtlchr3)), , drop = FALSE]
G167_regressG83qtlchr3 <- G167_regressG83qtl[grep("3_", rownames(G167_regressG83qtl)), , drop = FALSE]
G167_regressG83qtlchr3 <- G167_regressG83qtlchr3[-grep("13_", rownames(G167_regressG83qtlchr3)), , drop = FALSE]

rownames(G83_regressG33qtlchr3)[which(G83_regressG33qtlchr3[,1] == max(G83_regressG33qtlchr3[,1]))]    #"3_138209955"
rownames(G167_regressG33qtlchr3)[which(G167_regressG33qtlchr3[,1] == max(G167_regressG33qtlchr3[,1]))] #"3_135485216"
rownames(G167_regressG83qtlchr3)[which(G167_regressG83qtlchr3[,1] == max(G167_regressG83qtlchr3[,1]))] #"3_135428123"

G83_regressG33qtlchr3[which(rownames(G83_regressG33qtlchr3)== "3_135485216"),1]

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
peak <- "3_135485216"
peak3LODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:36){
  print(i)
  addcovar_genecovar1 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G33_ins_sec_norm + peak3intervalgenesdf[,i + 8], 
                                     data=peak3intervalgenesdf)[,-1]
  addcovar_genecovar2 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G83_ins_sec_norm + peak3intervalgenesdf[,i + 8], 
                                     data=peak3intervalgenesdf)[,-1]
  qtl1 <- scan1(genoprobs = probs, pheno = taki_phenos[,9,drop = FALSE],
               kinship = K, addcovar = addcovar_genecovar1, cores = 4, reml=TRUE)
  qtl2 <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar1, cores = 4, reml=TRUE)
  qtl3 <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar2, cores = 4, reml=TRUE)
  G83regG33_LODdrop <- G83_regressG33qtl[which(rownames(G83_regressG33qtl) == peak), 1] - qtl1[which(rownames(qtl1) == peak), 1]
  G167regG33_LODdrop <- G167_regressG33qtl[which(rownames(G167_regressG33qtl) == peak), 1] - qtl2[which(rownames(qtl2) == peak), 1]
  G167regG83_LODdrop <- G167_regressG83qtl[which(rownames(G167_regressG83qtl) == peak), 1] - qtl3[which(rownames(qtl3) == peak), 1]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak3intervalgenesdf)[i + 8], G83regG33_LODdrop, G167regG33_LODdrop, G167regG83_LODdrop)
  peak3LODdropDF <- rbind(peak3LODdropDF, geneLODrow)
}

write.csv(x=peak3LODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_14/peak3LODdropDF.csv", row.names = FALSE)
##################################################################################################################
# Chromosome 9 peak on G83_Ins_sec|G33_ins_sec, G167_ins_sec|G33_ins_sec, and G167_ins_sec|G83_ins_sec mediation #
##################################################################################################################

# Find confidence interval for chr 9 peak
G167_regressG33qtlchr9 <- G167_regressG33qtl[grep("9_", rownames(G167_regressG33qtl)), , drop = FALSE]
G167_regressG33qtlchr9 <- G167_regressG33qtlchr9[-grep("19_", rownames(G167_regressG33qtlchr9)), , drop = FALSE]
G167_regressG83qtlchr9 <- G167_regressG83qtl[grep("9_", rownames(G167_regressG83qtl)), , drop = FALSE]
G167_regressG83qtlchr9 <- G167_regressG83qtlchr9[-grep("19_", rownames(G167_regressG83qtlchr9)), , drop = FALSE]
    
rownames(G167_regressG33qtlchr9)[which(G167_regressG33qtlchr9[,1] == max(G167_regressG33qtlchr9[,1]))] #"9_67475399"
rownames(G167_regressG83qtlchr9)[which(G167_regressG83qtlchr9[,1] == max(G167_regressG83qtlchr9[,1]))] #"9_67672764" # superior peak

G167_regressG33qtlchr9[which(rownames(G167_regressG33qtlchr9)== "9_67672764"),1]
G167_regressG33qtlchr9[which(rownames(G167_regressG33qtlchr9)== "9_67475399"),1] 

# set interval interest to 65 - 70 mbp.
# Generate list of genes between  65 - 70 mbp on chromosome 9 using annot.mrna.
peak9intervalgenes <- annot.mrna[which(annot.mrna$chr == "9"), ]
peak9intervalgenes <- peak9intervalgenes[which(peak9intervalgenes$start > 65*10^6), ]
peak9intervalgenes <- peak9intervalgenes$symbol[which(peak9intervalgenes$end < 70*10^6)]

# generate mRNAseq dataframe for the list of genes in peak9intervalgenes using rankz.mrna
peak9intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak9intervalgenes)]
colnames(peak9intervalgenesdf) <- peak9intervalgenes

# combine peak11intervalgenesdataframe with annot.samples
peak9intervalgenesdf <- cbind(annot.samples, peak9intervalgenesdf)

# mediate chromosome 9 peak with gene expression of the genes that lie in the region of interest
peak <- "9_67672764"
peak9LODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:48){
  print(i)
  addcovar_genecovar1 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G33_ins_sec_norm + peak9intervalgenesdf[,i + 8], 
                                      data=peak9intervalgenesdf)[,-1]
  addcovar_genecovar2 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G83_ins_sec_norm + peak9intervalgenesdf[,i + 8], 
                                      data=peak9intervalgenesdf)[,-1]
  qtl1 <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar1, cores = 4, reml=TRUE)
  qtl2 <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar2, cores = 4, reml=TRUE)
  G167regG33_LODdrop <- G167_regressG33qtl[which(rownames(G167_regressG33qtl) == peak), 1] - qtl1[which(rownames(qtl1) == peak), 1]
  G167regG83_LODdrop <- G167_regressG83qtl[which(rownames(G167_regressG83qtl) == peak), 1] - qtl2[which(rownames(qtl2) == peak), 1]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak9intervalgenesdf)[i + 8], G167regG33_LODdrop, G167regG83_LODdrop)
  peak9LODdropDF <- rbind(peak9LODdropDF, geneLODrow)
}

write.csv(x=peak9LODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_14/peak9LODdropDF.csv", row.names = FALSE)
##################################################################################################################
# Chromosome 14 peak on G83_Ins_sec|G33_ins_sec, G167_ins_sec|G33_ins_sec, and G167_ins_sec|G83_ins_sec mediation #
##################################################################################################################

# Find confidence interval for chr 14 peak
G83_regressG33qtlchr14 <- G83_regressG33qtl[grep("14_", rownames(G83_regressG33qtl)), , drop = FALSE]
G167_regressG33qtlchr14 <- G167_regressG33qtl[grep("14_", rownames(G167_regressG33qtl)), , drop = FALSE]
G167_regressG83qtlchr14 <- G167_regressG83qtl[grep("14_", rownames(G167_regressG83qtl)), , drop = FALSE]

rownames(G83_regressG33qtlchr14)[which(G83_regressG33qtlchr14[,1] == max(G83_regressG33qtlchr14[,1]))]    #"14_45939142" # superior peak
rownames(G167_regressG33qtlchr14)[which(G167_regressG33qtlchr14[,1] == max(G167_regressG33qtlchr14[,1]))] #"14_47338172"
rownames(G167_regressG83qtlchr14)[which(G167_regressG83qtlchr14[,1] == max(G167_regressG83qtlchr14[,1]))] #"14_47338172" 

G167_regressG33qtlchr14[which(rownames(G167_regressG33qtlchr14)== "14_45939142"),1]

# set interval interest to 44 - 49 mbp
# Generate list of genes between  44 - 49 mbp on chromosome 14 using annot.mrna.
peak14intervalgenes <- annot.mrna[which(annot.mrna$chr == "14"), ]
peak14intervalgenes <- peak14intervalgenes[which(peak14intervalgenes$start > 44*10^6), ]
peak14intervalgenes <- peak14intervalgenes$symbol[which(peak14intervalgenes$end < 49*10^6)]

# generate mRNAseq dataframe for the list of genes in peak14intervalgenes using rankz.mrna
peak14intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak14intervalgenes)]
colnames(peak14intervalgenesdf) <- peak14intervalgenes

# combine peak11intervalgenesdataframe with annot.samples
peak14intervalgenesdf <- cbind(annot.samples, peak14intervalgenesdf)

# mediate chromosome 14 peak with gene expression of the genes that lie in the region of interest
peak <- "14_45939142"
peak14LODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:42){
  print(i)
  addcovar_genecovar1 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G33_ins_sec_norm + peak14intervalgenesdf[,i + 8], 
                                      data=peak14intervalgenesdf)[,-1]
  addcovar_genecovar2 <- model.matrix(~Sex + Generation + diet_days + taki_phenos$G83_ins_sec_norm + peak14intervalgenesdf[,i + 8], 
                                      data=peak14intervalgenesdf)[,-1]
  qtl1 <- scan1(genoprobs = probs, pheno = taki_phenos[,9,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar1, cores = 4, reml=TRUE)
  qtl2 <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar1, cores = 4, reml=TRUE)
  qtl3 <- scan1(genoprobs = probs, pheno = taki_phenos[,10,drop = FALSE],
                kinship = K, addcovar = addcovar_genecovar2, cores = 4, reml=TRUE)
  G83regG33_LODdrop <- G83_regressG33qtl[which(rownames(G83_regressG33qtl) == peak), 1] - qtl1[which(rownames(qtl1) == peak), 1]
  G167regG33_LODdrop <- G167_regressG33qtl[which(rownames(G167_regressG33qtl) == peak), 1] - qtl2[which(rownames(qtl2) == peak), 1]
  G167regG83_LODdrop <- G167_regressG83qtl[which(rownames(G167_regressG83qtl) == peak), 1] - qtl3[which(rownames(qtl3) == peak), 1]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak14intervalgenesdf)[i + 8], G83regG33_LODdrop, G167regG33_LODdrop, G167regG83_LODdrop)
  peak14LODdropDF <- rbind(peak14LODdropDF, geneLODrow)
}

write.csv(x=peak14LODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_14/peak14LODdropDF.csv", row.names = FALSE)

