###################################################################
# Chromosome 11 peak on Ins_tAUC_norm and Ins_0min_norm mediation #
###################################################################

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

# qtl mapping of Ins_tAUC_norm and Ins_0min_norm for peak on chromosome 11
qtl1.1 <- scan1(genoprobs = probs, pheno = taki_phenos[,c(3,5),drop = FALSE],
                          kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
par(mfrow = c(2,1))
for(i in 1:2){
  plot_scan1(x = qtl1.1, map = map, lodcolumn = i, chr = 11, main = colnames(qtl1.1)[i])
}

# Set Interval of interest to 79 - 84 mbp.
# Generate list of genes between 79 - 84 mbp using annot.mrna.
peak11intervalgenes <- annot.mrna[which(annot.mrna$chr == "11"), ]
peak11intervalgenes <- peak11intervalgenes[which(peak11intervalgenes$start > 79*10^6), ]
peak11intervalgenes <- peak11intervalgenes$symbol[which(peak11intervalgenes$end < 84*10^6)]

# generate mRNAseq dataframe for the list of genes in peak11intervalgenes using rankz.mrna
peak11intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak11intervalgenes)]
colnames(peak11intervalgenesdf) <- peak11intervalgenes

# combine peak11intervalgenesdataframe with annot.samples
peak11intervalgenesdf <- cbind(annot.samples, peak11intervalgenesdf)


# mediate chromosome 11 peak with gene expression of the genes that lie in the region of interest
tAUCpeak <- rownames(qtl1.1)[which(qtl1.1[,1] == max(qtl1.1[grep("11_", rownames(qtl1.1)),1]))]
Ins0minpeak <- rownames(qtl1.1)[which(qtl1.1[,2] == max(qtl1.1[grep("11_", rownames(qtl1.1)),2]))]

peak11geneLODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:79){
  print(i)
  addcovar_genecovar <- model.matrix(~Sex + Generation + diet_days + peak11intervalgenesdf[,i + 8], 
                                     data=peak11intervalgenesdf)[,-1]
  qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,c(3,5),drop = FALSE],
               kinship = K, addcovar = addcovar_genecovar, cores = 4, reml=TRUE)
  Ins_tAUC_lod_drop <- qtl1.1[which(rownames(qtl1.1) == tAUCpeak), 1] - qtl[which(rownames(qtl) == tAUCpeak), 1]
  Ins_0min_lod_drop <- qtl1.1[which(rownames(qtl1.1) == tAUCpeak), 2] - qtl[which(rownames(qtl) == tAUCpeak), 2]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak11intervalgenesdf)[i + 8], Ins_tAUC_lod_drop, Ins_0min_lod_drop)
  peak11geneLODdropDF <- rbind(peak11geneLODdropDF, geneLODrow)
}

# export dataframe into csv file
write.csv(x=peak11geneLODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_11/peak11geneLODdropDF.csv", row.names = FALSE)




