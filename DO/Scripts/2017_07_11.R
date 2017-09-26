install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl2)
library(qtl2convert)

load("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/DO378_islet.RData")
pheno_data <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", as.is = TRUE)
rownames(pheno_data) <- pheno_data[,1]

# creating improved dataframe:

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
  Ins_tAUC_lod_drop <- qtl1.1[which(rownames(qtl1.1) == Ins0minpeak), 1] - qtl[which(rownames(qtl) == Ins0minpeak), 1]
  Ins_0min_lod_drop <- qtl1.1[which(rownames(qtl1.1) == Ins0minpeak), 2] - qtl[which(rownames(qtl) == Ins0minpeak), 2]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak11intervalgenesdf)[i + 8], Ins_tAUC_lod_drop, Ins_0min_lod_drop)
  peak11geneLODdropDF <- rbind(peak11geneLODdropDF, geneLODrow)
}

# export dataframe into csv file
write.csv(x=peak11geneLODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_11/peak11geneLODdropDFv2.csv", row.names = FALSE)

peak11geneLODdropDFv2 <- peak11geneLODdropDFv2[order(-peak11geneLODdropDFv2$Ins_tAUC_lod_drop),]

###################################################################################################
# Chromosome 5 peak on Ins_tAUC_norm, Glu_tAUC_norm, Ins_0min_norm, and Weight_sac_norm mediation #
###################################################################################################

# qtl analysis on Ins_tAUC, Glu_tAUC, Ins_0min, weight_sac_norm
qtl1.1 <- scan1(genoprobs = probs, pheno = taki_phenos[,c(3:5,11),drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

# creating chromosome 5 subset of the qtl scan
chr5 <- qtl1.1[grep("5_", rownames(qtl1.1)), ]
chr5 <- chr5[- grep("15_", rownames(chr5)), ]

# subsetting chromosome 5 data to only include LOD scores that are between 135 - 150 mbp
chr5 <- chr5[which(as.numeric(substring(rownames(chr5), 3)) > 135*10^6), ]
chr5 <- chr5[which(as.numeric(substring(rownames(chr5), 3)) < 150*10^6), ]

# finding locations of max peaks in the subsetted chr 5 data
rownames(chr5)[which(chr5[,1] == max(chr5[,1]))]
rownames(chr5)[which(chr5[,2] == max(chr5[,2]))]
rownames(chr5)[which(chr5[,3] == max(chr5[,3]))]
rownames(chr5)[which(chr5[,4] == max(chr5[,4]))]

# Set Interval of interest to 145 - 150 mbp.
# Generate list of genes between 145 - 150 mbp using annot.mrna.
peak5intervalgenes <- annot.mrna[which(annot.mrna$chr == "5"), ]
peak5intervalgenes <- peak5intervalgenes[which(peak5intervalgenes$start > 145*10^6), ]
peak5intervalgenes <- peak5intervalgenes$symbol[which(peak5intervalgenes$end < 150*10^6)]

# generate mRNAseq dataframe for the list of genes in peak5intervalgenes using rankz.mrna
peak5intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak5intervalgenes)]
colnames(peak5intervalgenesdf) <- peak5intervalgenes

# combine peak11intervalgenesdataframe with annot.samples
peak5intervalgenesdf <- cbind(annot.samples, peak5intervalgenesdf)

# mediate chromosome 5 peak with gene expression of the genes that lie in the region of interest
Ins_tAUCpeak <- rownames(chr5)[which(chr5[,1] == max(chr5[,1]))]
Glu_tAUCpeak <- rownames(chr5)[which(chr5[,1] == max(chr5[,1]))]
Ins_0minpeak <- rownames(chr5)[which(chr5[,3] == max(chr5[,3]))]
weight_sacpeak <- rownames(chr5)[which(chr5[,4] == max(chr5[,4]))]

peak5geneLODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:63){
  print(i)
  addcovar_genecovar <- model.matrix(~Sex + Generation + diet_days + peak5intervalgenesdf[,i + 8], 
                                     data=peak5intervalgenesdf)[,-1]
  qtl <- scan1(genoprobs = probs, pheno = taki_phenos[,c(3:5,11),drop = FALSE],
               kinship = K, addcovar = addcovar_genecovar, cores = 4, reml=TRUE)
  Ins_tAUC_lod_drop <- qtl1.1[which(rownames(qtl1.1) == Ins_tAUCpeak), 1] - qtl[which(rownames(qtl) == Ins_tAUCpeak), 1]
  Glu_tAUC_lod_drop <- qtl1.1[which(rownames(qtl1.1) == Glu_tAUCpeak), 2] - qtl[which(rownames(qtl) == Glu_tAUCpeak), 2]
  Ins_0min_lod_drop <- qtl1.1[which(rownames(qtl1.1) == Ins_0minpeak), 3] - qtl[which(rownames(qtl) == Ins_0minpeak), 3]
  weight_sac_lod_drop <- qtl1.1[which(rownames(qtl1.1) == weight_sacpeak), 4] - qtl[which(rownames(qtl) == weight_sacpeak), 4]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak5intervalgenesdf)[i + 8], Ins_tAUC_lod_drop, Glu_tAUC_lod_drop, Ins_0min_lod_drop, weight_sac_lod_drop)
  peak5geneLODdropDF <- rbind(peak5geneLODdropDF, geneLODrow)
}

# export dataframe into csv file
write.csv(x=peak5geneLODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_11/peak5clinphenoLODdropDF.csv", row.names = FALSE)

##########################################################
# Chromosome 5 peak on various gene expression mediation #
##########################################################
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

Ube2i <- gene.expression.rankz("Ube2i")
Stard13 <- gene.expression.rankz("Stard13")
Gart <- gene.expression.rankz("Gart")
Rad9a <- gene.expression.rankz("Rad9a")
Os9<- gene.expression.rankz("Os9")

# creating dataframe with the selected gene expression
Genesdf <- cbind(pheno_data[,1:2], Ube2i, Stard13, Gart, Rad9a, Os9)

# qtl analysis on all gene expression
qtl2.1 <- scan1(genoprobs = probs, pheno = Genesdf[,c(3:7),drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

# creating chromosome 5 subset of the qtl scan
chr5 <- qtl2.1[grep("5_", rownames(qtl2.1)), ]
chr5 <- chr5[- grep("15_", rownames(chr5)), ]

# subsetting chromosome 5 data to only include LOD scores that are between 135 - 150 mbp
chr5 <- chr5[which(as.numeric(substring(rownames(chr5), 3)) > 135*10^6), ]
chr5 <- chr5[which(as.numeric(substring(rownames(chr5), 3)) < 149*10^6), ]

# finding locations of max peaks in the subsetted chr 5 data
rownames(chr5)[which(chr5[,1] == max(chr5[,1]))]
rownames(chr5)[which(chr5[,2] == max(chr5[,2]))]
rownames(chr5)[which(chr5[,3] == max(chr5[,3]))]
rownames(chr5)[which(chr5[,4] == max(chr5[,4]))]
rownames(chr5)[which(chr5[,5] == max(chr5[,5]))]

# Set Interval of interest to 145 - 150 mbp.
peak5intervalgenes <- annot.mrna[which(annot.mrna$chr == "5"), ]
peak5intervalgenes <- peak5intervalgenes[which(peak5intervalgenes$start > 145*10^6), ]
peak5intervalgenes <- peak5intervalgenes$symbol[which(peak5intervalgenes$end < 150*10^6)]

# generate mRNAseq dataframe for the list of genes in peak5intervalgenes using rankz.mrna
peak5intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak5intervalgenes)]
colnames(peak5intervalgenesdf) <- peak5intervalgenes

# combine peak11intervalgenesdataframe with annot.samples
peak5intervalgenesdf <- cbind(annot.samples, peak5intervalgenesdf)

# mediate chromosome 5 peak with gene expression of the genes that lie in the region of interest
peaklocation <- "5_146867112"

peak5geneLODdropDF <- data.frame(stringsAsFactors = FALSE)
for(i in 1:63){
  print(i)
  addcovar_genecovar <- model.matrix(~Sex + Generation + diet_days + peak5intervalgenesdf[,i + 8], 
                                     data=peak5intervalgenesdf)[,-1]
  qtl <- scan1(genoprobs = probs, pheno = Genesdf[,c(3:7),drop = FALSE],
               kinship = K, addcovar = addcovar_genecovar, cores = 4, reml=TRUE)
  Ube2i_lod_drop <- qtl2.1[which(rownames(qtl2.1) == peaklocation), 1] - qtl[which(rownames(qtl) == peaklocation), 1]
  Stard13_lod_drop <- qtl2.1[which(rownames(qtl2.1) == peaklocation), 2] - qtl[which(rownames(qtl) == peaklocation), 2]
  Gart_lod_drop <- qtl2.1[which(rownames(qtl2.1) == peaklocation), 3] - qtl[which(rownames(qtl) == peaklocation), 3]
  Rad9a_lod_drop <- qtl2.1[which(rownames(qtl2.1) == peaklocation), 4] - qtl[which(rownames(qtl) == peaklocation), 4]
  Os9_lod_drop <- qtl2.1[which(rownames(qtl2.1) == peaklocation), 5] - qtl[which(rownames(qtl) == peaklocation), 5]
  geneLODrow <- data.frame(Gene_symbol = colnames(peak5intervalgenesdf)[i + 8], Ube2i_lod_drop, Stard13_lod_drop, Gart_lod_drop, Rad9a_lod_drop, Os9_lod_drop)
  peak5geneLODdropDF <- rbind(peak5geneLODdropDF, geneLODrow)
}

# export dataframe into csv file
write.csv(x=peak5geneLODdropDF, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/2017_07_11/peak5genesLODdropDF.csv", row.names = FALSE)