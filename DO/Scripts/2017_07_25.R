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

# association mapping with chr 3 plot
library(RSQLite)
library(dplyr)

markers = read.delim("/Users/s-ishimt/Desktop/marker_grid_64K.txt")
assoc_mapping <- function(probs, pheno, idx, addcovar, intcovar = NULL, k, markers, 
                          chr, start, end, ncores, 
                          db.file = "/Users/s-ishimt/Desktop/ccfoundersnps.sqlite") {
  
  # Subset probs and K to keep only the current chromosome.
  probs = probs[,chr]
  k     = k[[chr]]
  
  # Split up markers into a vector of map positions.
  map = split(markers[,3] * 1e-6, markers[,2])
  nm  = split(markers[,1], markers[,2])
  map = mapply(function(x, y) { names(x) = y;x }, map, nm)
  map = map[order(as.numeric(names(map)))]
  
  # Extract SNPs from the database
  my_db = src_sqlite(db.file, create = FALSE)
  snpinfo = tbl(my_db, sql(paste0("SELECT * FROM snps WHERE chr='", 
                                  chr, "' AND pos_Mbp>=", start, " AND pos_Mbp<=", end))) %>%
    collect(n = Inf)
  
  # Names have to be replaced for future methods
  colnames(snpinfo)[c(1,3)] = c("snp", "pos")
  
  # Index groups of similar SNPs.
  snpinfo = index_snps(map = map, snpinfo)
  
  # Find which phenotype data actually exist
  sel = !is.na(pheno[,idx])
  
  # Convert genoprobs to snpprobs.
  snppr = genoprob_to_snpprob(probs[sel,], snpinfo)
  
  # Scan1.
  assoc = scan1(pheno = pheno[sel,idx, drop = FALSE], kinship = k[sel,sel],
                genoprobs = snppr, addcovar = addcovar[sel,], 
                intcovar = addcovar[sel,intcovar], cores = ncores)
  
  # Return the scan data.
  return(list(assoc, snpinfo))
  
} # assoc_mapping()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_24/chr3_pheno_associationmap.pdf", width = 18, height = 9)
for(i in c(3,6,7,8,11)){
  print(i)
  assoc.3 <- assoc_mapping(probs = probs, pheno = ins_test_clinphenos, idx = i, 
                           addcovar = addcovar, k = K, markers = markers, 
                           chr = 3, start = 133, end = 139, ncores = 4)
  
  # Create data structure for gene plot.
  genes = data.frame(chr   = annot.mrna$chr,
                     start = annot.mrna$start,
                     stop  = annot.mrna$end,
                     strand = annot.mrna$strand,
                     Name   = annot.mrna$symbol, stringsAsFactors = F)
  layout(matrix(1:2, 2, 1))
  par(plt = c(0.1, 0.99, 0.1, 0.9))
  plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
               drop.hilit = 1, xlim = c(133,139))
  
  plot_genes(genes[genes$chr == 3 & genes$start > 133e6 & genes$stop < 139e6,], 
             xlim = c(133, 139))
}
dev.off()


# correlation of gene expression traits within region of interest on chr 3 and clinical phenotypes:
# set interval interest to 133.5 - 138.5 mbp
# Generate list of genes between  133.5 - 138.5 mbp on chromosome 3 using annot.mrna.
peak3intervalgenes <- annot.mrna[which(annot.mrna$chr == "3"), ]
peak3intervalgenes <- peak3intervalgenes[which(peak3intervalgenes$start > 133.5*10^6), ]
peak3intervalgenes <- peak3intervalgenes$symbol[which(peak3intervalgenes$end < 138.5*10^6)]

# generate mRNAseq dataframe for the list of genes in peak3intervalgenes using rankz.mrna
peak3intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak3intervalgenes)]
colnames(peak3intervalgenesdf) <- peak3intervalgenes

# create dataframe with correlations between gene expression and clin phenos
Ins_tAUCcor <- vector("numeric")
G167_inscor <- vector("numeric")
G83_FC_inscor <- vector("numeric")
G167_FC_inscor <- vector("numeric")
G167_fract_inscor <- vector("numeric")
for(i in 1:36){
  print(i)
  Ins_tAUCcor <- c(Ins_tAUCcor, cor(peak3intervalgenesdf[,i], Ins_tAUC_norm, use = "complete.obs"))
  G167_inscor <- c(G167_inscor, cor(peak3intervalgenesdf[,i], G167_ins_sec_norm, use = "complete.obs"))
  G83_FC_inscor <- c(G83_FC_inscor, cor(peak3intervalgenesdf[,i], G83_FC_ins_secrete, use = "complete.obs"))
  G167_FC_inscor <- c(G167_FC_inscor, cor(peak3intervalgenesdf[,i], G167_FC_ins_secrete, use = "complete.obs"))
  G167_fract_inscor <- c(G167_fract_inscor, cor(peak3intervalgenesdf[,i], G167_fract_ins_secrete, use = "complete.obs"))
}

cordata <- cbind(Ins_tAUCcor, G167_inscor, G83_FC_inscor, G167_FC_inscor, G167_fract_inscor)
colnames(cordata) <- c("Ins_tAUC", "G167_ins_sec", "G83_FC_ins_sec", "G167_FC_ins_sec", "G167_fract_ins_sec")
rownames(cordata) <- colnames(peak3intervalgenesdf)

##################################################################################
# Investigations of chr 14 peak on multiple clinical phenotype scans of interest #
##################################################################################
# qtl effect plots for chr 14
ins_test_clinphenos_2 <- ins_test_clinphenos[,c(3,5,10,11), drop = FALSE]
ins_phenos_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos_2[,,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_24/Iso_Insulin_chr14.pdf", width = 20, height = 12)
for(i in 1:4){
  print(i)
  chr <- 14
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = ins_test_clinphenos_2[,i,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = ins_phenos_qtl[,i, drop = FALSE], main = colnames(ins_phenos_qtl)[i])
}
dev.off()

#################################
# Chr 14 peak mediation analysis #
#################################
G33_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,5,drop = FALSE],
                      kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G33_fract_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,10,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)
G83_fract_ins_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,11,drop = FALSE],
                         kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

# Find confidence interval for chr 14 peak
G33_ins_qtlchr14 <- G33_ins_qtl[grep("14_", rownames(G33_ins_qtl)), , drop = FALSE]
G33_fract_ins_qtlchr14 <- G33_fract_ins_qtl[grep("14_", rownames(G33_fract_ins_qtl)), , drop = FALSE]
G83_fract_ins_qtlchr14 <- G83_fract_ins_qtl[grep("14_", rownames(G83_fract_ins_qtl)), , drop = FALSE]

rownames(G33_ins_qtlchr14)[which(G33_ins_qtlchr14[,1] == max(G33_ins_qtlchr14[,1]))]    #"14_91585596"
rownames(G33_fract_ins_qtlchr14)[which(G33_fract_ins_qtlchr14[,1] == max(G33_fract_ins_qtlchr14[,1]))]    #"14_91792449"
rownames(G83_fract_ins_qtlchr14)[which(G83_fract_ins_qtlchr14[,1] == max(G83_fract_ins_qtlchr14[,1]))]    #"14_91585596"
# set peak position at 14_91585596

# set interval interest to 89 - 94 mbp
# Generate list of genes between  89 - 94 mbp on chromosome 14 using annot.mrna.
peak14intervalgenes <- annot.mrna[which(annot.mrna$chr == "14"), ]
peak14intervalgenes <- peak14intervalgenes[which(peak14intervalgenes$start > 85*10^6), ]
peak14intervalgenes <- peak14intervalgenes$symbol[which(peak14intervalgenes$end < 100*10^6)]

# generate mRNAseq dataframe for the list of genes in peak14intervalgenes using rankz.mrna
peak14intervalgenesdf <- rankz.mrna[,which(annot.mrna$symbol %in% peak14intervalgenes)]
colnames(peak14intervalgenesdf) <- peak14intervalgenes

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

get_genoprob <- function(chr, position){
  colnumber1 <- which.min(abs(snps$pos[which(snps$chr == chr)] - position))
  chrsnps <- snps$pos[which(snps$chr == chr)]
  snppos <- chrsnps[colnumber1]
  colnumber2 <- which(dimnames(genoprobs)[[3]] == paste(chr, "_", snppos, sep = ""))
  return(genoprobs[,,colnumber2])
}

qtl14 <- get_genoprob(14, 91585596)

colnumber1 <- which.min(abs(snps$pos[which(snps$chr == 14)] - 91585596))
dimnames(genoprobs)[[3]][grep("14_",dimnames(genoprobs)[[3]])]

head(dimnames(genoprobs)[[3]][grep("^1_",dimnames(genoprobs)[[3]])])
head(snps$marker[grep("^1_", snps$marker)])

for(i in 1:19){
  print("******************")
  print(i)
  print(identical(dimnames(genoprobs)[[3]][grep(paste("^", i, "_", sep = "")
                                          ,dimnames(genoprobs)[[3]])],snps$marker[grep(paste("^", i, "_", sep = ""), snps$marker)]))
}

qtl14 <- genoprobs[,,dimnames(genoprobs)[[3]] == "14_91585596"]
