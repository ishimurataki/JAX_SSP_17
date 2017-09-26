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

# confirming that the data is normal 
quartz()
par(mfrow = c(5,2))
for(i in 1:10){
  hist(ins_test_clinphenos[,i+2], main = colnames(ins_test_clinphenos)[i+2], breaks = 30)
}

# creating scatterplots with every pair of clinical phenotype
pheno_allpheno_cor <- function(colnumber){
  for(i in 3:12){
    if(i > colnumber){
      x <- ggplot(mapping = aes(x = ins_test_clinphenos[,colnumber], y = ins_test_clinphenos[,i], color = ins_test_clinphenos$sex)) + geom_point() + 
        xlab(colnames(ins_test_clinphenos)[colnumber]) + ylab(colnames(ins_test_clinphenos)[i]) + 
        ggtitle(paste(colnames(ins_test_clinphenos)[i], colnames(ins_test_clinphenos)[colnumber], sep = " vs "), 
                subtitle = cor(ins_test_clinphenos[,colnumber], ins_test_clinphenos[,i], use = "complete.obs"))
      pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_18/Iso_Insulin_pairwiseplots", 
                       paste(colnames(ins_test_clinphenos)[i], colnames(ins_test_clinphenos)[colnumber], sep = "_vs_"), ".pdf", sep = ""))
      plot(x)
      dev.off()
    }}}

for(i in 3:12){
  pheno_allpheno_cor(i)
}

# creating qtl mapping elements
probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

# performing qtl analysis
ins_phenos_qtl <- scan1(genoprobs = probs, pheno = ins_test_clinphenos[,3:12,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_18/Iso_nsulin_pheno_qtlmaps/allphenos_qtl.pdf"), width = 18, height = 20)
par(mfrow = c(5,2))
for(i in 1:10){
  plot_scan1(x = ins_phenos_qtl, map = map, lodcolumn = i, main = colnames(ins_phenos_qtl)[i])
}
dev.off()

for(i in 1:10){
  pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_18/Iso_nsulin_pheno_qtlmaps/", 
                   paste(colnames(ins_phenos_qtl)[i], "_qtl", ".pdf", sep = ""), sep = ""), width = 15, height = 6)
  plot_scan1(x = ins_phenos_qtl, map = map, lodcolumn = i, main = colnames(ins_phenos_qtl)[i])
  dev.off()
}

# association mapping with chr11 peak on Ins_tAUC and Ins_0min:
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
install.packages("RSQLite")
install.packages("dplyr")

library(qtl2)   # Loads qtl2geno, qtl2scan & qtl2plot.
library(qtl2convert)
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

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_18/Insulin_pheno_qtlmaps/chr11_tAUC_0min_associationmap.pdf", width = 18, height = 9)
for(i in 1:2){
  assoc.3 <- assoc_mapping(probs = probs, pheno = ins_test_clinphenos, idx = i+2, 
                           addcovar = addcovar, k = K, markers = markers, 
                           chr = 11, start = 80, end = 88, ncores = 4)

  # Create data structure for gene plot.
  genes = data.frame(chr   = annot.mrna$chr,
                     start = annot.mrna$start,
                     stop  = annot.mrna$end,
                     strand = annot.mrna$strand,
                     Name   = annot.mrna$symbol, stringsAsFactors = F)
  layout(matrix(1:2, 2, 1))
  par(plt = c(0.1, 0.99, 0.1, 0.9))
  plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
               drop.hilit = 1, xlim = c(80,88))
  
  plot_genes(genes[genes$chr == 11 & genes$start > 81e6 & genes$stop < 86e6,], 
             xlim = c(80, 88))
}
dev.off()

# qtl scans on all blood insulin levels at variable weeeks

Ins_6wk_norm <- log(pheno_data$Ins_6wk)
Ins_10wk_norm <- log(pheno_data$Ins_10wk)
Ins_14wk_norm <- log(pheno_data$Ins_14wk)
Ins_0min_norm <- log(pheno_data$Ins_0min)
Ins_sac_norm <- log(pheno_data$Ins_sac)

allinsulin_phenos <- cbind(pheno_data[,1:2], Ins_6wk_norm, Ins_10wk_norm, Ins_14wk_norm, Ins_0min_norm, Ins_sac_norm)

ins_phenos_qtl <- scan1(genoprobs = probs, pheno = allinsulin_phenos[,3:7,drop = FALSE],
                        kinship = K, addcovar = addcovar, cores = 4, reml = TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_18/insulin_weekly_qtlmaps/allphenos_qtl.pdf"), width = 18, height = 16)
par(mfrow = c(3,2))
for(i in 1:5){
  plot_scan1(x = ins_phenos_qtl, map = map, lodcolumn = i, main = colnames(ins_phenos_qtl)[i])
}
dev.off()

for(i in 1:5){
  pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_18/insulin_weekly_qtlmaps/", 
                   paste(colnames(ins_phenos_qtl)[i], "_qtl", ".pdf", sep = ""), sep = ""), width = 15, height = 6)
  plot_scan1(x = ins_phenos_qtl, map = map, lodcolumn = i, main = colnames(ins_phenos_qtl)[i])
  dev.off()
}


