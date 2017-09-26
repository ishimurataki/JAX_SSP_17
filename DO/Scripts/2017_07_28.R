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

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)
addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

qtl11 <- genoprobs[,,dimnames(genoprobs)[[3]] == "11_83591414"]

# association mapping with chr11 peak
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

assoc.3 <- assoc_mapping(probs = probs, pheno = ins_test_clinphenos, idx = 3, 
                         addcovar = addcovar, k = K, markers = markers, 
                         chr = 11, start = 82, end = 85, ncores = 4)
genes = data.frame(chr   = annot.mrna$chr,
                   start = annot.mrna$start,
                   stop  = annot.mrna$end,
                   strand = annot.mrna$strand,
                   Name   = annot.mrna$symbol, stringsAsFactors = F)
quartz()
layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
             drop.hilit = 1, xlim = c(82,85), cex = 0.75, ylim = c(0,9.95))

plot_genes(genes[genes$chr == 11 & genes$start > 82e6 & genes$stop < 85e6,], 
           xlim = c(82, 85))

# association mapping with only the snps that cause missense mutations
source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Missense_Association_function.R")

assoc.3 <- assoc_mapping_missense(probs = probs, pheno = ins_test_clinphenos, idx = 3, 
                                  addcovar = addcovar, k = K, markers = markers, 
                                  chr = 11, start = 82, end = 85, ncores = 4)
genes = data.frame(chr   = annot.mrna$chr,
                   start = annot.mrna$start,
                   stop  = annot.mrna$end,
                   strand = annot.mrna$strand,
                   Name   = annot.mrna$symbol, stringsAsFactors = F)
quartz()
layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
             drop.hilit = 1, xlim = c(82,85), cex = 0.8, ylim = c(0,9.95))

source("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Scripts/Snp_function_associationmapping.R")
scan1output = assoc.3[[1]]
snpinfo = assoc.3[[2]]
map <- snpinfo_to_map(snpinfo)
tmp <- expand_snp_results(scan1output, map, snpinfo)
scan1output <- tmp$lod
scan1output <- scan1output[which(rownames(scan1output) %in% snpinfo$snp),]

for(i in 1:179){
  if(scan1output[i] > 6){
    text(x = snpinfo$pos[snpinfo$snp == names(scan1output)[i]], y = scan1output[i], labels = names(scan1output)[i], pos = 4, srt = 270)
  }
}
plot_genes(genes[genes$chr == 11 & genes$start > 82e6 & genes$stop < 85e6,], 
           xlim = c(82, 85))


# Take control of the colors and plot the gene we mapped in red.
genes.ss = genes[genes$chr == 11 & genes$start > 82e6 & genes$stop < 85e6,]
# NOTE: You have to sort the genes by position for this to work.
genes.ss = genes.ss[order(genes.ss$start),]
gene.colors = rep("black", nrow(genes.ss))
gene.colors[genes.ss$Name %in% c("Hnf1b", "Synrg")] = "red"

quartz()
layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
             drop.hilit = 1, xlim = c(82,85), cex = 0.8, ylim = c(0,9.95))
plot_genes(genes.ss, xlim = c(82, 85), colors = gene.colors)


#######################################
# Creating Diabetes distribution Plot #
#######################################

Ins_sensitivity <- log(1/pheno_data$Ins_tAUC)
Ins_secretion <- G83_ins_sec_norm
Blood_Glucose <- log(pheno_data$Glu_tAUC)
Diabetes_spectrum_df <- cbind(pheno_data[,1:2, drop = FALSE], Ins_sensitivity,
                              Ins_secretion, Blood_Glucose)

quartz()
ggplot(data = Diabetes_spectrum_df, mapping = aes(x = Ins_secretion, y = Ins_sensitivity, col = log(pheno_data$Glu_tAUC))) + 
  geom_point(size = 3) + scale_color_gradient(low = "blue", high = "green") + 
  labs(title = "Insulin Secretiion vs Insulin Sensitivity vs Glucose Response", x = "Insulin Secretion", y = "Insulin Sensitivity", color = "Glucose Response")



