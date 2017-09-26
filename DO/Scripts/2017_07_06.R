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

# 9 phenotype p-value table cleanup 
Pvaluetable9pheno <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/Pvaluetable_9clinpheno.csv")

rowlist <- vector(mode = "numeric")
for(i in 1:21771){
  print(i)
  if(Pvaluetable9pheno[i,3] < 0.05 &
     Pvaluetable9pheno[i,4] < 0.05 &
     Pvaluetable9pheno[i,5] < 0.05 &
     Pvaluetable9pheno[i,6] < 0.05 &
     Pvaluetable9pheno[i,7] < 0.05 &
     Pvaluetable9pheno[i,8] < 0.05 &
     Pvaluetable9pheno[i,9] < 0.05 &
     Pvaluetable9pheno[i,10] < 0.05 &
     Pvaluetable9pheno[i,11] < 0.05){
    rowlist <- c(rowlist, i)
  }
}

Pvaluetable9pheno_significant <- Pvaluetable9pheno[rowlist, ]

# exportation of newly created csv file
write.csv(x=Pvaluetable9pheno_significant, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/Pvaluetable9pheno_significant.csv", row.names = FALSE)

# creating qtl mapping elements

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]
markers = read.delim("/Users/s-ishimt/Desktop/Jax_SSP'17/DO/DO_data/marker_grid_64K.txt")

# association mapping exploration
assoc_mapping <- function(probs, pheno, idx, addcovar, intcovar = NULL, k, markers, 
                         chr, start, end, ncores, 
                         db.file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/ccfoundersnps.sqlite") {
  
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

assoc.3 <- assoc_mapping(probs = probs, pheno = taki_phenos, idx = 3, 
                        addcovar = addcovar, k = K, markers = markers, 
                        chr = 11, start = 81, end = 86, ncores = 4)

# Association mapping plot with narrower boundaries.
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
             drop.hilit = 1)

# Create data structure for gene plot.
genes = data.frame(chr   = annot.mrna$chr,
                   start = annot.mrna$start,
                   stop  = annot.mrna$end,
                   strand = annot.mrna$strand,
                   Name   = annot.mrna$symbol, stringsAsFactors = F)


layout(matrix(1:2, 2, 1))
par(plt = c(0.1, 0.99, 0.1, 0.9))
plot_snpasso(scan1output = assoc.3[[1]], snpinfo = assoc.3[[2]], 
             drop.hilit = 1, xlim = c(81,86))

plot_genes(genes[genes$chr == 11 & genes$start > 81e6 & genes$stop < 86e6,], 
           xlim = c(81, 86))


# more narrowing down of the p-value genes

Pvaluetable9pheno_significant <- read.csv("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Derived_Data/Pvaluetable9pheno_significant.csv")
View(Pvaluetable9pheno_significant)

rowlist <- vector(mode = "numeric")
for(i in 1:327){
  if(i %in% order(Pvaluetable9pheno_significant[,3])[1:75] &
     i %in% order(Pvaluetable9pheno_significant[,4])[1:75] &
     i %in% order(Pvaluetable9pheno_significant[,5])[1:75] &
     i %in% order(Pvaluetable9pheno_significant[,7])[1:75] &
     i %in% order(Pvaluetable9pheno_significant[,8])[1:75]){
    rowlist <- c(rowlist, i)
  }
}

Pvaluetable9supersignificant <- Pvaluetable9pheno_significant[rowlist,]

genesofimportancelist <- as.character(Pvaluetable9supersignificant$annot.mrna.symbol)
genesofimportancedf <- rankz.mrna[,which(annot.mrna$symbol %in% genesofimportancelist)]
colnames(genesofimportancedf) <- genesofimportancelist

# qtl scans of the 12 genes of interest
qtlIns2 <- scan1(genoprobs = probs, pheno = genesofimportancedf[,1,drop = FALSE],
             kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlSsr4 <- scan1(genoprobs = probs, pheno = genesofimportancedf[,2,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlUbe2i <- scan1(genoprobs = probs, pheno = genesofimportancedf[,3,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlStard13 <- scan1(genoprobs = probs, pheno = genesofimportancedf[,4,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlGart <- scan1(genoprobs = probs, pheno = genesofimportancedf[,5,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlRad9a <- scan1(genoprobs = probs, pheno = genesofimportancedf[,6,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlPolr3a <- scan1(genoprobs = probs, pheno = genesofimportancedf[,7,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlFuom <- scan1(genoprobs = probs, pheno = genesofimportancedf[,8,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlKlf15 <- scan1(genoprobs = probs, pheno = genesofimportancedf[,9,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlAnp32a <- scan1(genoprobs = probs, pheno = genesofimportancedf[,10,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtlOs9 <- scan1(genoprobs = probs, pheno = genesofimportancedf[,11,drop = FALSE],
                 kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_06/11genes_qtl.pdf", width = 15, height = 20)
par(mfrow = c(6,2))
plot_scan1(x = qtlIns2, map = map, lodcolumn = 1, main = colnames(qtlIns2)[1])
plot_scan1(x = qtlSsr4, map = map, lodcolumn = 1, main = colnames(qtlSsr4)[1])
plot_scan1(x = qtlUbe2i, map = map, lodcolumn = 1, main = colnames(qtlUbe2i)[1])
plot_scan1(x = qtlStard13, map = map, lodcolumn = 1, main = colnames(qtlStard13)[1])
plot_scan1(x = qtlGart, map = map, lodcolumn = 1, main = colnames(qtlGart)[1])
plot_scan1(x = qtlRad9a, map = map, lodcolumn = 1, main = colnames(qtlRad9a)[1])
plot_scan1(x = qtlPolr3a, map = map, lodcolumn = 1, main = colnames(qtlPolr3a)[1])
plot_scan1(x = qtlFuom, map = map, lodcolumn = 1, main = colnames(qtlFuom)[1])
plot_scan1(x = qtlKlf15, map = map, lodcolumn = 1, main = colnames(qtlKlf15)[1])
plot_scan1(x = qtlAnp32a, map = map, lodcolumn = 1, main = colnames(qtlAnp32a)[1])
plot_scan1(x = qtlOs9, map = map, lodcolumn = 1, main = colnames(qtlOs9)[1])
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_06/11genes_qtl2.pdf", width = 20, height = 12)
plot_scan1(x = qtlIns2, map = map, lodcolumn = 1, main = colnames(qtlIns2)[1])
plot_scan1(x = qtlSsr4, map = map, lodcolumn = 1, main = colnames(qtlSsr4)[1])
plot_scan1(x = qtlUbe2i, map = map, lodcolumn = 1, main = colnames(qtlUbe2i)[1])
plot_scan1(x = qtlStard13, map = map, lodcolumn = 1, main = colnames(qtlStard13)[1])
plot_scan1(x = qtlGart, map = map, lodcolumn = 1, main = colnames(qtlGart)[1])
plot_scan1(x = qtlRad9a, map = map, lodcolumn = 1, main = colnames(qtlRad9a)[1])
plot_scan1(x = qtlPolr3a, map = map, lodcolumn = 1, main = colnames(qtlPolr3a)[1])
plot_scan1(x = qtlFuom, map = map, lodcolumn = 1, main = colnames(qtlFuom)[1])
plot_scan1(x = qtlKlf15, map = map, lodcolumn = 1, main = colnames(qtlKlf15)[1])
plot_scan1(x = qtlAnp32a, map = map, lodcolumn = 1, main = colnames(qtlAnp32a)[1])
plot_scan1(x = qtlOs9, map = map, lodcolumn = 1, main = colnames(qtlOs9)[1])
dev.off()

#qtl scans of the taki.phenos
qtl1.1 <- scan1(genoprobs = probs, pheno = taki_phenos[,3:11,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = paste("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_06/9phenotypes_qtl.pdf"), width = 18, height = 20)
par(mfrow = c(5,2))
for(i in 1:9){
  plot_scan1(x = qtl1.1, map = map, lodcolumn = i, main = colnames(qtl1.1)[i])
}
dev.off()

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_06/9phenotypes_qtl2.pdf", width = 20, height = 12)
for(i in 1:9){
  plot_scan1(x = qtl1.1, map = map, lodcolumn = i, main = colnames(qtl1.1)[i])
}
dev.off()

# qtl allele effects for chromosome 11 for taki.phenos
pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_06/9phenotypes_chr11.pdf", width = 20, height = 12)
for(i in 1:9){
  print(i)
  chr <- 11
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = taki_phenos[,i+2,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = qtl1.1[,i, drop = FALSE], main = colnames(qtl1.1)[i])
}
dev.off()

# qtl allele effects for chromosome 5 for 11 genes of interest
qtl2.1 <- scan1(genoprobs = probs, pheno = genesofimportancedf[,1:11,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_06/11genes_chr5.pdf", width = 20, height = 12)
for(i in 1:11){
  print(i)
  chr <- 5
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = genesofimportancedf[,i,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = qtl2.1[,i, drop = FALSE], main = colnames(qtl2.1)[i])
}
dev.off()



