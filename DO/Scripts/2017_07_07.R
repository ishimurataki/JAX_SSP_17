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
markers = read.delim("/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/marker_grid_64K.txt")

# qtl allele effects for chromosome 5 for 9 clinical phenotypes of interest
qtl3.1 <- scan1(genoprobs = probs, pheno = taki_phenos[,3:11,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

pdf(file = "/Users/s-ishimt/Desktop/Jax_SSP_17/DO/Plots/2017_07_07/9phenotypes_chr5.pdf", width = 20, height = 12)
for(i in 1:9){
  print(i)
  chr <- 5
  qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = taki_phenos[,i+2,drop = FALSE],
                        kinship = K[[chr]], addcovar = addcovar, cores = 4)
  plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
       scan1_output = qtl3.1[,i, drop = FALSE], main = colnames(qtl3.1)[i])
}
dev.off()

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
                         chr = 5, start = 81, end = 86, ncores = 4)

# searching for genes under chromosome 5 peak
chr5qtl3.1 <- qtl3.1[grep("5_", rownames(qtl3.1)),,drop = FALSE]
chr5qtl3.1 <- chr5qtl3.1[-grep("15_", rownames(chr5qtl3.1)),,drop = FALSE]
which(chr5qtl3.1[,5] == max(chr5qtl3.1[,5]))

chr5qtl2.1 <- qtl2.1[grep("5_", rownames(qtl2.1)),,drop = FALSE]
chr5qtl2.1 <- chr5qtl2.1[-grep("15_", rownames(chr5qtl2.1)),,drop = FALSE]
which(chr5qtl2.1[,2] == max(chr5qtl2.1[,2]))

# range: 145 - 149 mbp
chr5genesinterval <- annot.mrna[which(annot.mrna$chr == 5),]
chr5genesinterval <- chr5genesinterval[which(chr5genesinterval$middle_point > 145*10^6),]
chr5genesinterval <- chr5genesinterval$symbol[which(chr5genesinterval$middle_point < 149*10^6)]

chr5genesintervaldf <- rankz.mrna[,which(annot.mrna$symbol %in% chr5genesinterval)]
colnames(chr5genesintervaldf) <- chr5genesinterval

for(i in 1:49){
  print("**********************")
  print(cor(chr5genesintervaldf[,i], gene.expression.rankz("Stard13"), use = "complete.obs"))
  print(colnames(chr5genesintervaldf)[i])
}

