
# Test if a specific gene is an element of the dataset
gene.data.exists <- function(gene.name) {
  if (gene.name %in% annot$gene_symbol) {
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# Get gene expression data for a gene with a given name from a tissue
gene.exp <- function(gene.name, data.set) {
  if (gene.data.exists(gene.name)) {
  return(data.set[,annot$a_gene_id[which(annot$gene_symbol==gene.name)]])
  }
}

# Get clinical phenotype data from either phenotypes or phenotypes.rz
clinical <- function(clin.name, data.set) {
  if (clin.name %in% names(data.set)) {
    return(data.set[,clin.name])
  }
}

# get genotype given chromosome and cM position
genotype <- function(chr, pos) {
  if (chr > 0 && chr < 21 && pos >= 0 && pos <= 100) {
    return(f2g$geno[[chr]]$data[,find.marker(f2g, chr=chr, pos=pos)])
  }
}

quartz()
ggplot(aes(x=phenotypes.rz$insulin.10wk, y=Ucn3)) + geom_point(x=phenotypes.rz$insulin.10wk, y=Ucn3, color = phenotypes.rz$Weight)

for(i in 1:16677){
  if (cor(islet.rz[,i], Sst, use="complete.obs") < - 0.5 & cor(islet.rz[,i],phenotypes.rz$INS.10wk, use="complete.obs") > 0.5) {
    print("************************")
    print(subset(annot[i,], select="gene_symbol"))
    print(cor(islet.rz[,i], Sst, use="complete.obs"))
    print(cor(islet.rz[,i], phenotypes.rz$INS.10wk, use="complete.obs"))
  }
}

for(i in 1:16677){
  if (cor(islet.rz[,i], Sst, use="complete.obs") < - 0.5) {
    print("************************")
    print(subset(annot[i,], select="gene_symbol"))
    print(cor(islet.rz[,i], Sst, use="complete.obs"))
  }
}

for(i in 2:144){
  if (cor(phenotypes.rz[,i], Sst, use="complete.obs") < - 0.5) {
    print("************************")
    print(colnames(phenotypes.rz[i]))
    print(cor(phenotypes.rz[,i], Sst, use="complete.obs"))
  }
}

Sstr <- cbind(gene.exp("Sstr1", islet.rz), gene.exp("Sstr3", islet.rz), gene.exp("Sstr4", islet.rz), gene.exp("Sstr5", islet.rz))

Sstr[,5] <- rowSums(Sstr[,1:4], na.rm = TRUE)

Sst <- gene.exp("Sst", islet.rz)
Ghsr <- gene.exp("Ghsr", islet.rz)

INS.10wk <- phenotypes.rz$INS.10wk

f2g$pheno <- cbind(f2g$pheno[,c(2,6,7)], INS.10wk, Sst, Ghsr)

f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

Sst_islet_qtl <- scanone(f2g, pheno.col=5, method="hk")
Sst_islet_perm <- scanone(f2g, pheno.col=5, method="hk", n.perm=50, perm.Xsp = TRUE)

Ghsr_islet_qtl <- scanone(f2g, pheno.col=6, method="hk")
Ghsr_islet_perm <- scanone(f2g, pheno.col=6, method="hk", n.perm=50, perm.Xsp = TRUE)

quartz()
par(mfrow=c(1,2))
plot(Sst_islet_qtl,lodcolumn=1, main="QTL of Sst", ylab = "Sst", ylim=c(0,7))
add.threshold(Sst_islet_qtl, lodcolumn=1,
              perms = Sst_islet_perm, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(Sst_islet_qtl, lodcolumn=1,
              perms = Sst_islet_perm, alpha=0.20,lty="dashed",lwd=2,col="green")
plot(Ghsr_islet_qtl,lodcolumn=1, main="QTL of Ghsr", ylab = "Ghsr", ylim=c(0,7))
add.threshold(Ghsr_islet_qtl, lodcolumn=1,
              perms = Ghsr_islet_perm, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(Ghsr_islet_qtl, lodcolumn=1,
              perms = Ghsr_islet_perm, alpha=0.20,lty="dashed",lwd=2,col="green")

chr8qtl_Sst <- chr8qtl[which(!is.na(c(Sst,Ghsr)))]

chr8qlt2 <- chr8qtl[which(Ghsr[!is.na(Ghsr)==TRUE] | Sst[!is.na(Sst)==TRUE])]


