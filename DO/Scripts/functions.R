
# function to get gene expression given gene name from rankz.mrna data
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

gene.expression.raw <- function(gene.name){
  return(raw.mrna[,which(colnames(raw.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

    
# function to get clinical phenotype data given phenotype name
clin.pheno <- function(pheno.name){
  return(pheno_data[,which(colnames(pheno_data)==pheno.name)])
}

# typical ggplot formatting 
quartz()
ggplot(data=data, mapping = aes(y = y, x = x, group = group))+ geom_point() + geom_line()

# function to see correlation between gene expression and gene expresion
gene.gene.cor <- function(gene.name1, gene.name2){
  return(cor(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name1])], rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name2])], use = "complete.obs"))
}

# function to see correlation between gene expression and phenotype
gene.pheno.cor <- function(gene.name, pheno.name){
  return(cor(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])], log(pheno_data[,which(colnames(pheno_data)==pheno.name)]), use = "complete.obs"))
}

# function to see correlation between phenotype and phenotype
pheno.pheno.cor <- function(pheno.name1, pheno.name2){
  return(cor(log(pheno_data[,which(colnames(pheno_data)==pheno.name1)]), log(pheno_data[,which(colnames(pheno_data)==pheno.name2)]), use = "complete.obs"))
}

# function to compute pvalues adjusted for one covariate
pvalue_1cov <- function(x,f,c){
  fit1 <- lm(x ~ c+f)
  fit0 <- lm(x ~ c)
  anova(fit0, fit1)[2,6]
}

# function to insert row in a dataframe 
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

# function to retrieve genoprob data given chr number and position in bp
get_genoprob <- function(chr, position){
  colnumber1 <- which.min(abs(snps$pos[which(snps$chr == chr)] - position))
  chrsnps <- snps$pos[which(snps$chr == chr)]
  snppos <- chrsnps[colnumber1]
  colnumber2 <- which(dimnames(genoprobs)[[3]] == paste(chr, "_", snppos, sep = ""))
  return(genoprobs[,,colnumber2])
}


