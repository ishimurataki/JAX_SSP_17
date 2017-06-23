
# function to get gene expression given gene name from rankz.mrna data
gene.expression.rankz <- function(gene.name){
  return(rankz.mrna[,which(colnames(rankz.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

# function to get gene expression given gene name from raw.mrna data
gene.expression.raw <- function(gene.name){
  return(raw.mrna[,which(colnames(raw.mrna)==annot.mrna$id[annot.mrna$symbol == gene.name])])
}

# function to get clinical phenotype data given phenotype name
clin.pheno <- function(pheno.name){
  return(pheno_data[,which(colnames(pheno_data)==pheno.name)])
}

# see what clinical phenotypes are most correlated to G33_ins_secrete data
G33_ins_secrete <- log(pheno_data$G33_ins_secrete)

for(i in 1:157){
  if((class(pheno_data[,i])=="numeric") == TRUE){
    print(cor(pheno_data[,i], G33_ins_secrete, use = "complete.obs"))
    print(colnames(pheno_data)[i])
  }
}

# see what gene expression are most correlated to G33_ins_secrete data

for(i in 1:21771){
  if(abs(cor(G33_ins_secrete, rankz.mrna[,i], use = "complete.obs")) > 0.46) {
    print("*********************")
    print(annot.mrna$symbol[i])
    print(cor(G33_ins_secrete, rankz.mrna[,i], use = "complete.obs"))
  }
}

for(i in 1:157){
  if((class(pheno_data[,i])=="numeric") == TRUE ){
    print("****************")
    print(cor(pheno_data[,i], pheno_data$Ins_iAUC, use = "complete.obs"))
    print(colnames(pheno_data)[i])
  }
}

#create dataframe with insulin vs. time
insulin_data <- stack(pheno_data[,c(15,21,24,27)])

insulin_data[,2] <- as.character(insulin_data[,2])

colnames(insulin_data)[2] <- "num_week"

for(i in 1:1512){
  if(insulin_data[i,2] == "Ins_0min"){
    insulin_data[i,2] <- 0
  }
  if(insulin_data[i,2] == "Ins_6wk"){
    insulin_data[i,2] <- 6
  }
  if(insulin_data[i,2] == "Ins_10wk"){
    insulin_data[i,2] <- 10
  }
  if(insulin_data[i,2] == "Ins_14wk"){
    insulin_data[i,2] <- 14
  }
}

Mouse.ID <- as.character(pheno_data$Mouse.ID)
insulin_data <- cbind(Mouse.ID, insulin_data)

quartz()
ggplot(data=insulin_data, mapping = aes(y = values, x = num_week, group = Mouse.ID))+ geom_line()
