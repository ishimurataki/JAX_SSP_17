#####################################################################
# Script for creating dataframe of p-values when clinical phenotype #
#   G83_ins_secrete is regressed for all gene expression data       #
#####################################################################
# load library
library(dplyr)

# log transform the G83_ins_secrete to normalize the data
G83_ins_secrete <- log(pheno_data$G83_ins_secrete)

# create a vector containing all p-values when G83_ins_secrete is modeled by all gene expression data
for(i in 1:21771){
  vector <- as.numeric(as.character(c(vector,anova(lm(G83_ins_secrete ~ pheno_data$sex + rankz.mrna[,i]))[2,5])))
}

vector <- vector[-1]

# create vector containing adjusted p-values 

vector2 <- p.adjust(vector, method = "BH")

# create dataframe 
df <- as.data.frame(cbind(annot.mrna$id, annot.mrna$symbol, vector, vector2))
colnames(df) <- c("id", "symbol", "p-value", "p-value.adjusted")
df$`p-value` <- as.numeric(as.character(df$`p-value`))
df$`p-value.adjusted` <- as.numeric(as.character(df$`p-value.adjusted`))
df <- arrange(df, df$`p-value`)

