
pvalue_2cov <- function(x,f,c1, c2){
  fit1 <- lm(x ~ c1+ c2 + f)
  fit0 <- lm(x ~ c1 + c2)
  anova(fit0, fit1)[2,6]
}

vector <- vector('numeric')
for(i in 1:21771){
  print(i)
  vector <- c(vector, pvalue_2cov(clin_pheno11$G33_ins_secrete, rankz.mrna[,i], pheno_data$sex, gene.expression.rankz("Nfil3")))
}

vector <- p.adjust(vector, method = "BH")

df <- as.data.frame(cbind(annot.mrna$id, annot.mrna$symbol, vector))
colnames(df) <- c("id", "symbol", "pvalue.adjusted")
df$pvalue.adjusted <- as.numeric(as.character(df$pvalue.adjusted))
df <- arrange(df, df$pvalue.adjusted)

taki.cor <- function(x,y){
  
}

taki.sd <- function(x){
  avg <- sum(x)/length(x)
  vector <- vector('numeric')
  for(i in 1:length(x)){
    vector <- c(vector,(avg-x[i])^2)
  }
  standev <- sqrt(sum(vector)/(length(x)))
  return(standev)
}

taki.cor <- function(x,y){
  avg_x <- sum(x)/length(x)
  avg_y <- sum(y)/length(y)
  vector_x <- vector('numeric')
  vector_y <- vector('numeric')
  for(i in 1:length(x)){
    vector_x <- c(vector_x, (x[i] - avg_x)/taki.sd(x))
  }
  for(z in 1:length(y)){
    vector_y <- c(vector_y, (y[z] - avg_y)/taki.sd(y))
  }
  vector_xy <- vector('numeric')
  for(m in 1:length(x)){
    vector_xy <- c(vector_xy, vector_x[m]*vector_y[m])
  }
  correlation <- sum(vector_xy)/length(vector_xy)
  print(correlation)
}
