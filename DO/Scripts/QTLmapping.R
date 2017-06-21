library(tidyverse)
library(plotly)
load("~/data/DO378_islet.RData")

pheno <- read.csv("~/phenotypes/physiological_phenotypes.csv", as.is=TRUE)

# for mouse id < 100, add 0
pheno$mouse[nchar(pheno$mouse) == 5] <- sub("DO-", "DO-0", pheno$mouse[nchar(pheno$mouse) == 5])

names(pheno)[1] <- "Mouse.ID"
idx <- match(annot.samples$Mouse.ID, pheno$Mouse.ID)
pheno <- pheno[idx,]
stopifnot(annot.samples$Mouse.ID == pheno$Mouse.ID)

rownames(pheno) <- rownames(expr.mrna)

names(pheno)
for (i in 8:27) {
  hist(log10(pheno[,i]), main=names(pheno)[i])
}

pheno[,8:27] %>% 
  gather(var, value) %>%
  ggplot(aes(x=log10(value))) +
    geom_histogram() +
    facet_wrap(~var, scales="free")

ggsave("~/phenotypes/histograms.pdf")

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2plot)
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")


rownames(pheno) <- rownames(expr.mrna)
annot.samples$diet_days <- pheno$diet_days
covar2 <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)

rownames(pheno) <- pheno[,1]
fit2 <- scan1(genoprobs=probs, 
               kinship=Glist, 
               pheno=log10(pheno[,8:27]), 
               addcovar=covar2, 
               cores=10, 
               reml=TRUE)
maxlod <- apply(fit2$lod, 2, max)
maxlod[maxlod>7.18]

pdf("~/phenotypes/pheno_lodplots.pdf", width=10)
for (i in setdiff(1:ncol(fit2$lod),6)) { # problems with Ins_iAUC 
  plot(fit2, lodcolumn = i, main = colnames(pheno)[i+7], ylim=c(0,11))
  abline(h=7.18, col="red", lty=2)
}  
dev.off()

png(file = "~/phenotypes/lodplot_Ins_tAUC.png", width=8.06, height = 5.71, units = "in", res=300)
plot(fit2, lodcolumn = 5, main = "Ins_tAUC")
dev.off()

coef <- scan1coef(genoprobs=probs[,11], 
                  kinship=Glist[11], 
                  pheno=log10(pheno[,12,drop=FALSE]), 
                  addcovar=covar2, 
                  reml=TRUE)
plot_coefCC(coef)

annot <- annot.mrna
annot$chr[annot$chr=="MT"] <- "M"
names(annot)[names(annot) == "middle_point"] <- "pos"

qtl.pos = which.max(fit2$lod[,5])
qtl.geno = genoprobs[,,qtl.pos]
chr = modules$chr[i]

med <- mediation.scan(modexpr[,m], 
                      mediator = expr.mrna, 
                      annot, 
                      qtl.geno = qtl.geno, 
                      covar = as.matrix(covar))