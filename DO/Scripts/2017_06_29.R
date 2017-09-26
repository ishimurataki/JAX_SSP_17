Ins_tAUC_norm <- log(pheno_data$Ins_tAUC)

for(i in 1:157){
  if(class(pheno_data[,i]) == "numeric" | class(pheno_data[,i]) == "integer"){
    print("******************")
    print(colnames(pheno_data)[i])
    print(cor(log(pheno_data[,i]), Ins_tAUC_norm, use = "complete.obs"))
  }
}

# CLINICAL PHENOTYPES HIGHLY CORRELATED WITH Ins_tAUC_norm #
############################################################
# Ins_0min (0.7691675)
# oGTT_weight (0.6043955)
# G33/83/167_ins_secrete (0.45)
# Weight variables (0.6)
# HOMA_IR_0min (0.759)

cor(log(pheno_data$Ins_tAUC[
  which(pheno_data$sex == "F") 
]), log(pheno_data$oGTT_weight[
  which(pheno_data$sex == "F")
]), use = "complete.obs")

oGTT_weight_norm <- log(pheno_data$oGTT_weight)

oGTT_Ins_tAUC_plot <- ggplot(data = pheno_data, mapping = aes(x = log(oGTT_weight), y = log(Ins_tAUC), color = sex)) 
oGTT_Ins_tAUC_plot <- oGTT_Ins_tAUC_plot + geom_point() + geom_smooth(method = lm, se = FALSE, fullrange = TRUE)
oGTT_Ins_tAUC_plot <- oGTT_Ins_tAUC_plot + xlab("oGTT_weight") + ylab("Ins_tAUC") + ggtitle("Ins_tAUC vs. oGTT_weight")
quartz()
oGTT_Ins_tAUC_plot

anova(lm(Ins_tAUC_norm ~ oGTT_weight_norm + pheno_data$sex))
Ins_oGTT_lm <- lm(Ins_tAUC_norm ~ oGTT_weight_norm + pheno_data$sex)

Ins_oGTT_res <- Ins_oGTT_lm$residuals

for(i in 1:157){
  if(class(pheno_data[,i]) == "numeric" | class(pheno_data[,i]) == "integer"){
    print("******************")
    print(colnames(pheno_data)[i])
    print(cor(log(pheno_data[,i]), Ins_oGTT_res^0.57, use = "complete.obs"))
  }
}

for(i in 1:21771){
  if(abs(cor(Ins_tAUC_norm, rankz.mrna[,i], use ="complete.obs")) > 0.5){
    print("******************")
    print(annot.mrna$symbol[which(annot.mrna$id == colnames(rankz.mrna)[i])])
    print(cor(rankz.mrna[,i], Ins_tAUC_norm, use = "complete.obs"))
  }
}

G167_33_ratio <- log(pheno_data$G167_ins_secrete/pheno_data$G33_ins_secrete)
G167_83_ratio <- log(pheno_data$G167_ins_secrete/pheno_data$G83_ins_secrete)
G83_33_ratio <- log(pheno_data$G83_ins_secrete/pheno_data$G33_ins_secrete)

cor(G167_33_ratio,G167_83_ratio, use = "complete.obs")

cor(G83_33_ratio,Ins_tAUC_norm, use = "complete.obs")

Ins_Glu_stim_df <- cbind(pheno_data[,1:2], G167_33_ratio, G167_83_ratio, G83_33_ratio)
rownames(Ins_Glu_stim_df) <- pheno_data[,1]

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples)[,-1]

qtl1.1 <- scan1(genoprobs = probs, pheno = Ins_Glu_stim_df[,3,drop = FALSE],
              kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.2 <- scan1(genoprobs = probs, pheno = Ins_Glu_stim_df[,4,drop = FALSE],
              kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl1.3 <- scan1(genoprobs = probs, pheno = Ins_Glu_stim_df[,5,drop = FALSE],
              kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

quartz()
par(mfrow = c(3,1))
plot_scan1(x = qtl1.1, map = map, lodcolumn = 1, main = colnames(qtl1.1)[1])
plot_scan1(x = qtl1.2, map = map, lodcolumn = 1, main = colnames(qtl1.2)[1])
plot_scan1(x = qtl1.3, map = map, lodcolumn = 1, main = colnames(qtl1.3)[1])

quartz()
par(mfrow = c(3,1))
plot(log(pheno_data$G33_ins_secrete))
plot(log(pheno_data$G83_ins_secrete))
plot(log(pheno_data$G167_ins_secrete))

plot <- ggplot() + geom_point(mapping = aes(x = pheno_data$Mouse.ID, y = log(pheno_data$G33_ins_secrete), col = "red"))
plot <- plot + geom_point(mapping = aes(x = pheno_data$Mouse.ID, y = log(pheno_data$G83_ins_secrete), col = "orange"))
plot <- plot + geom_point(mapping = aes(x = pheno_data$Mouse.ID, y = log(pheno_data$G167_ins_secrete), col = "yellow"))
plot <- plot + scale_color_identity(guide = "legend", labels = c("G83_ins_secrete", "G33_ins_secrete", "G167_ins_secrete"))
plot(plot)

G33_ins_secrete_norm <- log(pheno_data$G33_ins_secrete)
G83_ins_secrete_norm <- log(pheno_data$G83_ins_secrete)
G167_ins_secrete_norm <- log(pheno_data$G167_ins_secrete)

df1 <- cbind(pheno_data[,1:2], G33_ins_secrete_norm, G83_ins_secrete_norm, G167_ins_secrete_norm, Ins_tAUC_norm)
rownames(df1) <- df1$Mouse.ID

qtl2.1 <- scan1(genoprobs = probs, pheno = df1[,3,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl2.2 <- scan1(genoprobs = probs, pheno = df1[,4,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl2.3 <- scan1(genoprobs = probs, pheno = df1[,5,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

quartz()
par(mfrow = c(3,1))
plot_scan1(x = qtl2.1, map = map, lodcolumn = 1, main = colnames(qtl2.1)[1])
plot_scan1(x = qtl2.2, map = map, lodcolumn = 1, main = colnames(qtl2.2)[1])
plot_scan1(x = qtl2.3, map = map, lodcolumn = 1, main = colnames(qtl2.3)[1])

rownames(pheno_data) <- pheno_data[,1]

qtl3.1 <- scan1(genoprobs = probs, pheno = log(pheno_data[,12,drop = FALSE]),
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
plot_scan1(x = qtl3.1, map = map, lodcolumn = 1, main = colnames(qtl3.1)[1])

chr <- 5
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno = log(pheno_data[,12,drop = FALSE]),
                     kinship = K[[chr]], addcovar = addcovar)

quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl3.1, main = colnames(qtl3.1)[1])

df2 <- df1[order(df1$G33_ins_secrete_norm), ]
df2 <- cbind(df1, c(1:nrow(df1)))
colnames(df2)[7] <- "MICE"

plot1 <- ggplot(data = df2, mapping = aes(x = MICE, y = G33_ins_secrete_norm, col = "green")) + geom_point()
plot1 <- plot1 + geom_point(data = df2, mapping = aes(x = MICE, y = G83_ins_secrete_norm, col = "blue"))
plot1 <- plot1 + geom_point(data = df2, mapping = aes(x = MICE, y = G167_ins_secrete_norm, col = "purple"))
plot1 <- plot1 + scale_color_identity(guide = "legend", labels = c("G83_ins_secrete", "G33_ins_secrete", "G167_ins_secrete"))

df3 <- df1[order(df1$G83_ins_secrete_norm), ]
df3 <- cbind(df3, c(1:nrow(df1)))
colnames(df3)[7] <- "MICE"

plot2 <- ggplot(data = df3, mapping = aes(x = MICE, y = G83_ins_secrete_norm, col = "green")) + geom_point()
plot2 <- plot2 + geom_point(data = df3, mapping = aes(x = MICE, y = G33_ins_secrete_norm, col = "blue"))
plot2 <- plot2 + geom_point(data = df3, mapping = aes(x = MICE, y = G167_ins_secrete_norm, col = "purple"))
plot2 <- plot2 + scale_color_identity(guide = "legend", labels = c("G33_ins_secrete", "G83_ins_secrete", "G167_ins_secrete"))

plot(Ins_tAUC_norm, G33_ins_secrete_norm, use = "complete.obs")

ggplot(data = df1, mapping = aes(x = G33_ins_secrete_norm, y = Ins_tAUC_norm, color = pheno_data$weight_sac, group = sex)) + geom_point()

tAUC_G33_lm <- lm(Ins_tAUC_norm ~ G33_ins_secrete_norm)
tAUC_G33_res <- tAUC_G33_lm$residuals

df4 <- cbind(pheno_data, tAUC_G33_res)

for(i in 1:157){
  if(class(pheno_data[,i]) == "numeric" | class(pheno_data[,i]) == "integer"){
    print("******************")
    print(colnames(pheno_data)[i])
    print(cor(log(pheno_data[,i]), df4$tAUC_G33_res, use = "complete.obs"))
  }
}

plot(log(pheno_data$weight_sac), G83_ins_secrete_norm, xlab = "sacrifice weight", ylab = "G83 Insulin Sec", main = "G83 Insulin Sec vs. Sacrifice Weight")
abline(lm(G83_ins_secrete_norm ~ log(pheno_data$weight_sac)), col = "black")

G83_sacweight_lm <- lm(G83_ins_secrete_norm ~ log(pheno_data$weight_sac), na.action = na.exclude)
G83_sacweight_res <- as.data.frame(G83_sacweight_lm$residuals)

which(c(1:378) %in%  rownames(G83_sacweight_res) == FALSE)

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

G83_sacweight_res <- insertRow(G83_sacweight_res, newrow, 66)
G83_sacweight_res <- insertRow(G83_sacweight_res, newrow, 67)
G83_sacweight_res <- insertRow(G83_sacweight_res, newrow, 68)

G83_sacweight_res <- G83_sacweight_res[,1]

hist(log(pheno_data$Ins_iAUC - min(pheno_data$Ins_iAUC, na.rm = TRUE)), breaks = 30)
hist(pheno_data$Ins_iAUC, breaks = 30)

df5 <- cbind(pheno_data[,1:2], G83_sacweight_res)

qtl4.1 <- scan1(genoprobs = probs, pheno = df5[,3,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)

quartz()
plot_scan1(x = qtl4.1, map = map, lodcolumn = 1, main = colnames(qtl4.1)[1])

chr <- 19
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno = df5[,3,drop = FALSE],
                     kinship = K[[chr]], addcovar = addcovar)

quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl4.1, main = colnames(qtl4.1)[1])

for(i in 1:21771){
  if(abs(cor(G83_sacweight_res, rankz.mrna[,i], use ="complete.obs")) > 0.5){
    print("******************")
    print(annot.mrna$symbol[which(annot.mrna$id == colnames(rankz.mrna)[i])])
    print(cor(rankz.mrna[,i], G83_sacweight_res, use = "complete.obs"))
  }
}

Nfil3 <- gene.expression.rankz("Nfil3")
Dnajc6 <- gene.expression.rankz("Dnajc6")
Gcc1 <- gene.expression.rankz("Gcc1")
Eif1 <- gene.expression.rankz("Eif1")
Siah2 <- gene.expression.rankz("Siah2")

weight_sac_norm <- log(pheno_data$weight_sac)
WPIC_norm <- log(pheno_data$WPIC)
Ins_per_islet_norm <- log(pheno_data$Ins_per_islet)

df6 <- cbind(pheno_data[,1:2], G33_ins_secrete_norm, G83_ins_secrete_norm, G167_ins_secrete_norm, weight_sac_norm, WPIC_norm, Ins_per_islet_norm, Nfil3, Dnajc6, Gcc1, Eif1, Siah2, residuals, Ins_tAUC_norm)

ggplot(data = df6, mapping = aes(x = weight_sac_norm, y = G83_ins_secrete_norm, color = Nfil3)) + geom_point(size = 5) +scale_color_gradient(low = "blue", high = "green")

Ins_iAUC_norm <- (pheno_data$Ins_iAUC - min(pheno_data$Ins_iAUC, na.rm = TRUE))^3

annot.samples2 <- cbind(annot.samples, weight_sac_norm, Ins_per_islet_norm)

addcovar <- model.matrix(~Sex + Generation + diet_days + weight_sac_norm + Ins_per_islet_norm, data=annot.samples2)[,-1]

qtl5.1 <- scan1(genoprobs = probs, pheno = df6[,3,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl5.2 <- scan1(genoprobs = probs, pheno = df6[,4,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
qtl5.3 <- scan1(genoprobs = probs, pheno = df6[,5,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
par(mfrow = c(3,1))
plot_scan1(x = qtl5.1, map = map, lodcolumn = 1, main = colnames(qtl5.1)[1])
plot_scan1(x = qtl5.2, map = map, lodcolumn = 1, main = colnames(qtl5.2)[1])
plot_scan1(x = qtl5.3, map = map, lodcolumn = 1, main = colnames(qtl5.3)[1])

anova(lm(G83_ins_secrete_norm ~ annot.samples2$Sex + annot.samples2$Generation + annot.samples2$diet_days + Ins_per_islet_norm + weight_sac_norm))

equation1 <- lm(G83_ins_secrete_norm ~ annot.samples2$Sex + annot.samples2$Generation + annot.samples2$diet_days + Ins_per_islet_norm + weight_sac_norm)
residuals <- as.data.frame(equation1$residuals)

residuals <- insertRow(residuals, newrow, 66)
residuals <- insertRow(residuals, newrow, 67)
residuals <- insertRow(residuals, newrow, 68)
residuals <- residuals[,1]

qtl6.1 <- scan1(genoprobs = probs, pheno = df6[,14,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
plot_scan1(x = qtl6.1, map = map, lodcolumn = 1, main = colnames(qtl6.1)[1])

addcovar <- model.matrix(~Sex + Generation + diet_days, data=annot.samples2)[,-1]

qtl6.2 <- scan1(genoprobs = probs, pheno = df6[,8,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
plot_scan1(x = qtl6.2, map = map, lodcolumn = 1, main = colnames(qtl6.2)[1])

chr <- 5
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno =  df6[,8,drop = FALSE],
                     kinship = K[[chr]], addcovar = addcovar)

quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl6.2, main = colnames(qtl6.2)[1])

qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = df6[,8,drop = FALSE],
                      kinship = K[[chr]], addcovar = addcovar, cores = 4)

quartz()
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl6.2, main = colnames(qtl6.2)[1])

qtl6.3 <- scan1(genoprobs = probs, pheno = df6[,15,drop = FALSE],
                kinship = K, addcovar = addcovar, cores = 4, reml=TRUE)
quartz()
plot_scan1(x = qtl6.3, map = map, lodcolumn = 1, main = colnames(qtl6.3)[1])

chr <- 11
qtl.coef = scan1coef(genoprobs = probs[,chr], pheno =  df6[,15,drop = FALSE],
                     kinship = K[[chr]], addcovar = addcovar)

quartz()
plot(x = qtl.coef, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl6.3, main = colnames(qtl6.3)[1])

qtl.blup <- scan1blup(genoprobs = probs[,chr], pheno = df6[,15,drop = FALSE],
                      kinship = K[[chr]], addcovar = addcovar, cores = 4)

quartz()
plot(x = qtl.blup, map = map[[chr]], columns = 1:8, col = CCcolors,
     scan1_output = qtl6.3, main = colnames(qtl6.3)[1])

