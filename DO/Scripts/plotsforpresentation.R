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
Ins_per_islet <- log(pheno_data$Ins_per_islet)

probs <- probs_doqtl_to_qtl2(probs = genoprobs, map = snps, chr_column = "chr", pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, chr_column = "chr", pos_column = "bp")
annot.samples$diet_days <- pheno_data$diet_days
K = calc_kinship(probs = probs, type = "loco", cores = 4)

## Insulin secretion dependent on glucose stimulation plot
G33_ins <- cbind(pheno_data[,1:2],G33_ins_sec_norm, rep("3.3 mM", 378))
colnames(G33_ins) <- c("Mouse_ID", "Sex","Secretion", "Stimulation_levels")
G83_ins <- cbind(pheno_data[,1:2],G83_ins_sec_norm, rep("8.3 mM", 378))
colnames(G83_ins) <- c("Mouse_ID", "Sex","Secretion", "Stimulation_levels")
G167_ins <- cbind(pheno_data[,1:2],G167_ins_sec_norm, rep("16.7 mM", 378))
colnames(G167_ins) <- c("Mouse_ID", "Sex","Secretion", "Stimulation_levels")

Insulin_sec_df <- rbind(G33_ins, G83_ins, G167_ins)

ggplot(data = Insulin_sec_df, mapping = aes(x = Stimulation_levels, y = Secretion, color = Sex, group = Mouse_ID)) + geom_jitter(width = 0.1) + 
  geom_line(data = Insulin_sec_df[Insulin_sec_df$Mouse_ID %in% c("DO-279", "DO-254"),], color = "#8488F9", size = 1) +
  labs(title = "Insulin Secretion vs Glucose Stimulation Levels", x = "Glucose Stimulation Levels", y = "Insulin Secretion") +
  scale_x_discrete(labels = c("3.3 mM", "8.3 mM", "16.7 mM"))

# Insulin Secretion is mediated stoichiometrically by current insulin levels in the islet
Ins_per_islet_df <- cbind(pheno_data[,1:2], G83_ins_sec_norm, Ins_per_islet)
quartz()
ggplot(data = Ins_per_islet_df, mapping = aes(x = Ins_per_islet, y = G83_ins_sec_norm, color = sex)) + geom_point(size = 1.4) + 
  labs(x = "Islet Insulin Levels", y = "Insulin Secretion (8.3 mM Glucose Stim.)") +
  ggtitle("Insulin Secretion vs Islet Insulin Levels", subtitle = paste("Pearson's Correlation:", 
                                round(cor(Ins_per_islet, G83_ins_sec_norm, use = "complete.obs"), digits = 3), sep = " "))



