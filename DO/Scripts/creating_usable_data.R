##############################################################################
# Creating Usable Data from Updated DO data provided by Carl Broman 06/27/17 #
##############################################################################

# creating subset of pheno_data that matches the RNAseq data
pheno_clin$mouse[nchar(pheno_clin$mouse)==5] <- sub(pattern = "DO-", replacement = "DO-0", x = pheno_clin$mouse[nchar(pheno_clin$mouse)==5])
colnames(pheno_clin)[1] <- "Mouse.ID"

pheno_clin <- pheno_clin[match(x=annot.samples$Mouse.ID, table=pheno_clin$Mouse.ID),]
stopifnot(annot.samples$Mouse.ID == pheno_clin$Mouse.ID)

# exportation of newly created csv file
write.csv(x=pheno_clin, file="/Users/s-ishimt/Desktop/Jax_SSP_17/DO/DO_data/usable_clin_pheno_data.csv", row.names = FALSE)