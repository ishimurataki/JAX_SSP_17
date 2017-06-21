########################################################
# DO Data Importation and Preliminary Investigations #
########################################################

# remove all preexisiting objects and variables in the environmnet
rm(list=ls())

# set working directory to the physical location of the data
setwd("/Users/s-ishimt/Desktop/DO_data")

#importation of the DO data
load("DO378_islet.RData")
pheno_data <- read.csv("pheno_clin.csv", as.is = TRUE)

# creating subset of pheno_data that matches the RNAseq data
pheno_data$mouse[nchar(pheno_data$mouse)==5] <- sub(pattern = "DO-", replacement = "DO-0", x = pheno_data$mouse[nchar(pheno_data$mouse)==5])
colnames(pheno_data)[1] <- "Mouse.ID"

pheno_data <- pheno_data[match(x=annot.samples$Mouse.ID, table=pheno_data$Mouse.ID),]
stopifnot(annot.samples$Mouse.ID == pheno_data$Mouse.ID)

# exportation of newly created csv file
write.csv(x=pheno_data, file="usable_clin_pheno_data", row.names = FALSE)

# create data frame with just food consumption
food_columns <- which(substr(colnames(pheno_data), 1, 5) == "food_")[
  -length(which(substr(colnames(pheno_data), 1, 5) == "food_"))]

food_data <- pheno_data[,food_columns]
food_data2 <- stack(food_data)
week_numbers <- as.numeric(sub(pattern = "wk", replacement = "", x = 
  sub(pattern = "food_", replacement = "", food_data2$ind)
))
food_data2 <- cbind(food_data2, week_numbers)

# creating cool graphics with food consumption
ggplot(data=food_data2, mapping = aes(y = values, x = week_numbers)) + geom_point() 

# create data frame with just weight consumption
weight_columns <- which(substr(colnames(pheno_data), 1, 7) == "weight_")[
  -length(which(substr(colnames(pheno_data), 1, 7) == "weight_"))]

weight_data <- pheno_data[,weight_columns]
weight_data2 <- stack(weight_data)
week_weight_numbers <- as.numeric(sub(pattern = "wk", replacement = "", x = 
                                 sub(pattern = "weight_", replacement = "", weight_data2$ind)
))
weight_data2 <- cbind(weight_data2, week_weight_numbers)

Mouse.ID <- pheno_data$Mouse.ID
weight_data2 <- cbind(Mouse.ID, weight_data2)

weight_data2 <- cbind(weight_data2, pheno_data$coat_color)
colnames(weight_data2)[5] <- "coat_color"

weight_data2 <- cbind(weight_data2, pheno_data$Glu_tAUC)
colnames(weight_data2)[7] <- "Glu_tAUC"



# creating cool graphics with weight
ggplot(data=weight_data2, mapping = aes(y = values, x = week_weight_numbers, by = Mouse.ID, color = Glu_tAUC))+ geom_line(scale_fill_distiller(pallete="Blues"))



for(i in 1:157){
  if((class(pheno_data[,i])=="numeric") == TRUE){
    print(cor(pheno_data[,i], pheno_data[,104], use = "complete.obs"))
    print(colnames(pheno_data)[i])
  }
}





