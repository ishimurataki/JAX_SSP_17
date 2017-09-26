padams_data <- read.csv(file = "/Users/s-ishimt/Desktop/Workbook3.csv", as.is = TRUE)
View(padams_data)
colnames(padams_data) <- c("Sample", "Value", "SE", "Genotype")

padams_data <- padams_data[which(padams_data$Value != "x"), ]
padams_data$SE <- as.numeric(padams_data$SE)
padams_data$Value <- as.numeric(padams_data$Value)

padams_plot <- ggplot(data = padams_data, mapping = aes(x = Value, y = Genotype)) +
  geom_point() + geom_line(x = 
                             c(mean(padams_data$Value[padams_data$Genotype == "B"]), 
                               mean(padams_data$Value[padams_data$Genotype == "D"])), 
                           y = padams_data$Genotype)

quartz()
plot(padams_plot)

quartz()
par(mfrow = c(2,1))
hist(padams_data$Value[padams_data$Genotype == "D"])
hist(padams_data$Value[padams_data$Genotype == "B"])

quartz()
par(mfrow = c(2,1))
hist(log2(padams_data$Value[padams_data$Genotype == "D"]))
hist(log2(padams_data$Value[padams_data$Genotype == "B"]))
