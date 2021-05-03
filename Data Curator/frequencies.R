#Importing necessary packages.
library(dplyr)
library(ggplot2)

#Reading in Individual counts for barcodes.
df_604 <- read.table("SRR3879604_1_bc_COUNTS")
df_605 <- read.table("SRR3879605_1_bc_COUNTS")
df_606 <- read.table("SRR3879606_1_bc_COUNTS")

#Reading in Combined counts for barcodes.
df_combined <- read.table("Counts/combined_counts")

#Generating summary statistics for Combined barcodes.
summary(df_combined$V1)

#Filtering counts < mean counts.
df_filtered <- filter(df_combined,df_combined$V1>mean(df_combined$V1))
write.table(df_filtered$V2,file="whitelist_from_mean_counts.txt",row.names = FALSE,col.names = FALSE,quote=FALSE)

#Plotting CDF for individual samples.
par(mfrow=c(1,3))
plot(ecdf(df_604$V1),xlab="Barcodes",ylab="Density",main="SRR3879604_1")
plot(ecdf(df_605$V1),xlab="Barcodes",ylab="Density",main="SRR3879605_1")
plot(ecdf(df_606$V1),xlab="Barcodes",ylab="Density",main="SRR3879606_1")