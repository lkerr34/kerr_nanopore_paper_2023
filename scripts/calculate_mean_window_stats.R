# Calculate mean statistics for each window

# Packages ####

library(readr)

args <- commandArgs(trailingOnly=TRUE);
input <- args[1];
output <- args[2];

# Load in required data ####

data <- read.table(paste(input), sep="\t", header=FALSE);
colnames(data)<-c("chromosome","start","end","name","num_CpGs","mean_meth","variance","CV","correlation","RTS","SE","mean_run_length","window_name")

# Calculate average stats for each window ####

# Remove all reads with less than 100 CpGs

data_excl <- data[which(data$num_CpGs>=100),];

# Split data based on windows

split_data <- split(data_excl,data_excl$window_name);

window_mean_stats <- sapply(split_data, function(df)
{
  window_name <- df$window_name[[1]];
  mean_meth <- mean(df$mean_meth, na.rm=TRUE);
  mean_var <- mean(df$variance, na.rm=TRUE);
  mean_CV <- mean(df$CV, na.rm=TRUE);
  mean_cor <- mean(df$correlation, na.rm=TRUE);
  mean_RTS <- mean(df$RTS, na.rm=TRUE);
  mean_SE <- mean(df$SE, na.rm=TRUE);
  c(window_name, mean_meth, mean_var, mean_CV, mean_cor, mean_RTS, mean_SE);
}
)
rownames(window_mean_stats)<-c("window_name","mean_meth","mean_variance","mean_CV","mean_correlation","mean_RTS","mean_SE");
window_mean_stats <- as.data.frame(t(window_mean_stats))

# Save output

write_tsv(window_mean_stats, file=paste(output,"gz",sep="."))
