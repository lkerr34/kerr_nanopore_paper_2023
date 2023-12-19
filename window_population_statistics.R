# In this script we analyse the methylation levels for CpGs in genomic 
# windows of the genome

# Libraries/directories ----
# Load (and install if necessary) the packages required to read the data

# install.packages("readr")
gc()
library(readr)

args <- commandArgs(trailingOnly=TRUE);
input <- args[1];
output <- args[2];

# Load the CpG methylation data

data <- read.table(paste(input), sep="\t", header=FALSE)
colnames(data)<-c("chromosome","start","end","coverage","methylated","meth_levels","window_name")

# Extract the methylation levels corresponding to each individual window
split_window <- split(data,data$window_name);

# Calculate statistics associated with the methylation level
# In particular, we calculate the number of CpGs as well as the mean, variance, 
# coefficient of variation and correlation

window_stats <- sapply(split_window, function(df)
{
  window <- df$window_name[[1]];
  meth <- df$meth_levels;
  num <- length(meth);
  mean_dis <- (df$end[length(df$end)]-df$start[1])/(num-1);
  mu <- mean(meth);
  var <- var(meth);
  CV <- var^0.5/mu;
  if(num >= 3 && var(meth[1:length(meth)-1])>0 && var(meth[2:length(meth)])>0)
  {cors <- cor(meth[1:length(meth)-1],meth[2:length(meth)])}
  else {cors <- NA};
  c(num, mean_dis, mu, var, CV, cors, window);
}
)
rownames(window_stats)<-c("CpGs","mean_distance","mean","variance","CV","correlation","window_name");
window_stats <- as.data.frame(t(window_stats));

# If we want to visualise data using other software (e.g. Mathematica) then we 
# can save statistics as tsv files

write_tsv(window_stats, file=paste(output,"gz",sep="."))


