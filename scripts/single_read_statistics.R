# Conduct analysis of single-molecule DNA methylation patterns

# Packages and arguments ####

library(readr)

args <- commandArgs(trailingOnly=TRUE);
input <- args[1];
output <- args[2];
threshold <- as.numeric(args[3]);


# Load in required data ####

data <- read.table(paste(input), sep="\t", header=FALSE);
colnames(data)<-c("chromosome","start","end","name","meth_sequence","log_like_ratio","context")



# Calculate statistics for individual reads ####

stats <- apply(data, 1, function(df)
{
  loglikeratio <- as.numeric(unlist(strsplit(df[[6]],",")));
  
  binary <- rep(NA,length(loglikeratio));
  
  methind <- which(loglikeratio>=threshold);
  unmethind <- which(-(loglikeratio)>=threshold);
  binary[methind] <- 1;
  binary[unmethind] <- 0;
  binary <- na.omit(binary);
  num <- length(binary);
  mu <- mean(binary);
  var <- var(binary);
  CV <- sqrt(var)/mu;
  
  if(mu!=0 & mu!=1 & is.na(mu)==FALSE){
    SE <- (1-mu)*log2(1/(1-mu))+mu*log2(1/mu)  
  }
  else{if(is.na(mu)==FALSE){
    SE <- 0
  }
    else{
      SE <- NA_real_
    }
  };
  
  if(num>=2){
    changes <- sum(abs(binary[1:(num-1)]-binary[2:num]))/(num-1);
    runs <- rle(as.integer(binary));
    average_run_length <- sum(runs[[1]])/length(runs[[1]]);
  }
  else{
    changes <- NA;
    average_run_length <- NA
  }
  
  if(num>=3 && var(binary[1:(num-1)])>0 && var(binary[2:num])>0){
    corr <- cor(binary[1:(num-1)],binary[2:num]);
  }
  else{
    corr <- NA
  }
  
  c(num, mu, var, CV, corr, changes, SE, average_run_length)
  
}
)

data$num_CpGs <- stats[1,];
data$mean_meth_level <- stats[2,];
data$var_meth_level <- stats[3,];
data$cv_meth_level <- stats[4,];
data$cor_meth_level <- stats[5,];
data$changes <- stats[6,];
data$SE <- stats[7,];
data$average_run_length <- stats[8,];

write_tsv(data, file=gzfile(paste(output,"gz",sep=".")))

