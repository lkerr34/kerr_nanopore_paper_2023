# Calculate single molecule distance-dependent correlations between neighbouring CpGs

args <- commandArgs(trailingOnly=TRUE);
input <- args[1];
output <- args[2];
threshold <- as.numeric(args[3]);
maxdist <- as.numeric(args[4]);

# Packages ####

library(readr)

# Load in required data ####

df <- read.table(paste(input), sep="\t", header=FALSE);
colnames(df)<-c("chromosome","start","end","read_name","log_like_ratios_1","log_like_ratios_2","distances")

# Extract the required data

required_data <- df[,c(4,5,6,7)];


results <- apply(required_data, 1, function(data)
{
  read_name <- data[[1]];
  log_lik_1 <- as.numeric(unlist(strsplit(data[[2]],",")));
  log_lik_2 <- as.numeric(unlist(strsplit(data[[3]],",")));
  distances <- as.numeric(unlist(strsplit(data[[4]],",")));
  
  print(read_name);
  
  # Create new variables (1=methylated; 0=unmethylated; NA=unassigned)
  
  methind1 <- which(log_lik_1>=threshold);
  unmethind1 <- which(-(log_lik_1)>=threshold);

  binary1 <- rep(NA,length(log_lik_1));
  binary1[methind1] <- 1;
  binary1[unmethind1] <- 0;
  
  methind2 <- which(log_lik_2>=threshold);
  unmethind2 <- which(-(log_lik_2)>=threshold);
  
  binary2 <- rep(NA,length(log_lik_2));
  binary2[methind2] <- 1;
  binary2[unmethind2] <- 0;
  
  # Combine new variables into a dataframe with distance info and remove missing data
  
  new_mat <- matrix(nrow=,ncol=3);
  new_df <- data.frame(new_mat);
  
  for (i in 1:length(binary1))
  {
    for (j in i:length(binary2))
    {
      d <- sum(distances[i:j]);
      if(d>maxdist)
      {break}
        else
          {
      new_df[nrow(new_df) + 1,] = c(binary1[i],binary2[j],d)
      }
    }
  }
  new_df_excl <- data.frame(na.omit(new_df));
  rm(new_df);
  colnames(new_df_excl) <- c("state_1","state_2","distances");
  

  # Split the new dataframe based on the distance between adjacent CpGs and calculate correlations where possible
  
  split_df <- split(new_df_excl,new_df_excl$distances);

    correlations <- sapply(split_df, function(df)
    {
      if (dim(df)[1]>1 && var(df[,1]!=0) && var(df[,2])!=0)
      {
        corr <- cor(df[,1],df[,2]);
      }
      else
      {
        corr <- NA;
      }
      corr
    }
    );

  cbind(read_name,correlations);

}
);


for(i in 1:length(results)){
  read_results <- data.frame(results[[i]]);
  read_name <- read_results[1,1];
  read_results_excl <- read_results[which(read_results$correlations!="NA"),];
  if(length(read_results)>1){
    if(length(read_results_excl[[1]])>0){
    results_to_save <- data.frame(cbind(rownames(read_results_excl),read_results_excl[,2]));
    colnames(results_to_save)=c("distance","correlation");
    write_tsv(results_to_save,file=paste0(output,read_name,".tsv"));
  }
}
}

