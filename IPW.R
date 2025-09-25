IPW.obj <- function(data_df, model = 'non-spatial'){
  
  print(paste('before:', dim(data_df)[1]))
  
  if (model == 'non-spatial'){
    IPW.fit <- bam(treatment ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean), 
                   family = binomial,
                   data = data_df)
  }else if (model == 'spatial'){
    IPW.fit <- bam(treatment ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean)+
                     s(lon, lat, k = k), 
                   family = binomial,
                   data = data_df)
  }else if (model == 'temporal'){
    IPW.fit <- bam(treatment ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean), 
                   family = binomial,
                   data = data_df)
  }else{
    IPW.fit <- bam(treatment ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean)+
                     s(lon, lat, k = k), 
                   family = binomial,
                   data = data_df)
    
  }
  
  IPW.Trimming.index <- !((data_df$prop.scores.IPW > 0.9) * !(data_df$prop.scores.IPW < 0.1))
  print(paste('after:', sum(IPW.Trimming.index)))
  data_df <- data_df[IPW.Trimming.index,]
  
  n <- nrow(data_df)
  
  z <- data_df$treatment
  
  p <- IPW.fit$fitted.values
  
  denominator <- sum(z/p + (1 - z)/(1 - p))
  
  data_df$weights <- (z/p + (1 - z)/(1 - p)) / denominator *n
  
  return(data_df)
}

###non spatial-