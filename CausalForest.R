#memory.limit(10 * 10^10)
#Idea 2: Causal Forest
CausalForest <- function(data_df, tr = 2e3){
  
  data_df<-data_df[!is.na(data_df$yield),]
  
  result_grf <- list()
  result_grf[[1]] <- list()
  result_grf[[2]] <- list()
  
  W <- data_df$treatment
  Y <- data_df$yield
  
  C <- as.matrix(dplyr::select(data_df,
                               PC.W1,
                               PC.W2,
                               PC.W3,
                               PC.W4,
                               PC.L1,
                               PC.L2))
                               
  
  #single year
  X1 <- C
  X2 <- as.matrix(cbind(data_df$lon, data_df$lat, C)); colnames(X2)[1:2] <- c('lon','lat')
  #multiple year
  X3 <- as.matrix(cbind(data_df$year, C)); colnames(X3)[1] <- 'year'
  X4 <- as.matrix(cbind(data_df$year, X2)); colnames(X4)[1] <- 'year'
  
  for (j in 1:13){
    
    ind <- data_df$year == (j + 2007)
    result_grf[[1]][[j]] <- causal_forest(X1[ind,], Y[ind], W[ind], num.trees = tr)
    result_grf[[2]][[j]] <- causal_forest(X2[ind,], Y[ind], W[ind], num.trees = tr)
    
  }
  
  result_grf[[3]] <- causal_forest(X3, Y, W, num.trees = tr)
  result_grf[[4]] <- causal_forest(X4, Y, W, num.trees = tr)
  
  return(result_grf)
}


ate.CFF.abject <- function(res){
  
  ate_non_sptial <- rep(0, 14)
  ate_non_sptial.se <- rep(0, 14)
  ate_sptial <- rep(0, 14)
  ate_sptial.se <- rep(0, 14)
  
  ate_non_sptial <- rep(0, 13)
  ate_non_sptial.se <- rep(0, 13)
  ate_sptial <- rep(0, 13)
  ate_sptial.se <- rep(0, 13)
  
  for (j in 1:13){
    
    ate_cf = average_treatment_effect(res[[1]][[j]], target.sample = 'control')
    ate_cf_sp = average_treatment_effect(res[[2]][[j]], target.sample = 'control')
    
    ate_non_sptial[j] <- ate_cf[1]
    ate_non_sptial.se[j] <- ate_cf[2]
    ate_sptial[j] <- ate_cf_sp[1]
    ate_sptial.se[j] <- ate_cf_sp[2]
    
  }
  
  ate_cf = average_treatment_effect(res[[3]], target.sample = 'control')
  ate_cf_sp = average_treatment_effect(res[[4]], target.sample = 'control')

  ate_non_sptial[14] <- ate_cf[1]
  ate_non_sptial.se[14] <- ate_cf[2]
  ate_sptial[14] <- ate_cf_sp[1]
  ate_sptial.se[14] <- ate_cf_sp[2]
  
  ate_repot <- rbind(ate_non_sptial, 
                     ate_non_sptial.se,
                     ate_sptial,
                     ate_sptial.se)  
  return(ate_repot)
}

ate.CFF.abject.matching <- function(res, data_df_matching){
  
  ate_non_sptial <- rep(0, 14)
  ate_non_sptial.se <- rep(0, 14)
  ate_sptial <- rep(0, 14)
  ate_sptial.se <- rep(0, 14)
  
  ate_non_sptial <- rep(0, 13)
  ate_non_sptial.se <- rep(0, 13)
  ate_sptial <- rep(0, 13)
  ate_sptial.se <- rep(0, 13)
  
  for (j in 1:13){
    
    data_df_matching_fyear <- data_df %>% filter(year == (j + 2007))
    
    data_df_matching_fyear$ID <- 1:dim(data_df_matching_fyear)[1]
    
    data_df_matching_fyear <- match.obj(data_df_matching_fyear, 'spatial')
    
    ate_cf = average_treatment_effect(res[[1]][[j]], subset = data_df_matching_fyear$ID)
    ate_cf_sp = average_treatment_effect(res[[2]][[j]], subset = data_df_matching_fyear$ID)
    
    ate_non_sptial[j] <- ate_cf[1]
    ate_non_sptial.se[j] <- ate_cf[2]
    ate_sptial[j] <- ate_cf_sp[1]
    ate_sptial.se[j] <- ate_cf_sp[2]
    
  }
  
  ate_cf = average_treatment_effect(res[[3]], subset = data_df_matching$ID)
  ate_cf_sp = average_treatment_effect(res[[4]], subset = data_df_matching$ID)

  ate_non_sptial[14] <- ate_cf[1]
  ate_non_sptial.se[14] <- ate_cf[2]
  ate_sptial[14] <- ate_cf_sp[1]
  ate_sptial.se[14] <- ate_cf_sp[2]
  
  ate_repot <- rbind(ate_non_sptial, 
                     ate_non_sptial.se,
                     ate_sptial,
                     ate_sptial.se)  
  return(ate_repot)
}



