GAM_function <- function(data_df, method = 'GAM', k = k){
  
  #####PCA component
  spatial <- ' s(lon, lat, k = k) +'
  temporal <- ' fyear +'
  spatial_tem <- ' s(lon, lat, k = k) + fyear +'

  treatment1 <- paste(c('treatment + s(PC.W.1xTRT) + s(PC.L.1xTRT)+'))
  treatment2 <- paste(c('s(lon, lat, k = k, by = treatment)',
                        's(lon, lat, k = k, by = PC.W.1xTRT)',
                        's(lon, lat, k = k, by = PC.L.1xTRT) +'),
                      collapse = " + ")
  
  covariate_n <- paste(c('s(PC.W1)',
                         's(PC.W2)',
                         's(PC.W3)',
                         's(PC.W4)',
                         's(PC.L1)',
                         's(PC.L2)',
                         's(PC.W.1xPC.L.1)'),
                       collapse = " + ")

  covariate_s <- paste(c('s(lon, lat, k = k, by = PC.W1)',
                         's(lon, lat, k = k, by = PC.W2)',
                         's(lon, lat, k = k, by = PC.W3)',
                         's(lon, lat, k = k, by = PC.W4)',
                         's(lon, lat, k = k, by = PC.L1)',
                         's(lon, lat, k = k, by = PC.L2)',
                         's(lon, lat, k = k, by = PC.W.1xPC.L.1)'),
                         collapse = " + ")
  
  
  #single year
  #non-spatial
  formula1 <- as.formula(paste0("yield ~ ", treatment1, covariate_n))
  #Spatial
  formula2 <- as.formula(paste0("yield ~ ", treatment2, spatial, covariate_s))
  #Multiple year
  #non-spatial
  formula3 <- as.formula(paste0("yield ~ ", treatment1, temporal, covariate_n))
  #Spatial
  formula4 <- as.formula(paste0("yield ~ ", treatment2, spatial_tem, covariate_s))
  
  
  fit <- list()
  fit[[1]] <- list()
  fit[[2]] <- list()
    
  for (j in 1:13){
    
    data_df_fyear <- data_df %>% filter(year == (j + 2007))
      
    if (method == 'IPW'){
        data_df_tmp <- IPW.obj(data_df_fyear, 'non-spatial')
        fit[[1]][[j]] <- bam(formula = formula1, data = data_df_tmp, weights = weights)
        data_df_tmp <- IPW.obj(data_df_fyear, 'spatial')
        fit[[2]][[j]] <- bam(formula = formula2, data = data_df_tmp, weights = weights)
    }else if(method == 'match'){
        data_df_tmp <- match.obj(data_df_fyear, 'non-spatial')
        fit[[1]][[j]] <- bam(formula = formula1, data = data_df_tmp, weights = weights)
        data_df_tmp <- match.obj(data_df_fyear, 'spatial')
        fit[[2]][[j]] <- bam(formula = formula2, data = data_df_tmp, weights = weights)
    }else{
      print('gam b')
      fit[[1]][[j]] <- bam(formula = formula1, data = data_df_fyear)
      print('gam m')
      fit[[2]][[j]] <- bam(formula = formula2, data = data_df_fyear)
    }
    
    print(j)
  }
  
  if (method == 'IPW'){
    data_df_tmp <- IPW.obj(data_df, 'temporal')
    fit[[3]] <- bam(formula = formula1, data = data_df_tmp, weights = weights)
    data_df_tmp <- IPW.obj(data_df, 'spatio-temporal')
    fit[[4]] <- bam(formula = formula2, data = data_df_tmp, weights = weights)
  }else if(method == 'match'){
    data_df_tmp <- match.obj(data_df, 'temporal')
    fit[[3]] <- bam(formula = formula1, data = data_df_tmp, weights = weights)
    data_df_tmp <- match.obj(data_df, 'spatio-temporal')
    fit[[4]] <- bam(formula = formula2, data = data_df_tmp, weights = weights)
  }else{
    fit[[3]] <- bam(formula = formula3, data = data_df)
    fit[[4]] <- bam(formula = formula4, data = data_df)
  }
  
  return(fit)
}

ate.GAM.object <- function(res, method){
  
  ate_non_sptial <- rep(0, 14)
  ate_non_sptial.se <- rep(0, 14)
  ate_sptial <- rep(0, 14)
  ate_sptial.se <- rep(0, 14)
  
  ate_non_sptial <- rep(0, 13)
  ate_non_sptial.se <- rep(0, 13)
  ate_sptial <- rep(0, 13)
  ate_sptial.se <- rep(0, 13)
  
  for (j in 1:13){
    
    data_df_fyear <- data_df %>% filter(year == (j + 2007))
    
    if (method == 'IPW'){
      data_df_fyear <- IPW.obj(data_df_fyear)
    }else if(method == 'match'){
      data_df_fyear <- match.obj(data_df_fyear)
    }else{data_df_fyear$weights = 1}
    
    ate_non_tmp <- ate_cal(res[[1]][[j]], data_df_fyear, 'non-spatial')
    ate_spatial_tmp <- ate_cal(res[[2]][[j]], data_df_fyear, 'spatial')
    
    ate_non_sptial[j] <- ate_non_tmp$fit
    ate_non_sptial.se[j] <- ate_non_tmp$se.fit
    ate_sptial[j] <- ate_spatial_tmp$fit
    ate_sptial.se[j] <- ate_spatial_tmp$se.fit
  }
  
  ate_non_tmp <- ate_cal(res[[3]], data_df, 'non-spatial')
  ate_spatial_tmp <- ate_cal(res[[4]], data_df, 'spatial')

  ate_non_sptial[14] <- ate_non_tmp$fit
  ate_non_sptial.se[14] <- ate_non_tmp$se.fit
  ate_sptial[14] <- ate_spatial_tmp$fit
  ate_sptial.se[14] <- ate_spatial_tmp$se.fit

  ate_repot <- rbind(ate_non_sptial, 
                     ate_non_sptial.se,
                     ate_sptial,
                     ate_sptial.se)
  
  return(ate_repot)
}


ate.GAM.object.matching <- function(res, data_df_matching, method){
  
  ate_non_sptial <- rep(0, 14)
  ate_non_sptial.se <- rep(0, 14)
  ate_sptial <- rep(0, 14)
  ate_sptial.se <- rep(0, 14)
  
  ate_non_sptial <- rep(0, 13)
  ate_non_sptial.se <- rep(0, 13)
  ate_sptial <- rep(0, 13)
  ate_sptial.se <- rep(0, 13)
  
  for (j in 1:13){
    
    year <- j + 2007
    data_df_fyear <- data_df %>% filter(year == year)
    
    data_df_fyear <- match.obj(data_df_fyear, 'spatial')
    
    ate_non_tmp <- ate_cal(res[[1]][[j]], data_df_fyear, 'non-spatial')
    ate_spatial_tmp <- ate_cal(res[[2]][[j]], data_df_fyear, 'spatial')
    
    ate_non_sptial[j] <- ate_non_tmp$fit
    ate_non_sptial.se[j] <- ate_non_tmp$se.fit
    ate_sptial[j] <- ate_spatial_tmp$fit
    ate_sptial.se[j] <- ate_spatial_tmp$se.fit
  }
  
  ate_non_tmp <- ate_cal(res[[3]], data_df_matching, 'non-spatial')
  ate_spatial_tmp <- ate_cal(res[[4]], data_df_matching, 'spatial')

  ate_non_sptial[14] <- ate_non_tmp$fit
  ate_non_sptial.se[14] <- ate_non_tmp$se.fit
  ate_sptial[14] <- ate_spatial_tmp$fit
  ate_sptial.se[14] <- ate_spatial_tmp$se.fit
  
  ate_repot <- rbind(ate_non_sptial, 
                     ate_non_sptial.se,
                     ate_sptial,
                     ate_sptial.se)
  
  return(ate_repot)
}

##other function for GAM object
ate_cal <- function(res_tmp, data_df_new, method = 'non-spatial'){
  
  ate_report <- list()
  
  data_df_new$treatment <- 1
  
  if (method == 'non-spatial'){
    predict.tmp <- predict(res_tmp, 
                           newdata = data_df_new, 
                           type = "terms", 
                           terms = "treatment",
                           se.fit = TRUE)
    ate_report$fit <- mean(predict.tmp$fit)
    ate_report$se.fit <- predict.tmp$se.fit[1]
  }else{
    predict.tmp <- predict(res_tmp, 
                                  newdata = data_df_new, 
                                  type = "terms", 
                                  terms = "s(lon,lat):treatment",
                                  se.fit = TRUE)
    ate_report$fit <- mean(predict.tmp$fit)
    ate_report$se.fit <- ate_se_spatial(res_tmp, data_df_new)
    
  }
  
  return(ate_report)
}


ate_se_spatial <- function(res_tmp, data_df_new){
  
  predict.tmp <- predict(res_tmp, 
                         newdata = data_df_new, 
                         type="lpmatrix", 
                         terms = "s(lon,lat):treatment")
  
  sp.index <- predict.tmp[1,] != 0
  n <- dim(data_df_new)[1]
  X.lp <- predict.tmp[, sp.index]
  X <- apply(X.lp, MARGIN = 2, FUN = mean)
  Vp <- res_tmp$Vp[sp.index, sp.index]
  tr.Vp <- t(X) %*% Vp %*% X
  
  return(sqrt(as.numeric(tr.Vp)))
}












