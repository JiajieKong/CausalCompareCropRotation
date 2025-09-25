add.covariates <- function(data_df, seed.num = 94703){
  
  set.seed(seed.num)
  
  C.1 <- as.matrix(dplyr::select(data_df,
                                 soil_6,
                                 soil_7,
                                 soil_8,
                                 def_6,
                                 def_7,
                                 def_8,
                                 pr_6,
                                 pr_7,
                                 pr_8,
                                 tmmn_6,
                                 tmmn_7,
                                 tmmn_8,
                                 tmmx_6,
                                 tmmx_7,
                                 tmmx_8,
                                 vpdmax_6,
                                 vpdmax_7,
                                 vpdmax_8,
                                 cGDD_6m,
                                 cGDD_7m,
                                 cGDD_8m))
  
  C.2 <- as.matrix(dplyr::select(data_df,
                                 rootznaws_mean,
                                 soc0_100_mean,
                                 nccpi3corn_mean))
  
  
  pc.1 <- prcomp(C.1, center = TRUE, scale. = TRUE) 
  pc.2 <- prcomp(C.2, center = TRUE, scale. = TRUE)
  
  print(pc.1[[1]]^2/sum(pc.1[[1]]^2))
  print(pc.2[[1]]^2/sum(pc.2[[1]]^2))
  
  colnames(pc.1$x) <- paste('PC.W', 1:dim(pc.1$x)[2], sep = '')
  colnames(pc.2$x) <- paste('PC.L', 1:dim(pc.2$x)[2], sep = '')
  
  pc.i <- cbind(pc.1$x[,1] * pc.2$x[,1],
                pc.1$x[,1] * data_df$treatment,
                pc.2$x[,1] * data_df$treatment)
  
  colnames(pc.i) <- c('PC.W.1xPC.L.1', 
                      'PC.W.1xTRT', 
                      'PC.L.1xTRT')
  
  data_df <- cbind(data_df, pc.1$x, pc.2$x, pc.i)
  
  return(data_df)
}
