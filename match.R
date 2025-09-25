# library(devtools)
# devtools::install_github("gpapadog/DAPSm")
# Sys.setenv('R_MAX_VSIZE'=32000000000)
# data_df$prop.scores <- glm(treatment ~ NCCPI + AVG_AWC + AVG_SAND + AVG_CLAY, family = binomial,
#                             data = data_df)$fitted.values
# 
# daps2 <- DAPSest(data_df, 
#                  out.col = 8, 
#                  trt.col = 131, 
#                  caliper = 0.3,
#                  weight = 'optimal', 
#                  coords.columns = c(4, 5),
#                  pairsRet = TRUE, 
#                  cov.cols = c(81, 126, 127, 129), 
#                  cutoff = 0.15,
#                  w_tol = 0.001, 
#                  coord_dist = TRUE, 
#                  caliper_type = 'DAPS',
#                  matching_algorithm = 'greedy')

match.obj <- function(data_df, model = 'non-spatial', seeds = 94703){
  
  set.seed(seeds)
  
  data_df$treatment.matching <- as.integer(!as.logical(data_df$treatment))
  
  if (model == 'non-spatial'){
    m.out <- matchit(treatment.matching ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean), 
                       data = data_df, 
                       method = "nearest", 
                       distance = "gam")
  }else if (model == 'spatial'){
    m.out <- matchit(treatment.matching ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean) +
                       s(lon, lat, k = k), 
                       data = data_df, 
                       method = "nearest", 
                       distance = "gam")
  }else if (model == 'temporal'){
    m.out <- matchit(treatment.matching ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean), 
                     data = data_df, 
                     method = "nearest", 
                     distance = "gam")
  }else{
    m.out <- matchit(treatment.matching ~ s(nccpi3corn_mean) + s(rootznaws_mean) + s(soc0_100_mean) +
                       s(lon, lat, k = k), 
                     data = data_df, 
                     method = "nearest", 
                     distance = "gam")
    
  }

  # plot(m.out.1, type = "jitter", interactive = FALSE)
  # plot(m.out.1, type = "density", interactive = FALSE, 
  #      which.xs = ~ NCCPI + AVG_AWC + AVG_SAND + AVG_CLAY)
  # summary(m.out.1, un = FALSE)
  
  m.data <- match.data(m.out)
  return(m.data)
}

# functional diversity (couple way) yield touch weather
# functional diversity: where we categroy to rotation to small number of group according 
# their ecolocal function. estalish way to doing it --- C3 C4 graph broadly plan
# 3 or 4 5 come up with. What's the in the current year? what the history function
# diversity affect yield in the current year, especially yield from the current
# what's relative of ---- building up specific --- 
# 1. hegieotory coefficient effect
# 2. the stressful weather, any interaction the weather --- truning weather in
# PCA, interpret the PCA, 
# 3. data-driven --- generic ml method like  --- yield - covariate + previous crop
# hot-one-coding, 






# 1 issue matching include distance
# 2 sptial propensiity score for matching, if no reason, eason when -----


# ggplot(m.data.1,
#        aes(lon, lat, col = fannual_RCI)) + geom_point(size = .5) + coord_equal()
# 
# # head(m.data.1)
# 
# fit.1 <- lm(yield ~ fannual_RCI * NCCPI, data = m.data.1, weights = weights)
# 
# att.1 = avg_comparisons(fit.1,
#                         variables = "fannual_RCI",
#                         vcov = ~subclass,
#                         newdata = subset(m.data.1, fannual_RCI == 1),
#                         wts = "weights")

# data('toyData2')
# 
# toyData2$prop.scores <- glm(Z ~ X1 + X2 + X3 + X4, family = binomial,
#                             data = toyData2)$fitted.values
# 
# daps <- DAPSest(toyData2, out.col = 2, trt.col = 1, caliper = 0.3,
#                 weight = 0.7, coords.columns = c(4, 5),
#                 pairsRet = TRUE, cov.cols = 6:9, cutoff = 0.1,
#                 coord_dist = TRUE, caliper_type = 'DAPS',
#                 matching_algorithm = 'greedy')
# 
# bal <- CalcDAPSWeightBalance(toyData2, weights = seq(0, 1, length.out = 40),
#                              cov.cols = 6:9, trt.col = 1,
#                              coords.columns = c(4, 5), caliper = 0.3,
#                              matching_algorithm = 'greedy')
# PlotWeightBalance(bal$balance, weights = seq(0, 1, length.out = 40), cutoff = 0.15)
# DAPS <- DAPSchoiceModel(toyData2, trt.col = 1, balance = bal$balance,
#                         cutoff = 0.15, pairs = bal$pairs,
#                         weights = seq(0, 1, length.out = 40))
# f
# MatchedDataMap(x = bal$full_pairs[[10]], trt_coords = c(3, 4),
#                con_coords = c(7, 8))




