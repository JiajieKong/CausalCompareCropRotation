#datahub berkeley
rm(list=ls())
library(cowplot)
library(tidyverse)
library(ggplot2)
library(mgcv)
library(grf)
library(MatchIt)
library(marginaleffects)
library(arrow)
library(metafor)

current_path = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_path)

source('GAM.R')
source('Covariates.R')
source('IPW.R')
source('match.R')
source('CausalForest.R')

corn_df <- read_parquet(paste(current_path, '/../data/iGIS/d_igis12.parquet', sep = ''))

corn_df$treatment <- NaN
soybeans_ind <- corn_df$prioryr_crop == 'Soybeans'
corn_ind <- corn_df$prioryr_crop == 'Corn'
corn_df$treatment[soybeans_ind] <- 1
corn_df$treatment[corn_ind] <- 0
data_df_st<-corn_df[!is.na(corn_df$treatment),]
data_df_st<-data_df_st[!is.na(data_df_st$corn_yield),]
data_df_st$fyear = as.factor(data_df_st$year)

data_df_st <- data_df_st[!is.na(data_df_st$vpdmax_4),]
data_df_st <- as.data.frame(data_df_st)

n_all <- dim(data_df_st)[1]

normalize_col <- c(5:10, 15:64, 80:94)

re_mean <- apply(data_df_st[,normalize_col], MARGIN = 2, FUN = mean)
re_sd <- apply(data_df_st[,normalize_col], MARGIN = 2, FUN = sd)
mean_mat <- matrix(rep(re_mean, each = n_all), nrow = n_all)
sd_mat <- matrix(rep(re_sd, each = n_all), nrow = n_all)

data_df_st[, normalize_col] <- (data_df_st[, normalize_col] - mean_mat)/sd_mat
data_df_st <- as.data.frame(data_df_st)

data_df_st$tile_field_ID = paste(data_df_st[,1], data_df_st[,2], sep="")

#choose number of tile/field
set.seed(94703)
m <- 5e3
#randomly select tile
selected_tile <- sample(row.names(table(data_df_st$tile_field_ID)), m)
data_df_tmp <- data_df_st[data_df_st$tile_field_ID %in% selected_tile,]
data_df_tmp$ID <- 1:dim(data_df_tmp)[1]

data_df <- add.covariates(data_df_tmp)
data_df$yield <- data_df$corn_yield 
rm(corn_df)
rm(data_df_st)
# saveRDS(data_df, file = 'data_df_0917.rds')
# data_df <- readRDS('data_df_0917.rds')

########Algorithm begin######
##############PCA############
result.pca <- list()
k <- 10
result.pca[[1]] <- GAM_function(data_df, method = 'GAM', k = k)
result.pca[[2]] <- GAM_function(data_df, method = 'IPW', k = k)
result.pca[[3]] <- GAM_function(data_df, method = 'match', k = k)
result.pca[[4]] <- CausalForest(data_df, tr = 4000)


###step 2: get ate and s.e.
res.all.pca <- list()

res.all.pca[[1]] <- ate.GAM.object(result.pca[[1]], method = 'GAM')
res.all.pca[[2]] <- ate.GAM.object(result.pca[[2]], method = 'IPW')
res.all.pca[[3]] <- ate.GAM.object(result.pca[[3]], method = 'match')
res.all.pca[[4]] <- ate.CFF.abject(result.pca[[4]])

# result.pca <- readRDS('result_pca_new_1-4.rds')
# res.all.pca <- readRDS('result_pca_new_all.rds')
# 
# saveRDS(result.pca, file = 'result_pca_new_1-3-0918.rds')
# saveRDS(res.all.pca, file = 'result_pca_new_all.rds')

##step 3: got o plots_new for all plots

