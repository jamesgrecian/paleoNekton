####################################################
### Ensemble predictions for PMIP4 LGM scenarios ###
####################################################

# 2024-11-27

# This is based on a revised version of the futures predictions from nektonAES

# There are 28 species
# There are 5 GCMs for the LGM
# Each ensemble model is 4 models x 10 folds
# (4 x 10) x 5 = 200 surfaces per species...

# Calculate LGM suface for each fold and each model
# Don't present individual folds, average over all
# Also average over ensemble of GLM, GAM, RF and BRT models
# Then each e-sdm is only one field
# Generate that field for each GCM
# Save a raster stack with the LGM prediction for each of 5 GCMs, and an ensemble

# function operates the same as the general prediction function from nektonAES
# but no plotting
# input is:
# the model_df to scrape which covariates are needed for each species
# folds subset for focal species
# file path and name to output the LGM raster prediction stack

# load libraries
require(tidyverse)
require(raster)
require(tidymodels)
require(sf)
sf::sf_use_s2(FALSE)
require(spatialsample)
require(DALEXtra)

#Mikes hack function
quantile.hardhat_importance_weights <- \(x, ...) rep(NA, length(x))

# helper function
source("R/predictLGM.R")

# load folds data
folds_all <- readRDS("~/nektonAES/data/folds_weights_v2.rds")

# use custom function to generate a prediction stack for LGM
# containing estimated distributions for each of the 5 GCMs from the ensemble model
# don't store the individual folds or individual SDM models

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_01_v2.rds"),
            folds = folds_all[[1]],
            out_path = "data/paleo_raster_01_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_02_v2.rds"),
            folds = folds_all[[2]],
            out_path = "data/paleo_raster_02_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_03_v2.rds"),
            folds = folds_all[[3]],
            out_path = "data/paleo_raster_03_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_04_v2.rds"),
            folds = folds_all[[4]],
            out_path = "data/paleo_raster_04_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_05_v2.rds"),
            folds = folds_all[[5]],
            out_path = "data/paleo_raster_05_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_06_v2.rds"),
           folds = folds_all[[6]],
           out_path = "data/paleo_raster_06_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_07_v2.rds"),
           folds = folds_all[[7]],
           out_path = "data/paleo_raster_07_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_08_v2.rds"),
           folds = folds_all[[8]],
           out_path = "data/paleo_raster_08_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_09_v2.rds"),
           folds = folds_all[[9]],
           out_path = "data/paleo_raster_09_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_10_v2.rds"),
           folds = folds_all[[10]],
           out_path = "data/paleo_raster_10_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_11_v2.rds"),
           folds = folds_all[[11]],
           out_path = "data/paleo_raster_11_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_12_v2.rds"),
           folds = folds_all[[12]],
           out_path = "data/paleo_raster_12_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_13_v2.rds"),
           folds = folds_all[[13]],
           out_path = "data/paleo_raster_13_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_14_v2.rds"),
           folds = folds_all[[14]],
           out_path = "data/paleo_raster_14_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_15_v2.rds"),
           folds = folds_all[[15]],
           out_path = "data/paleo_raster_15_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_16_v2.rds"),
           folds = folds_all[[16]],
           out_path = "data/paleo_raster_16_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_17_v2.rds"),
           folds = folds_all[[17]],
           out_path = "data/paleo_raster_17_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_18_v2.rds"),
           folds = folds_all[[18]],
           out_path = "data/paleo_raster_18_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_19_v2.rds"),
           folds = folds_all[[19]],
           out_path = "data/paleo_raster_19_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_20_v2.rds"),
           folds = folds_all[[20]],
           out_path = "data/paleo_raster_20_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_21_v2.rds"),
           folds = folds_all[[21]],
           out_path = "data/paleo_raster_21_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_22_v2.rds"),
           folds = folds_all[[22]],
           out_path = "data/paleo_raster_22_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_23_v2.rds"),
           folds = folds_all[[23]],
           out_path = "data/paleo_raster_23_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_24_v2.rds"),
           folds = folds_all[[24]],
           out_path = "data/paleo_raster_24_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_25_v2.rds"),
           folds = folds_all[[25]],
           out_path = "data/paleo_raster_25_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_26_v2.rds"),
           folds = folds_all[[26]],
           out_path = "data/paleo_raster_26_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_27_v2.rds"),
           folds = folds_all[[27]],
           out_path = "data/paleo_raster_27_lgm.rds")

predictLGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_28_v2.rds"),
           folds = folds_all[[28]],
           out_path = "data/paleo_raster_28_lgm.rds")

# ends
