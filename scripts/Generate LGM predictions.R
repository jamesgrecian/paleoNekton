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

predict_LGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_01_v2.rds"),
            folds = folds_all[[1]],
            out_path = "data/paleo_raster_01_lgm.rds")

predict_LGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_02_v2.rds"),
            folds = folds_all[[2]],
            out_path = "data/paleo_raster_02_lgm.rds")

predict_LGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_03_v2.rds"),
            folds = folds_all[[3]],
            out_path = "data/paleo_raster_03_lgm.rds")

predict_LGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_04_v2.rds"),
            folds = folds_all[[4]],
            out_path = "data/paleo_raster_04_lgm.rds")

predict_LGM(model_df = readRDS("~/nektonAES/data/ensemble model outputs/model_df_05_v2.rds"),
            folds = folds_all[[5]],
            out_path = "data/paleo_raster_05_lgm.rds")
