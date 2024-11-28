######################################################
### Ensemble predictions for PMIP4/CMIP6 scenarios ###
######################################################

# this function generates output predictions for LGM

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

predictLGM <- function(model_df, folds, out_path){
  
  # identify which covariates required using fitted model
  model_covs <- names(model_df$workflows[[1]]$pre$mold$predictors)
  
  # paths to LGM covariate stacks from Rahul
  lgm_fn <- list.files("data/", full.names = T)
  
  model_list <- c("AWI-ESM-1-1-LR", "CESM2-FV2", "CESM2-WACCM-FV2", "MIROC-ES2L", "MPI-ESM1-2-LR")
  
  # Loop through each LGM covariate stack
  # Save only the ensemble mean (across models and folds) for each LGM
  
  output_raster <- stack()
  
  for(i in 1:length(model_list)){
    
    cat(paste0("\nprocessing ", model_list[i], "..."))
    
    # load lgm covariates
    covs <- readRDS(lgm_fn[i])
    covs <- subset(covs, model_covs) # subset to only covariates in model formula
    
    ####################################
    ### Generate spatial predictions ###
    ####################################
    
    # generate predictions
    # https://github.com/tidymodels/planning/issues/26
    
    ### now need to generate 10 spatial predictions
    ### could simply pass the covariate data to the predict function
    ### or should this be some other weird explainer thing?
    explainer_generator <- function(input_model, input_data) {
      pred <-
        input_data |> analysis() |> st_drop_geometry() |> dplyr::select(all_of(model_covs))
      resp <-
        input_data |> analysis() |> st_drop_geometry() |> pull(PresAbs)
      explainer <-
        explain_tidymodels(input_model, data = pred, y = resp, verbose = F)
      return(explainer)
    }
    
    # create explainers for each model type
    cat("\n...processing glm explainers")
    glm_explainers <- map2(model_df$workflows[1:10], folds$splits, explainer_generator)
    cat("\n...processing gam explainers")
    gam_explainers <- map2(model_df$workflows[11:20], folds$splits, explainer_generator)
    cat("\n...processing rf explainers")
    rf_explainers <- map2(model_df$workflows[21:30], folds$splits, explainer_generator)
    cat("\n...processing brt explainers")
    brt_explainers <- map2(model_df$workflows[31:40], folds$splits, explainer_generator)
    
    # generate predictions for covariate stack
    cat("\n...generating ensemble predictions")
    glm_preds_stack <- lapply(glm_explainers, function(x) terra::predict(covs, x)) |> raster::stack()
    gam_preds_stack <- lapply(gam_explainers, function(x) terra::predict(covs, x)) |> raster::stack()
    rf_preds_stack <- lapply(rf_explainers, function(x) terra::predict(covs, x)) |> raster::stack()
    brt_preds_stack <- lapply(brt_explainers, function(x) terra::predict(covs, x)) |> raster::stack()
    
    # need final workflows from each model saved so that weights can be calculated...
    glm_auc_weights <- model_df |> filter(model == "GLM") |> pull(AUC)
    gam_auc_weights <- model_df |> filter(model == "GAM") |> pull(AUC)
    rf_auc_weights <- model_df |> filter(model == "RF") |> pull(AUC)
    brt_auc_weights <- model_df |> filter(model == "BRT") |> pull(AUC)
    
    # Create AUC weighted average of each model for the ensemble
    glm_preds_raster <- weighted.mean(glm_preds_stack, w = glm_auc_weights)
    gam_preds_raster <- weighted.mean(gam_preds_stack, w = gam_auc_weights)
    rf_preds_raster <- weighted.mean(rf_preds_stack, w = rf_auc_weights)
    brt_preds_raster <- weighted.mean(brt_preds_stack, w = brt_auc_weights)
    
    # ensemble mean
    ensemble_stack <- stack(glm_preds_raster,
                            gam_preds_raster,
                            rf_preds_raster,
                            brt_preds_raster)
    
    auc_w <- c(mean(glm_auc_weights),
               mean(gam_auc_weights),
               mean(rf_auc_weights),
               mean(brt_auc_weights))
    
    ensemble_raster <- weighted.mean(ensemble_stack, w = auc_w)
    names(ensemble_raster) <- model_list[i]
    
    output_raster <- stack(output_raster, ensemble_raster)
  }
  
  saveRDS(output_raster, out_path)
  
}

# ends