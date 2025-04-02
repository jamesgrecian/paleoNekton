########################
### Raster threshold ###
########################

# helper function to calculate Youden's index
# use for thresholding rasters based on observations
# and model input raster (ensemble average)
# really only need to do this once and save the thresholds...

# function to calculate raster thresholds from ROC curves
raster_threshold <- function(input_raster_path, sp){
  
  # load raster
  input_raster <- readRDS(input_raster_path)
  
  # load original data
  dat <- readRDS("~/nektonAES/data/presence_absence_data_10k_with_covariates_2024-06-17.rds")
  dat <- dat |> filter(species == unique(species)[sp]) # focal species
  
  # extract probabilities from raster and append to data
  est_prob <- raster::extract(input_raster, cbind(dat$x, dat$y))
  
  roc_obj <- pROC::roc(dat$PresAbs, est_prob)   # calculate ROC curve
  roc_threshold <- pROC::coords(roc_obj, x = "best")$threshold  # pull out threshold based on Youden
  
  return(roc_threshold)
}

# ends
