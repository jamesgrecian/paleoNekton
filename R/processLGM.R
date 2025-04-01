####################################
### Function to process LGM data ###
####################################

# 2024-11-27
# 2025-04-01
# edit to ensure output raster has same extent as contemporary covariates
# probably should clip water temperature...

# Given a GCM model name
# Load the covariates matching the model
# calculate the seasonal mean
# output a raster stack

# this may need modifying for the sea ice data...

processLGM <- function(model){
  
  # define pathnames for lgm model covariates
  fn <- list.files("data/LGM covariates", full.names = T, pattern = model)
  
  # load contemporary covariate stack
  covs <- readRDS("~/nektonAES/data/covariate_stack.rds")
  
  # define projection
  prj <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  
  # create a stack of the monthly covariate and calculate the mean...
  mld <- stack(fn[1])
  mld <- calc(mld, mean)
  
  sos <- stack(fn[3])
  sos <- calc(sos, mean)
  
  tos <- stack(fn[4])
  tos <- calc(tos, mean)
  
  zos <- stack(fn[5])
  zos <- rotate(zos) # for some reason zos is 0-360 not -180 to 180
  zos <- calc(zos, mean)
  
  r <- stack(mld, sos, tos, zos)
  names(r) <- c("mld", "sos", "tos", "zos")
  
  r <- raster::crop(r, raster::extent(-180, 180, -90, -30)) # crop to southern hemisphere
  
  r <- raster::projectRaster(r,
                             covs, # align to contemporary covariates
                             method = "bilinear",
                             crs = "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  
  r <- mask(r, subset(covs, 1)) # crop to same extent as contemporary covariates
  
  r$tos[r$tos < -1.8] <- -1.8 # force minimum temperature
  
  # sea ice concentration is a bit more complicated
  sic <- ncdf4::nc_open(fn[2])          # open ncdf
  tmp_array <- ncvar_get(sic, "siconc") # create array of sea ice concentration values
  e <- ncvar_get(sic, "easting")        # extract eastings
  n <- ncvar_get(sic, "northing")       # extract northings
  coords <- expand.grid(e,n)            # combine grid cell values together to get coordinates
  d1 <- as.vector(tmp_array[,,1])       # vectorise each sea ice month layer
  d2 <- as.vector(tmp_array[,,2])
  d3 <- as.vector(tmp_array[,,3])
  d4 <- as.vector(tmp_array[,,4])
  d5 <- as.vector(tmp_array[,,5])
  d6 <- as.vector(tmp_array[,,6])
  ncdf4::nc_close(sic)                  # close ncdf
  
  # from these objects recreate a raster stack
  sic_stack <- rasterFromXYZ(cbind(coords, d1, d2, d3, d4, d5, d6))
  projection(sic_stack) <- CRS(prj)
  sic_stack <- extend(sic_stack, r, value = 0)                      # extend to same area as other covariates
  sic_stack <- resample(sic_stack, r, method="bilinear")            # resample to the same grid resolution
  sic_stack <- mask(sic_stack, mask = subset(r, 1), maskvalue = NA) # mask the continent
  sic_stack <- calc(sic_stack, mean)                                # calculate the mean
  sic_stack[sic_stack < 0] <- 0                                     # force negative values to be 0
  names(sic_stack) <- "sic"
  
  # calculate the gradient covariates
  sst_grad <- terrain(r$tos, opt = "slope", unit = "radians")
  names(sst_grad) <- "sst_grad"
  ssh_grad <- terrain(r$zos, opt = "slope", unit = "radians")
  names(ssh_grad) <- "ssh_grad"
  
  # download some bathymetry data as an example from marmap
  # use 10 minute resolution as proxy of <0.25 degrees
  # this can be downgraded to 25 km x 25 km
  bathy <- marmap::getNOAA.bathy(lon1 = -180,
                                 lon2 = 180,
                                 lat1 = -90,
                                 lat2 = -40,
                                 resolution = 10)
  bat <- marmap::as.raster(bathy)
  bat <- projectRaster(from = bat, to = r)
  bat <- mask(bat, mask = subset(r, 1), maskvalue = NA)
  bat <- bat + 125 # sea levels were 125 metres lower at the LGM
  names(bat) <- "bat"
  
  # combine ice with the rest of the covariates in a stack and save output
  stack_out <- stack(r$tos, sst_grad, r$sos, r$zos, ssh_grad, r$mld, sic_stack, bat)
  names(stack_out) <- c("sst", "sst_grad", "sal", "ssh", "ssh_grad", "mld", "sic", "bat")
  
  return(stack_out)
}
