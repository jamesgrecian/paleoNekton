###################################
### Process LGM data from Rahul ###
###################################

# 2024-11-27

# need to load in the models from the nektonAES folder
# then generate predictions from the paleo model outputs

# 5 GCMs
# create a raster stack for each...
# save these as rds files

# then should 'just' be able to generate the predictions using the same function as before....

require(tidyverse)
require(sf)
sf::sf_use_s2(FALSE)
require(raster)
require(ncdf4)

# custom helper function
source("R/processLGM.R")

# process files for each GCM in turn
covs_lgm_1 <- processLGM(model = "AWI-ESM-1-1-LR")
covs_lgm_2 <- processLGM(model = "CESM2-FV2")
covs_lgm_3 <- processLGM(model = "CESM2-WACCM-FV2")
covs_lgm_4 <- processLGM(model = "MIROC-ES2L")
covs_lgm_5 <- processLGM(model = "MPI-ESM1-2-LR")

saveRDS(covs_lgm_1, "data/covariate_stack_lgm_AWI-ESM-1-1-LR.rds")
saveRDS(covs_lgm_2, "data/covariate_stack_lgm_CESM2-FV2.rds")
saveRDS(covs_lgm_3, "data/covariate_stack_lgm_CESM2-WACCM-FV2.rds")
saveRDS(covs_lgm_4, "data/covariate_stack_lgm_MIROC-ES2L.rds")
saveRDS(covs_lgm_5, "data/covariate_stack_lgm_MPI-ESM1-2-LR.rds")

# ends

# legacy code below 
# code explores initial issue with sea ice extent

# define a list of models
# models <- c("AWI-ESM-1-1-LR", "CESM2-FV2", "CESM2-WACCM-FV2", "MIROC-ES2L", "MPI-ESM1-2-LR")
# 
# sic_stack <- stack()
# 
# for (i in 1:5){
#   # define pathnames for model covariates
#   fn <- list.files("data/LGM covariates", full.names = T, pattern = models[i])
#   
#   sic <- ncdf4::nc_open(fn[2])
#   tmp_array <- ncvar_get(sic, "siconc")
#   e <- ncvar_get(sic, "easting")
#   n <- ncvar_get(sic, "northing")
#   coords <- expand.grid(e,n)
#   d1 <- as.vector(tmp_array[,,1])
#   d2 <- as.vector(tmp_array[,,2])
#   d3 <- as.vector(tmp_array[,,3])
#   d4 <- as.vector(tmp_array[,,4])
#   d5 <- as.vector(tmp_array[,,5])
#   d6 <- as.vector(tmp_array[,,6])
# 
#   foo <- rasterFromXYZ(cbind(coords, d1, d2, d3, d4, d5, d6))
#   projection(foo) <- CRS(prj)
#   foo <- extend(foo, r, value = 0)
#   foo <- resample(foo, r, method="bilinear")
#   foo <- mask(foo, mask = subset(r, 1), maskvalue = NA)
#   foo <- calc(foo, mean)
#   names(foo) <- models[i]
#   
#   sic_stack <- stack(sic_stack, foo)
# }
# 
#   
#   
# plot(sic_stack)
# 
# sic_stack_df <- sic_stack |> 
#   rasterToPoints() |> 
#   as_tibble() |>
#   pivot_longer(3:7, names_to = "model", values_to = "concentration")
# 
# 
# require(sf)
# sf::sf_use_s2(FALSE)
# 
# source("~/nektonAES/R/polar_mask.R")
# polar_buffer <- polar_mask(radius_size = 5750000)
# 
# # load in shapefile for background mask, clip and project
# world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
# CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) |> sf::st_as_sfc()
# 
# world_shp <- world_shp |> sf::st_crop(CP)
# world_shp <- world_shp |> sf::st_transform(prj)
# world_shp <- world_shp |> st_union()
# 
# p1 <- ggplot() +
#   theme_void(base_size = 12) +
#   geom_raster(aes(x = x, y = y, fill = concentration), data = sic_stack_df) +
#   geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") + 
#   geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = NA, size = NA) +
#   coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = prj, ndiscr = 1000) +
#   facet_wrap(~model) +
#   scale_fill_binned("", 
#                     type = "viridis",
#                     n.breaks = 10,
#                     limits = c(0, 100),
#                     guide = guide_colourbar(
#                       barwidth = 1,
#                       barheight = 15))
# 
# 
# ggsave(filename = "lgm sea ice.jpeg", 
#        plot = p1,
#        width = 7,
#        height = 5,
#        dpi = 500)
# 
# 
# # ends