###########################
### Polar Mask function ###
###########################

# 2022-04-26

# Function based on Claus Wilke code
# https://gist.github.com/clauswilke/783e1a8ee3233775c9c3b8bfe531e28a

# create circular buffer mask to make ggplot maps look prettier

polar_mask <- function(radius_size = 5750000){
  
  # define projection
  crs_polar <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" 
  # projection outline
  polar_outline <- st_buffer(st_point(x = c(0, 0)), dist = radius_size) |>
    st_sfc(crs = crs_polar)
  
  # bounding box in transformed coordinates
  xlim <- c(-6500000, 6500000)
  ylim <- c(-6500000, 6500000)
  polar_bbox <- 
    list(
      cbind(
        c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]), 
        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
      )
    ) |>
    st_polygon() |>
    st_sfc(crs = crs_polar)
  
  # area outside the earth outline
  polar_without <- st_difference(polar_bbox, polar_outline)
  
  # return object
  p_mask <- list(
    buffer = polar_outline, # buffer shape
    mask = polar_without # mask shape
  )
  return(p_mask)
  
}
