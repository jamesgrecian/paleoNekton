##################################################################
### How to assess movement in species distributions at the LGM ###
##################################################################

# most dicussions around climate velocities seem to be focused on climate/ environment
# how to assess changes in distributions (in response to climate) in this paradigm

# can we calculate the distance and bearing
# from contemporary distributions to LGM?

# 2025-04-02

# need contemporary distribution
# need historical (or future)
# calculate centroid of contemporary
# calculate centroid of past
# bearing and distance between the two...?
require(tidyverse)
require(raster)
require(circular)
require(patchwork)
source("R/raster_threshold.R")

### Load and process contemporary data ###
files_now <- list.files("~/nektonAES/data/ensemble model outputs", pattern = "ensemble", full.names = T)
files_now <- files_now[grep("v2", files_now)]
stack_now <- lapply(files_now, readRDS) |> stack()
names(stack_now) <- c("species01", "species02", "species03", "species04", "species05",
                      "species06", "species07", "species08", "species09", "species10",
                      "species11", "species12", "species13", "species14", "species15",
                      "species16", "species17", "species18", "species19", "species20",
                      "species21", "species22", "species23", "species24", "species25",
                      "species26", "species27", "species28")

# apply function to each species
thresholds <- map2(files_now, 1:28, raster_threshold)
thresholds <- thresholds |> unlist()

# create binary raster stack based on species specific thresholding
x <- stack_now > thresholds

# calculate the distance to the south pole for each raster cell
distance_now <- x |> 
  rasterToPoints() |> 
  as_tibble() |> 
  pivot_longer(3:30, names_to = "species", values_to = "preds") |>
  filter(preds > 0) |>
  mutate(distance = sqrt(x^2 + y^2),
         bearing = atan2(y, x) * (180 / pi)) |>
  group_by(species) |>
  summarise(distance = mean(distance)/1000,
            bearing = mean(circular(bearing, units = "degrees", template = "geographics")))

# raster thresholds for LGM
paleo_files <- list.files("data/", pattern = "paleo_raster", full.names = T)
paleo_stack <- lapply(paleo_files, readRDS) 
paleo_stack_threshold <- paleo_stack

# threshold each species and GCM combo
for(i in 1:length(thresholds)){
  paleo_stack_threshold[[i]] <- paleo_stack[[i]] > thresholds[i]
}

# convert stack of rasters and GCMS to a tibble
distance_paleo <- lapply(paleo_stack_threshold, function(x){ x |> rasterToPoints() |> as_tibble() |> 
    pivot_longer(3:7, names_to = "model", values_to = "preds") |>
    filter(preds > 0) |>
    group_by(model)
}
)

distance_paleo <- distance_paleo |>
  tibble() |>
  mutate(species = names(x)) |>
  unnest(cols = c(distance_paleo))

distance_paleo <- distance_paleo |>
  mutate(distance = sqrt(x^2 + y^2),
         bearing = atan2(y, x) * (180 / pi)) |>
  group_by(species, model) |>
  summarise(distance = mean(distance)/1000,
            bearing = mean(circular(bearing, units = "degrees", template = "geographics")))



distance_shift <- distance_paleo |> left_join(distance_now, by = "species")

distance_shift <- distance_shift |> mutate(shift_d = distance.x - distance.y,
                                           shift_b = (bearing.x - bearing.y + 180) %% 360 - 180)

distance_shift


diff_deg <- as.numeric(diff)
if (diff_deg > 180) {
  diff_deg <- diff_deg - 360
} else if (diff_deg < -180) {
  diff_deg <- diff_deg + 360
}


distance_shift <- distance_shift |>
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                        "species05", "species06", "species07", "species08",
                                        "species09", "species10", "species11", "species12",
                                        "species13") ~ "Fish",
                         species == "species14" ~ "Krill",
                         species %in% c("species15", "species16", "species17", "species18",
                                        "species19", "species20", "species21", "species22",
                                        "species23", "species24", "species25", "species26", 
                                        "species27", "species28") ~ "Squid"))

distance_shift <- distance_shift |> mutate(species = factor(species))
species_names <- c("B. antarcticus", "N. coatsorum", "P. antarctica", "E. antarctica", "E. carlsbergi", "G. bolini", "G. braueri",
                   "G. fraseri", "G. nicholsi", "G. opisthopterus", "K. anderssoni", "P. bolini", "P. tenisoni", "E. superba",
                   "A. antarcticus", "B. abyssicola", "G. glacialis", "G. antarcticus", "H. atlantica", "H. eltaninae", "K. longimana",
                   "M. hyadesi", "M. hamiltoni", "M. ingens", "M. robsoni", "P. glacialis", "S. circumantarctica", "T. filippovae")     

levels(distance_shift$species) <- species_names


p1 <- ggplot() + 
  theme_minimal(base_size = 8) +
  geom_boxplot(aes(x = shift_d, y = species), colour = "#D3DDDC", fill = "#D3DDDC", 
               width = .4, outliers = F,
               data = distance_shift) +
  geom_point(aes(x = shift_d, y = species), 
               data = distance_shift, size = .25) +
  scale_x_continuous(name = "Range shift (km)", limits = c(-1000, 2500)) +
  scale_y_discrete(name = NULL, limits = rev) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(panel.border = element_rect(color = "gray 50", fill = NA)) +
  facet_grid(rows = vars(guild), scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_text(face = "italic"))


ggplot() + 
  theme_minimal(base_size = 8) +
  geom_histogram(aes(x = shift_b), 
             data = distance_shift) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks =  seq(-180, 180, 90)) +
  coord_polar(start = pi) +
  facet_wrap(~guild)

    geom_boxplot(aes(x = shift_d, y = species), colour = "#D3DDDC", fill = "#D3DDDC", 
               width = .4, outliers = F,
               data = distance_shift) +
  geom_point(aes(x = shift_d, y = species), 
             data = distance_shift, size = .25) +
  scale_x_continuous(name = "Range shift (km)", limits = c(-1000, 2500)) +
  scale_y_discrete(name = NULL, limits = rev) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(panel.border = element_rect(color = "gray 50", fill = NA)) +
  facet_grid(rows = vars(guild), scales = "free_y", space = "free_y") +
  theme(axis.text.y = element_text(face = "italic"))




ggsave(filename = "LGM shift boxplot.jpeg",
       plot = p1,
       path = "figures/",
       width = 140,
       height = 140,
       units = "mm",
       dpi = 500)



paleo_centroids <- paleo_centroids |>
  tibble() |>
  mutate(species = names(x)) |>
  unnest(cols = c(paleo_centroids))

paleo_centroids <- paleo_centroids |>
  dplyr::select("species", "model", "centroid_x", "centroid_y")

centroid_now[c("lon_now", "lat_now")] <- centroid_now |>
  st_as_sf(coords = c("centroid_x", "centroid_y")) |>
  sf::st_set_crs(projection(stack_now)) |>
  st_transform(4326) |>
  st_coordinates()

paleo_centroids[c("lon_paleo", "lat_paleo")] <- paleo_centroids |>
  st_as_sf(coords = c("centroid_x", "centroid_y")) |>
  sf::st_set_crs(projection(stack_now)) |>
  st_transform(4326) |>
  st_coordinates()


centroids <- centroid_now |> left_join(paleo_centroids, by = "species")


geo_dist <- geosphere::distGeo(p1 = centroids[ , c("lon_now", "lat_now")], p2 = centroids[ , c("lon_paleo", "lat_paleo")])/1000
centroids$alt_dist <- geo_dist
geo_bearing <- geosphere::bearing(p1 = centroids[ , c("lon_now", "lat_now")], p2 = centroids[ , c("lon_paleo", "lat_paleo")])
centroids$alt_bearing <- geo_bearing

centroids <- centroids |> mutate(distance = sqrt(((centroid_x.y - centroid_x.x)^2) + ((centroid_y.y - centroid_y.x)^2)),
                    bearing = atan2(centroid_y.y - centroid_y.x, centroid_x.y - centroid_x.x)* (180 / pi))

centroids <- centroids |> 
  mutate(guild = case_when(species %in% c("species01", "species02", "species03", "species04",
                                          "species05", "species06", "species07", "species08",
                                          "species09", "species10", "species11", "species12",
                                          "species13") ~ "Fish",
                           species == "species14" ~ "Krill",
                           species %in% c("species15", "species16", "species17", "species18",
                                          "species19", "species20", "species21", "species22",
                                          "species23", "species24", "species25", "species26", 
                                          "species27", "species28") ~ "Squid"))

ggplot() +
  theme_minimal() +
  geom_point(aes(x = alt_bearing, y = alt_dist, colour = guild), data = centroids) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks =  seq(-180, 180, 90)) +
  coord_polar(start = pi) +
  facet_wrap(~model)

p1 <- ggplot() + 
  ggtitle("cartesian bearing") +
  theme_minimal() +
  geom_histogram(aes(x = bearing), bins = 50, data = centroids) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks =  seq(-180, 180, 90)) +
  coord_polar(start = pi) +
  facet_wrap(~model)

p2 <- ggplot() + 
  ggtitle("geographic bearing") +
  theme_minimal() +
  geom_histogram(aes(x = alt_bearing), bins = 50, data = centroids) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks =  seq(-180, 180, 90)) +
  coord_polar(start = pi) +
  facet_wrap(~model)

p1 + p2


species_names <- c("B. antarcticus", "N. coatsorum", "P. antarctica", "E. antarctica", "E. carlsbergi", "G. bolini", "G. braueri",
                   "G. fraseri", "G. nicholsi", "G. opisthopterus", "K. anderssoni", "P. bolini", "P. tenisoni", "E. superba",
                   "A. antarcticus", "B. abyssicola", "G. glacialis", "G. antarcticus", "H. atlantica", "H. eltaninae", "K. longimana",
                   "M. hyadesi", "M. hamiltoni", "M. ingens", "M. robsoni", "P. glacialis", "S. circumantarctica", "T. filippovae")     

distance_now <- distance_now |> mutate(species = factor(species))
range_shift_245 <- range_shift_245 |> mutate(species = factor(species))
range_shift_585 <- range_shift_585 |> mutate(species = factor(species))

levels(distance_now$species) <- species_names
levels(range_shift_245$species) <- species_names
levels(range_shift_585$species) <- species_names

range_shift_245 <- range_shift_245 |> mutate(scenario = "SSP2-4.5")
range_shift_585 <- range_shift_585 |> mutate(scenario = "SSP5-8.5")





x
# Shannon diversity function
shannon_fun <- function(x) {
  if (all(is.na(x))) return(NA)
  total <- sum(x, na.rm = TRUE)
  if (total == 0) return(0)
  p_norm <- x / total
  p_nonzero <- p_norm[p_norm > 0]
  -sum(p_nonzero * log(p_nonzero))
}

# Apply function over the stack
shannon_raster <- calc(stack_now, fun = shannon_fun)

# Simpson diversity function
simpson_fun <- function(x) {
  if (all(is.na(x))) return(NA)
  total <- sum(x, na.rm = TRUE)
  if (total == 0) return(0)
  p_norm <- x / total
  return(1 - sum(p_norm^2, na.rm = TRUE))
}

# Apply to all raster cells
simpson_raster <- calc(stack_now, fun = simpson_fun)

plot(simpson_raster)


diversity_df <- x |> rasterToPoints() |> as_tibble()

diversity_df$H <- vegan::diversity(diversity_df[,3:30])
diversity_df$simp <- vegan::diversity(diversity_df[,3:30], "simpson")


ggplot() + 
  geom_raster(aes(x, y, fill = simp), data = diversity_df) +
  coord_fixed()



# should start simple...
# map the three guilds at the LGM in the same way as for other paper

guilds <- c(rep(1, times = 13), 2, rep(3, times = 14))

paleo_stack_1 <- lapply(paleo_stack, subset, "AWI.ESM.1.1.LR") |> stack() |> stackApply(guilds, mean)
paleo_stack_2 <- lapply(paleo_stack, subset, "CESM2.FV2") |> stack() |> stackApply(guilds, mean)
paleo_stack_3 <- lapply(paleo_stack, subset, "CESM2.WACCM.FV2") |> stack() |> stackApply(guilds, mean)
paleo_stack_4 <- lapply(paleo_stack, subset, "MIROC.ES2L") |> stack() |> stackApply(guilds, mean)
paleo_stack_5 <- lapply(paleo_stack, subset, "MPI.ESM1.2.LR") |> stack() |> stackApply(guilds, mean)

lgm_fish <- stack(subset(paleo_stack_1, 1),
                  subset(paleo_stack_2, 1),
                  subset(paleo_stack_3, 1),
                  subset(paleo_stack_4, 1),
                  subset(paleo_stack_5, 1))

lgm_krill <- stack(subset(paleo_stack_1, 2),
                  subset(paleo_stack_2, 2),
                  subset(paleo_stack_3, 2),
                  subset(paleo_stack_4, 2),
                  subset(paleo_stack_5, 2))

lgm_squid <- stack(subset(paleo_stack_1, 3),
                  subset(paleo_stack_2, 3),
                  subset(paleo_stack_3, 3),
                  subset(paleo_stack_4, 3),
                  subset(paleo_stack_5, 3))

lgm_fish <- mean(lgm_fish)
lgm_krill <- mean(lgm_krill)
lgm_squid <- mean(lgm_squid)

lgm_guilds <- stack(lgm_fish, lgm_krill, lgm_squid)
names(lgm_guilds) <- c("fish", "krill", "squid")
plot(lgm_guilds)

lgm_guilds_df <- lgm_guilds |> 
  rasterToPoints() |> 
  as_tibble() |>
  pivot_longer(3:5, names_to = "guild", values_to = "preds")

# libraries
require(sf)
sf::sf_use_s2(FALSE)
require(patchwork)

# pretty gradient and polar mask functions
source("R/discrete_gradient.R")
source("R/polar_mask.R")
polar_buffer <- polar_mask(radius_size = 5750000)

# define projection
crs_polar <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# load in shapefile for background mask, clip and project
world_shp <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf")
CP <- sf::st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 0), crs = 4326) |> sf::st_as_sfc()

world_shp <- world_shp |> sf::st_crop(CP)
world_shp <- world_shp |> sf::st_transform(crs_polar)
world_shp <- world_shp |> st_union()


p1 <- ggplot() + 
  theme_void() +
  geom_raster(aes(x = x, y = y, fill = preds), data = lgm_guilds_df) +
  geom_sf(aes(), data = world_shp, colour = "grey60", fill = "grey60") +
  geom_sf(aes(), data = polar_buffer$mask, fill = "white", color = "grey40", size = 0.5 / .pt) +
  coord_sf(xlim = c(-6400000, 6400000), ylim = c(-6400000, 6400000), expand = FALSE, crs = crs_polar, ndiscr = 1000) +
  facet_grid(rows = vars(guild), switch = "y") +
  scale_fill_discrete_gradient("Mean Habitat Importance",
                               colours = viridis::viridis(10),
                               bins = 10,
                               limits = c(0, 1),
                               breaks = seq(0, 1, 0.2),
                               labels = seq(0, 1, 0.2),
                               guide = guide_colourbar(
                                 nbin = 500,
                                 raster = T,
                                 frame.colour = "grey40",
                                 ticks.colour = "grey40",
                                 frame.linewidth = .1,
                                 barwidth = 0.5,
                                 barheight = 15,
                                 direction = "vertical",
                                 title.position = "right",
                                 title.theme = element_text(
                                   hjust = 0.5,
                                   size = 10,
                                   angle = 90))) +
  theme(strip.text = element_text(angle = 90))
                                 


