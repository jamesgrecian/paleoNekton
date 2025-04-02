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
centroid_now <- x |> 
  rasterToPoints() |> 
  as_tibble() |> 
  pivot_longer(3:30, names_to = "species", values_to = "preds") |>
  filter(preds > 0) |>
  group_by(species) |>
  mutate(centroid_x = mean(x),
         centroid_y = mean(y)) |>
  slice(1) |>
  ungroup()
  
centroid_now <- centroid_now |>
  dplyr::select("species", "centroid_x", "centroid_y")


# raster thresholds for LGM
paleo_files <- list.files("data/", pattern = "paleo_raster", full.names = T)
paleo_stack <- lapply(paleo_files, readRDS) 
paleo_stack_threshold <- paleo_stack

# threshold each species and GCM combo
for(i in 1:length(thresholds)){
  paleo_stack_threshold[[i]] <- paleo_stack[[i]] > thresholds[i]
}

# convert stack of rasters and GCMS to a tibble
paleo_centroids <- lapply(paleo_stack_threshold, function(x){ x |> rasterToPoints() |> as_tibble() |> 
    pivot_longer(3:7, names_to = "model", values_to = "preds") |>
    filter(preds > 0) |>
    group_by(model) |>
    mutate(centroid_x = mean(x),
           centroid_y = mean(y)) |>
    slice(1) |>
    ungroup()
}
)

paleo_centroids <- paleo_centroids |>
  tibble() |>
  mutate(species = names(x)) |>
  unnest(cols = c(paleo_centroids))

paleo_centroids <- paleo_centroids |>
  dplyr::select("species", "model", "centroid_x", "centroid_y")

centroids <- centroid_now |> left_join(paleo_centroids, by = "species")

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
  geom_point(aes(x = bearing, y = distance/1000, colour = guild), data = centroids) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks =  seq(-180, 180, 90)) +
  coord_polar(start = pi) +
  facet_wrap(~model)





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
