#################################################################
### Full data Bray Curtis and clustering calculation at 25 km ###
#################################################################
#
# 2026-05-22
#
# Compute the definitive 25 km hierarchical clustering structure for the
# contemporary Southern Ocean micronekton dataset (n = 116k cells) using Bray-
# Curtis dissimilarity. This produces the hc object that downstream k-selection
# (bootstrap ROC and silhouette) and analogue matching steps will rely on.
#
# At 25 km resolution this is impossible to run on local machine
# the distance object is ~54GB.
#
# Instead using a Google cloud instance e2-highmem-32 (32 vCPUs, 256GB RAM)
# 
# Also at 25 km analogue::distance() will hit R's .C long vector limit of ~65k cells.
# Instead parallelDist::parDist sidesteps the .C path entirely (uses Rcpp), and
# fastcluster::hclust supports long vectors where base hclust() does not.
# 
# Together they make full-data clustering tractable at this resolution given
# enough RAM.
#
# Don't save the bc object - where can you open a ~54GB object in local RAM
#
# Google cloud will create both objects in under 3 min each
#

# libraries
library(tidyverse)
library(parallelDist)
library(fastcluster)
library(tictoc)

# load data
x_df <- readRDS("/home/rstudio/x_df_25km.rds")

tic()
# parallel distance calculation
bc <- parallelDist::parDist(as.matrix(x_df[,3:30]), method = "bray", threads = 10)
toc()

tic()
# hclust but with support for long vectors
hc <- fastcluster::hclust(bc, method = "ward.D2")
toc()

plot(hc)
saveRDS(hc, "hc_fastcluster_25km.rds")

# ends
