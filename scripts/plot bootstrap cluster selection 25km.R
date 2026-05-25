###############################################
### Plot bootstraps for k cluster selection ###
###############################################

# 2026-05-25

# Plot the output from the bootstrap k-selection script
#
# Generates three panels showing:
# 1. the bootstrap distribution of the optimal dissimilarity threshold
# 2. the AUC
# 3. the silhouette width
# All plotted across candidate cluster k values 2:28. Each panel overlays
# per-iteration points with the mean and 95% bootstrap CI.
#
# Where the optimal threshold asymptotes and the silhouette width peaks is
# the chosen k
#
# Reads either the partial output (to follow a run in progress) or the final
# output (when revisiting). Comment out whichever isn't being used.
#

# load libraries
library(tidyverse)
library(patchwork)

# If bootstrap script is running then load partial to follow progress
# If revisiting then load final
#dat <- readRDS("data/roc_silhouette_bootstrap_partial.rds")
dat <- readRDS("data/roc_silhouette_bootstrap_kselection.rds")

# Summary statistics per k cluster: mean and 95% bootstrap CI
summary_df <- dat |>
  group_by(k) |>
  summarise(
    optimal_mean = mean(optimal),
    optimal_lo   = quantile(optimal, 0.025),
    optimal_hi   = quantile(optimal, 0.975),
    auc_mean     = mean(auc),
    auc_lo       = quantile(auc, 0.025),
    auc_hi       = quantile(auc, 0.975),
    sil_mean     = mean(silhouette),
    sil_lo       = quantile(silhouette, 0.025),
    sil_hi       = quantile(silhouette, 0.975),
  ) |>
  ungroup()

# Optimal dissimilarity threshold
p1 <- ggplot() + 
  theme_bw() +
  geom_point(aes(x = k, y = optimal), data = dat, alpha = 0.05) +
  geom_ribbon(aes(x = k, ymin = optimal_lo, ymax = optimal_hi), 
              data = summary_df, alpha = 0.3, fill = "steelblue") +
  geom_line(aes(x = k, y = optimal_mean), 
            data = summary_df, color = "steelblue", linewidth = 1) +
  xlim(2, 30) +
  geom_vline(xintercept = 14) +
  labs(x = "Number of clusters (k)", y = "Optimal dissimilarity threshold")

# Area under the curve
p2 <- ggplot() + 
  theme_bw() +
  geom_point(aes(x = k, y = auc), data = dat, alpha = 0.05) +
  geom_ribbon(aes(x = k, ymin = auc_lo, ymax = auc_hi), 
              data = summary_df, alpha = 0.3, fill = "steelblue") +
  geom_line(aes(x = k, y = auc_mean), 
            data = summary_df, color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = 14) +
  xlim(2, 30) +
  labs(x = "Number of clusters (k)", y = "Area Under the Curve")

# Silhouette width
p3 <- ggplot() + 
  theme_bw() +
  geom_point(aes(x = k, y = silhouette), data = dat, alpha = 0.05) +
  geom_ribbon(aes(x = k, ymin = sil_lo, ymax = sil_hi), 
              data = summary_df, alpha = 0.3, fill = "steelblue") +
  geom_line(aes(x = k, y = sil_mean), 
            data = summary_df, color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = 14) +
  xlim(2, 30) +
  labs(x = "Number of clusters (k)", y = "Silhouette width")

# Stack vertically with patchwork
p <- p1 / p2 / p3

# output
ggsave(filename = "bootstrap cluster selection 25 km.jpeg",
       plot = p,
       path = "figures/",
       width = 200,
       height = 300,
       units = "mm",
       dpi = 300)

# ends
