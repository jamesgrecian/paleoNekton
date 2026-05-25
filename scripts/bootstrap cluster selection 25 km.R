##################################################################
### Bootstrap 25km cluster selection based on ROC + silhouette ###
##################################################################

# 2026-05-22

# Impossible to run cluster assignment on full dataset at 25 km following
# Christine's approach. Instead use the cloud cluster object and bootstrap the
# data frame.
#
# For each candidate cluster set, cut the cloud tree and characterise the bootstrap
# distribution of two independent metrics: the analogue::roc optimal dissimilarity
# threshold (used by Christine) and the mean silhouette width (used by Ryan).
#
# Find the k cluster where dissimilarity asymptotes and the silhouette width
# peaks.
#
# Each iteration of the bootstrap only recomputes Bray-Curtis on a small stratified
# sample (500 cells from each of k samples) and so the full 116k x 116k distance
# matrix is never created. This will therefore run locally on an M1 Pro with 16GB RAM.
#
# Note analogue::roc's k argument is the number of nearest analogues to compare
# per sample NOT the k clusters to cut the tree at. Original code from Christine
# based roc k on the size of the smallest cluster but this saturates when 
# subsampling and so each iteration returns the same value. Instead use default
# k = 1, so that function can calculate true AUC
#
# This loop will iteratively save as it progresses through k clusters as whole
# process takes around 24 hours -- run time increases as k increases
#

# load libraries
library(parallelDist)
library(analogue)
library(cluster)
library(tictoc)

# load 25 km resolution species dataset
x_df <- readRDS("data/x_df_25km.rds")
# load cloud computed hierarchical cluster tree
hc   <- readRDS("data/hc_fastcluster_25km.rds")

# Check minimum cluster size for each tree branch
# If size falls below 500 then cluster would be fully sampled not subsampled
min_cluster_sizes <- sapply(2:28, function(k) min(table(cutree(hc, k = k))))
plot(2:28, min_cluster_sizes, log = "y", type = "b")
abline(h = 500, lty = 2, col = "red")

################################
### Set bootstrap parameters ###
################################
k_values           <- 2:28 # number of candidate clusters to check
n_iter             <- 500  # number of iterations for the bootstrap
sample_per_cluster <- 500  # number of cells to subsample
roc_k              <- 1    # analogue::roc nearest analogues NOT cutree k
n_threads          <- 5    # threads for parDist

output_path_partial <- "data/roc_silhouette_bootstrap_partial.rds"
output_path_final   <- "data/roc_silhouette_bootstrap_kselection.rds"

#################
### Bootstrap ###
#################

# run through k 2:28
# cut cluster tree at k clusters
# sample 500 cells from each of the k clusters
# repeat 500 times for each k clusters

results <- list()

for (k_cut in k_values) {
  
  message(sprintf("Running k = %d (%d/%d)", k_cut, which(k_values == k_cut), length(k_values)))
  tic()
  
  clusters_full <- cutree(hc, k = k_cut)
  
  optimal_thresholds <- numeric(n_iter)
  auc_values         <- numeric(n_iter)
  silhouette_values  <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    
    # Stratified sample up to sample_per_cluster cells per cluster; each cell
    # retains its full-clustering label.
    idx <- unlist(lapply(unique(clusters_full), function(cl) {
      members <- which(clusters_full == cl)
      sample(members, min(sample_per_cluster, length(members)))
    }))
    
    d_sub        <- parallelDist::parDist(as.matrix(x_df[idx, 3:30]), method = "bray", threads = n_threads)
    clusters_sub <- clusters_full[idx]
    
    # Silhouette directly on the dist; ROC needs a matrix.
    sil                  <- cluster::silhouette(clusters_sub, d_sub)
    silhouette_values[i] <- mean(sil[, "sil_width"])
    
    p_roc                 <- analogue::roc(as.matrix(d_sub), groups = clusters_sub, k = roc_k)
    optimal_thresholds[i] <- p_roc$roc$Combined$optimal
    
    # Per-group AUCs only; Combined AUC is NA when classification saturates.
    n_groups      <- nrow(p_roc$statistics) - 1L
    auc_values[i] <- mean(p_roc$statistics$AUC[seq_len(n_groups)], na.rm = TRUE)
  }
  
  results[[as.character(k_cut)]] <- data.frame(
    k          = k_cut,
    iteration  = seq_len(n_iter),
    optimal    = optimal_thresholds,
    auc        = auc_values,
    silhouette = silhouette_values
  )
  
  toc()
  
  # Progressive save: if R crashes mid-sweep, only the in-progress k is lost.
  saveRDS(do.call(rbind, results), output_path_partial)
}


# Final save 
results_df <- do.call(rbind, results)
saveRDS(results_df, output_path_final)

# ends
