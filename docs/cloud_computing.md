# Cloud computing

This project required a cloud compute step for the contemporary bioregionalisation
at 25 km resolution. Everything else — the bootstrap k-selection, analogue matching
to paleo grids, visualisation — runs locally on an M1 Pro (16 GB RAM).

This document covers what the cloud step was, why it exists, how to set one up,
and what comes out of it.

## Why cloud

The contemporary Southern Ocean micronekton dataset at 25 km has 116,341 grid
cells across 28 species. Computing the full Bray-Curtis dissimilarity matrix and 
Ward.D2 hierarchical clustering on this dataset needs two things that aren't
simultaneously available locally:

**A C++/Rcpp-backed distance function.** `analogue::distance()` — the function 
in Christine's original methodology — calls R's `.C` interface, which is subject
to a long-vector limit at ~65k cells. At 25 km resolution we're well over this,
so the function errors out regardless of RAM.

`parallelDist::parDist` uses Rcpp instead and produces a bit-identical result for
Bray-Curtis without the limit. We verified this at 50 km, where both functions
can be run for comparison.

**Enough RAM to hold the distance object plus working memory.** The Bray-Curtis
dist object at 25 km is ~54 GB just for the values (116,341 choose 2 × 8 bytes).
Peak RAM during the `fastcluster::hclust` step climbs to ~90 GB.

The intermediate alternative — aggregating to 50 km — was the original workaround
and produces ~29k cells, which fits under the `.C` limit and into local RAM. The
downstream paleo deposit matching, however, is more informative at 25 km when zooming in
on snow petrel foraging ranges. Swapping to `parallelDist` + cloud RAM made this
feasible. So: 25 km, with the trade that the full-data clustering step runs on GCP.

The bootstrap k-selection step does *not* need the cloud. It subsamples per
iteration (~ k × 500 cells) and never materialises the full distance matrix.

## Instance specification

- **Name:** `bioregion-25km`
- **Region:** `europe-west2` (London)
- **Machine type:** `e2-highmem-32` (32 vCPUs, 256 GB RAM)
- **Boot disk:** Ubuntu 22.04 LTS, ≥ 200 GB
- **Cost:** ~£1.30/hr when running; <10p/day when stopped but holding memory

Stop the instance when not in use — packages, data, and the home directory all persist across stop/start cycles.

## Initial setup

Create an instance in the GCP console (Compute Engine → VM Instances → Create), 
SSH in via the browser SSH button. The SSH user has sudo; the RStudio login user does not,
so all system-level installs happen here.

```bash
# R and geospatial system libraries
sudo apt update
sudo apt install -y r-base r-base-dev libproj-dev libgdal-dev libgeos-dev

# RStudio Server
sudo apt install -y gdebi-core
wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.04.2-764-amd64.deb
sudo gdebi rstudio-server-2024.04.2-764-amd64.deb

# RStudio login user (no sudo)
sudo adduser rstudio
```

Open firewall port 8787: VPC Network → Firewall → Create Rule, ingress, TCP 8787, source `0.0.0.0/0`.

Inside RStudio Server (access via `http://[external-ip]:8787` as `rstudio`), set the
Posit binary mirror to avoid compiling packages from source:

```r
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
# add this line to ~/.Rprofile to make it permanent

install.packages(c("parallelDist", "fastcluster", "tidyverse", "tictoc"))
```

Upload `x_df_25km.rds` via the RStudio Files panel (or use a Cloud Storage
bucket + `gsutil cp` for larger files).

## Day-to-day workflow

1. VM Instances → **Start** the instance
2. Note the **new external IP** (it changes every restart)
3. Open `http://[external-ip]:8787` and log in as `rstudio`
4. Run the script
5. Close the browser tab, then VM Instances → **Stop**

## Cloud-generated artefacts

| File | Script | Date | Notes |
|------|--------|------|-------|
| `data/hc_fastcluster_25km.rds` | `scripts/cloud_compute_hclust_25km.R` | 2026-05-22 | Ward.D2 hclust on 116,341 cells. ~5 min total runtime. |

To regenerate any artefact: spin up the instance, upload `x_df_25km.rds` if not already there, and run the corresponding script.

## Gotchas

- **External IP changes on every restart** — grab the new one from the VM Instances list before opening RStudio.
- **System dependencies need the SSH user**, not the rstudio user (no sudo). If an R package fails with a missing `.so`, 
install the underlying library via `sudo apt install` from the browser SSH terminal.
- **Memory peaks late in long runs** — monitor with `free -h` from the RStudio Terminal during execution; usage can climb to 90+ GB even on jobs that start lightly.
- **Stopping ≠ deleting.** Stop to pause billing; only disk costs persist (pennies/day). Deleting the instance loses everything.
- **Set a billing alert** (Billing → Budgets & Alerts) at a sensible threshold so a forgotten running instance doesn't surprise you.

## Paper methods citation

> Analyses were conducted on a Google Cloud Compute Engine instance (e2-highmem-32, 32 vCPUs, 256 GB RAM) 
running R version X.X.X. Bray-Curtis dissimilarities were computed via `parallelDist::parDist` (vX.X.X) 
and hierarchical clustering via `fastcluster::hclust` (vX.X.X) with Ward.D2 linkage, to circumvent R's `.C` 
long-vector limit which prevents `analogue::distance()` from running at this resolution.
