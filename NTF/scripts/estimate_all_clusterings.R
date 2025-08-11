library(tidyverse)
library(parallel)
library(BNPdensity)

try_to_read <- function(path) {
  res <- try(readRDS(path))

  if (inherits(res, "try-error")) res <- get(load(path))

  return(res)
}

saved_mcmc_samples <- list.files("NTF/saves/posterior_samples_rivm")

dir.create("NTF/saves/clustering_estimates/", recursive = TRUE, showWarnings = FALSE)


estimate_optimal_clustering <- function(saved_mcmc_sample) {
  saved_mcmc_sample %>%
    try_to_read() %>%
    compute_optimal_clustering(method = "SALSO")
}

saved_mcmc_samples %>%
  mclapply(FUN = function(saved_mcmc_sample) {
    optimal_clustering_estimate <- saved_mcmc_sample %>%
      paste("NTF/saves/posterior_samples_rivm/", ., sep = "") %>%
      estimate_optimal_clustering()
    fname_save <- saved_mcmc_sample %>%
      str_split_1(pattern = "c_") %>%
      first() %>%
      paste("NTF/saves/clustering_estimates/", ., "minVI.rds", sep = "")

    saveRDS(optimal_clustering_estimate, file = fname_save)
  }, mc.cores = detectCores() - 1, mc.preschedule = F)
