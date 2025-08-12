library(tidyverse)
library(parallel)
library(abind)
library(pbapply)


source("NTF/scripts/fit_all_contaminants.R") # This may have to be run multiple times, depending on how stable is your system
source("NTF/scripts/estimate_all_clusterings.R")

source("NTF/scripts/load_rivm_db.R")

get_all_species <- function() {
  get(load("data/rivm_db.Rdata")) %>%
    dplyr::select(species) %>%
    unique() %>%
    unlist() %>%
    as.character()
}

get_species_names <- function(CAS) {
  get_log_dat(CAS, cens = T, centred_scaled = TRUE, geomean = T, filt_hickey = T, names = T)$species
}


get_all_tested_contaminants <- function() {
  list.files("NTF/saves/clustering_estimates/") %>%
    sapply(FUN = function(x) {
      x %>%
        str_split_1(pattern = "minVI") %>% 
        first
    }) %>%
    return()
}

get_species_to_idx_converter <- function() {
  sp <- get_all_species()
  n_sp <- sp %>% length()

  cv_dic <- c(1:n_sp, sp) %>%
    setNames(c(sp, 1:n_sp %>% as.character()))

  cv_fun <- function(labels) {
    labels %>%
      as.character() %>%
      sapply(function(label) cv_dic[[label]])
  }

  return(cv_fun)
}

species_to_idx_converter <- get_species_to_idx_converter()


get_clustering <- function(CAS) {
  fname <- paste("NTF/saves/clustering_estimates/", CAS, "minVI.rds", sep = "")

  fit1.VI <- readRDS(fname)

  data.frame(
    nm = get_species_names(CAS),
    alloc = fit1.VI %>% 
      as.vector()
  ) %>%
    mutate(CAS = CAS) %>%
    return()
}


create_species_association_matrix_for_one_c <- function(CAS, species_to_idx_converter, get_all_species) {
  n_species <- get_all_species() %>% length()

  association_matrix <- matrix(data = 0, nrow = n_species, ncol = n_species)

  clustering <- CAS %>%
    get_clustering()

  clusters <- unique(clustering$alloc)

  # Fill the diagonal elements
  for (sp in clustering$nm) {
    idx <- sp %>%
      species_to_idx_converter() %>%
      as.numeric()
    association_matrix[idx, idx] <- 1
  }

  # Build a block for the cluster
  for (cluster in clusters) {
    species_in_the_cluster <- clustering %>%
      subset(alloc == cluster) %>%
      .$nm

    if (length(species_in_the_cluster) > 1) {
      duplets <- species_in_the_cluster %>%
        species_to_idx_converter() %>%
        combn(2, simplify = F) %>%
        lapply(as.numeric)

      for (duplet in duplets) {
        association_matrix[duplet[1], duplet[2]] <- 1
        association_matrix[duplet[2], duplet[1]] <- 1
      }
    }
  }
  return(association_matrix)
}


create_full_association_tensor_noNA <- function() {
  
  dir.create("saves", recursive = TRUE, showWarnings = FALSE)
  
  get_all_tested_contaminants() %>%
    pbapply::pblapply(FUN = function(CAS) {create_species_association_matrix_for_one_c(CAS = CAS, 
                                                                            get_all_species = get_all_species, 
                                                                            species_to_idx_converter = species_to_idx_converter)}) %>%
    Reduce(function(x, y) abind::abind(x, y, along = 3), .) %>%
    saveRDS("saves/species_association_tensor_noNA.Rdata")
}
create_full_association_tensor_noNA()
