library(tidyverse)
library(parallel)
library(abind)
library(rTensor)

source("~/Documents/Github/BNP_SSD/NTF/scripts/load_rivm_db.R")

get_all_species <- function() {
  get(load("rivm_db.Rdata")) %>%
    dplyr::select(species) %>%
    unique() %>%
    unlist() %>%
    as.character()
}

get_species_names <- function(CAS) {
  get_log_dat(CAS, cens = T, centred_scaled = TRUE, geomean = T, filt_hickey = T, names = T)$species
}


get_all_tested_contaminants <- function() {
  # list.files('minVIs/') %>%
  #   c(list.files('minVIs_new/')) %>%
  list.files("minVIs_new//") %>%
    (function(v) v[!grepl("nc", v)]) %>%
    unique() %>%
    gsub("c.Rdata", "", .)
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
  fname <- paste("minVIs_new/", CAS, "c.Rdata", sep = "")

  fit1.VI <- readRDS(fname)

  nclust <- length(unique(fit1.VI$cl))

  data.frame(
    nm = get_species_names(CAS),
    alloc = fit1.VI$cl
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
  get_all_tested_contaminants() %>%
    pbapply::pblapply(FUN = function(CAS) {create_species_association_matrix_for_one_c(CAS = CAS, 
                                                                            get_all_species = get_all_species, 
                                                                            species_to_idx_converter = species_to_idx_converter)}) %>%
    Reduce(function(x, y) abind(x, y, along = 3), .) %>%
    saveRDS("species_association_tensor_noNA.Rdata")
}
