library(tidyverse)
library(broom)
library(Ckmeans.1d.dp)

get_kmeans_clustering<- function(tensor_factor, select = NULL, unsep = NULL, name = FALSE){
  if(name){
    res <- (tensor_factor) %>% 
      data.frame
    names(res)[which(names(res) %in% select)] <- LETTERS[1:length(select)]
    names(res)[which(names(res) %in% unsep)] <- LETTERS[length(select)+(1:length(unsep))]
    res <- res %>%
      gather(variable, value) %>%
      group_by(variable) %>% 
      mutate(cluster_id = Ckmeans.1d.dp(value, k=2)[['cluster']])
  }
  else
    res <- (tensor_factor) %>% 
      data.frame %>%
      gather(variable, value) %>%
      group_by(variable) %>% 
      mutate(cluster_id = Ckmeans.1d.dp(value, k=2)[['cluster']])
  return(res)
}

get_species_membership <- function(tensor_factor, select = NULL, unsep = NULL, name = FALSE){
  kmeans_clustering = get_kmeans_clustering(tensor_factor, select, unsep, name)
  res = kmeans_clustering %>%
    group_by(variable) %>%
    mutate(sp_id = seq_along(variable), cluster_id = cluster_id-1) %>%
    dplyr::select(-value) 
  col <- res %>%
    spread(variable, cluster_id) %>%
    stack(select = c(-sp_id)) %>%
    group_by(ind)
  res1 <- tibble(res,col) %>%
    mutate(variable = ind) %>%
    dplyr::select(-cluster_id, -variable)
  return(res1)
}

plot_kmeans_decomposition <- function(tensor_factor, select = NULL, unsep = NULL, name = FALSE){
  
  to_plot = get_kmeans_clustering(tensor_factor, select, unsep, name)
  if(!is.null(select) && name){
    to_plot = to_plot[which(to_plot$variable %in% LETTERS[1:(length(select)+length(unsep))]),]
  }else if(!is.null(select) && !name){
    print("coucou")
    to_plot = to_plot[which(to_plot$variable %in% c(select, unsep)),]
  }
  
  to_plot %>% 
    ggplot(aes(x = value, fill = factor(cluster_id))) +
    geom_histogram() +
    facet_wrap(~variable) + 
    theme_few() + 
    scale_fill_ptol(labels = c('out', 'in'), name = 'Clustered')
}

source("~/Documents/Github/BNP_SSD/NTF/create_species_association_tensor.R")
get_column = function(column){
  get(load('rivm_db.Rdata')) %>% 
    .[,column] %>% 
    unlist() %>% 
    as.character()
}
major <- get_column('major')
species <- get_column('species')

species_major <- c(major, species) %>%
  setNames(c(species, major))

major_codes <- read.csv("major_codes.csv")
major_phylum <- c(major_codes$Major_code, major_codes$Phylum.Division) %>%
  setNames(c(major_codes$Phylum.Division, major_codes$Major_code))
major_group <- c(major_codes$Major_code, major_codes$SSD_book_group) %>%
  setNames(c(major_codes$SSD_book_group, major_codes$Major_code))


get_idx_to_species_converter <- function() {
  sp <- get_all_species()
  n_sp <- sp %>% length()
  
  cv_dic <- c(1:n_sp, sp) %>%
    setNames(c(sp, 1:n_sp %>% as.character()))
  
  species <- 1:n_sp %>%
    as.character() %>%
    sapply(function(label) cv_dic[[label]])
  species <- species[-species_rm]
  n <- length(species)
  species_name <- unname(species)
  
  cv_dic <- c(1:n, species_name) %>%
    setNames(c(species_name, 1:n %>% as.character()))
  
  cv_fun <- function(labels) {
    labels %>%
      as.character() %>%
      sapply(function(label) cv_dic[[label]])
  }
  
  return(cv_fun)
}
idx_to_species_converter <- get_idx_to_species_converter()

cont_code <- read.csv("chem_name_and_index.csv")
idx_to_cont <- c(cont_code$name, cont_code$index) %>%
  setNames(c(cont_code$index, cont_code$name))


plot_species_repart <- function() {
  species <- get_all_species()
  n_sp <- species %>% length()
  species <- species[-species_rm]
  n <- length(species)
  
  comp_species <- tibble(sp_id = species, ind = 1:n) %>%
    mutate(sp_id = idx_to_species_converter(ind)) %>%
    mutate(major = sp_id %>% sapply(function(sp) species_major[[sp]])) %>%
    mutate(Group = major %>% sapply(function(maj) major_group[[maj]])) %>%
    mutate(phylum = major %>% sapply(function(maj) major_phylum[[maj]]))
  
  p = comp_species %>%
    ggplot(aes(x = Group, fill = Group, group = Group)) + 
    geom_bar(position = 'dodge') + 
    theme_few() + 
    scale_fill_ptol() + 
    ylab('Count') + 
    xlab('Group')
  plot(p)
  
  p = comp_species %>%
    ggplot(aes(x = major, fill = major, group = major)) + 
    geom_bar(position = 'dodge') + 
    theme_few() + 
    scale_fill_ptol() + 
    ylab('Count') + 
    xlab('major')
  plot(p)
  
  p = comp_species %>%
    ggplot(aes(x = phylum, fill = phylum, group = phylum)) + 
    geom_bar(position = 'dodge') + 
    theme_few() + 
    scale_fill_ptol() + 
    ylab('Count') + 
    xlab('phylum')
  plot(p)
  return(comp_species)
}


# plot_cont_species_tested <- function() {
#   comp_cont <- tibble(sp_id = count_species, ind = 1:) %>%
#   comp_species <- tibble(sp_id = species, ind = 1:n) %>%
#     mutate(sp_id = idx_to_species_converter(ind)) %>%
#     mutate(major = sp_id %>% sapply(function(sp) species_major[[sp]])) %>%
#     mutate(Group = major %>% sapply(function(maj) major_group[[maj]])) %>%
#     mutate(phylum = major %>% sapply(function(maj) major_phylum[[maj]]))
#   
#   p = comp_species %>%
#     ggplot(aes(x = Group, fill = Group, group = Group)) + 
#     geom_bar(position = 'dodge') + 
#     theme_few() + 
#     scale_fill_ptol() + 
#     ylab('Count') + 
#     xlab('Group')
#   plot(p)
# }


