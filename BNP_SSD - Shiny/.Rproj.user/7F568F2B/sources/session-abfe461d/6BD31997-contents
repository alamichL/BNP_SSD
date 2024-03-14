library(tidyverse)
library(ggthemes)
source("~/Documents/Github/BNP_SSD/NTF/helper_functions.R")
library(multiway)
library(patchwork)
library(foreach)
library(doParallel)
registerDoParallel(cores=32)

set.seed(197)

## Loading of the data
tensor_species_noNA = readRDS('~/Documents/Github/BNP_SSD/data/species_association_tensor_noNA.Rdata')

tensor_species_NA = tensor_species_noNA
for(i in 1:dim(tensor_species_noNA)[3]){
  tensor_species_NA[diag(tensor_species_noNA[,,i])==0,,i] = NA
  tensor_species_NA[,diag(tensor_species_noNA[,,i])==0,i] = NA
}
dim(tensor_species_NA)

count_species = apply(tensor_species_noNA, MARGIN = c(1,2), FUN = function(m){sum(m)})
species_rm = which(diag(count_species) <= 12)

tensor_species_NA1 = tensor_species_NA[-species_rm,-species_rm,]
tensor_species_noNA1 = tensor_species_noNA[-species_rm,-species_rm,]

V = c()
for(c in (1:dim(tensor_species_NA1)[3])){
  v_c = sum(is.na(diag(tensor_species_NA1[,,c])))
  V = c(V,v_c)
}

X_na = tensor_species_NA1#[,,-which(V==0)]


 count_tests = apply(X_na, MARGIN = 3, FUN = function(m){rowSums(!is.na(m))})
 species_na = which(rowSums(count_tests) == 0)
 X_na <- X_na[-species_na,-species_na,]

# cross validation
cross <- cross_validation(X_na, X_na, p=0.1, rank=c(seq(1,29, by = 2), 30:50), save=T)

saveRDS(cross,file=paste0("~/Documents/Github/BNP_SSD/NTF/save/crossV_data.rds"))