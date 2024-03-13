library(tidyr)
library(broom)
library(Ckmeans.1d.dp)
library(fields)
library(ggthemes)
library(NMF)
source("~/Documents/Github/BNP_SSD/NTF/scripts/helper_functions.R")
library(nnTensor)
library(multiway)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(ComplexHeatmap)
library("ggpubr")

fig_path <- "~/Documents/Github/BNP_SSD/NTF/figures/"

scale_rows = function(M){
  M  %>% 
    apply(MARGIN = 1, FUN = function(x) x/sum(x))
}

scale_columns = function(M){
  M  %>% 
    apply(MARGIN = 2, FUN = function(x) x/sum(x))
}


tensor_species_noNA = readRDS('~/Documents/Github/BNP_SSD/data/species_association_tensor_noNA.Rdata')
count_tests = apply(tensor_species_noNA, MARGIN = 3, FUN = function(m){rowSums(m)})
species_na = which(rowSums(count_tests) == 0) # species not tested for any contaminant
length(species_na)
#diag(count_tests) # Number of species tested for each contaminants

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

X_na = tensor_species_NA1

print(paste0("With ", dim(X_na)[2], " species and ", dim(X_na)[3], " contaminants."))
print(paste0("With ",round(as.numeric(table(is.na(X_na))['TRUE'])/prod(dim(X_na))*100,2), " % of NA in the data."))

source("~/Documents/Github/BNP_SSD/NTF/scripts/gathering_function.R")
pdf(paste0(fig_path, "species_data.pdf"), width = 10, height= 5)
repart_species <- plot_species_repart()
dev.off()
# pdf(paste0(fig_path, "contaminant_data.pdf"), width = 10, height= 5)
# plot_cont_species_tested()
# dev.off()


### Plot cross_V ###
crossV <- readRDS("~/Documents/Github/BNP_SSD/NTF/saves/crossV_data.rds")
pdf(paste0(fig_path, "cross_validation_figure.pdf"), width = 12, height= 6)
cross_validation_plot_label(crossV, 47, 5)
dev.off()

cross_validation_table <- crossV %>% as.data.frame() %>%
  group_by(rank) %>%
  summarize(min=t.test(error_mean, conf.level = 0.95)$conf.int[1], max=t.test(error_mean, conf.level = 0.95)$conf.int[2],
            error_mean = mean(error_mean))

i <- which(cross_validation_table$error_mean==min(cross_validation_table$error_mean))
error_sup <- cross_validation_table$max[i]
rank <- min(cross_validation_table[which(cross_validation_table$error_mean<=error_sup),]$rank)
print(rank)

### Select rank ###
ntf_rm <- parafac(X_na, nfac=rank, nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F)
saveRDS(ntf,file=paste0("~/Documents/Github/BNP_SSD/NTF/saves/ntf_chosen_rank_", rank,".rds"))

ntf_rm$C %>% 
  scale_columns() %>% 
  plot_kmeans_decomposition()

pdf(paste0(fig_path, "contaminants_heatmap.pdf"), width = 10, height= 16)
ntf_rm$C %>% 
  scale_columns() %>% 
  aheatmap()
dev.off()

ntf_rm$B %>% 
  scale_columns() %>% 
  plot_kmeans_decomposition()

pdf(paste0(fig_path, "species_heatmap.pdf"), width = 10, height= 16)
ntf_rm$B %>% 
  scale_columns() %>% 
  aheatmap()
dev.off()

# chosen_comp = c("X1", "X2", "X9", "X4", "X5", "X8", "X26")
chosen_comp = c("X1", "X19", "X26", "X9", "X2", "X21", "X24")
un_sep = c("X20", "X27")

pdf(paste0(fig_path, "contaminants_kmeans.pdf"), width = 6, height= 6)
ntf_rm$C %>% 
  scale_columns() %>% 
  plot_kmeans_decomposition(select = chosen_comp, unsep = un_sep, name = T)
dev.off()

pdf(paste0(fig_path, "species_kmeans.pdf"), width = 6, height= 6)
ntf_rm$B %>% 
  scale_columns() %>% 
  plot_kmeans_decomposition(select = chosen_comp, unsep = un_sep, name = T)
dev.off()

### CONTAMINANTS
comp_cont <- data.frame(ntf_rm$C) 
names(comp_cont)[which(names(comp_cont) %in% chosen_comp)] <- LETTERS[1:length(chosen_comp)]
comp_cont <- comp_cont %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  mutate(cluster_id = Ckmeans.1d.dp(value, k=2)[['cluster']]) %>%
  mutate(sp_id = seq_along(variable), cluster_id = cluster_id-1) %>%
  mutate(id = sp_id, name = sp_id %>% sapply(function(sp) idx_to_cont[[sp]]))
comp = comp_cont[which(comp_cont$variable %in% LETTERS[1:length(chosen_comp)]),]

set.seed(3)
p = tibble(Component = comp$variable %>% as.character, Activity = comp$value, Contaminant_name = comp$name, Contaminant_index = comp$id %>% as.numeric, Contaminant_comp = comp$cluster_id %>% as.numeric) %>%
  group_by(Component) %>%
  mutate(indmax = Contaminant_index[which.max(Activity)], valmax = max(Activity)) %>%
  mutate(Chem_max = Contaminant_name[indmax %>% unique],
         jitter_amount = 0.1*runif(1),
         is_in = Contaminant_comp>0) %>%
  (function(df){
    ggplot(df, aes(x = Contaminant_index, y = Activity, colour = Component, group = Component)) +
      geom_point(size = 1+df$Activity*5) +
      #geom_text(aes(x = indmax + (23 - indmax) * 0.2, y = valmax - jitter_amount, label = Chem_max)) +
      geom_text_repel(data = subset(df, is_in), aes(label = Contaminant_name)) +
      theme_few() +
      scale_colour_ptol()
  })
pdf(paste0(fig_path, "contaminants.pdf"), width = 10, height=6)
plot(p)
dev.off()

### SPECIES
comp_species <- get_species_membership(ntf_rm$B, select = chosen_comp, unsep = un_sep, name = T) %>%
  mutate(sp_id = idx_to_species_converter(sp_id)) %>%
  mutate(major = sp_id %>% sapply(function(sp) species_major[[sp]])) %>%
  mutate(Group = major %>% sapply(function(maj) major_group[[maj]])) %>%
  mutate(phylum = major %>% sapply(function(maj) major_phylum[[maj]]))
comp_species = comp_species[which(comp_species$ind %in% LETTERS[1:length(chosen_comp)]),]

p = tibble(feature = comp_species$ind %>% as.character, presence_index = comp_species$values, Group = comp_species$Group %>% as.character) %>%
  subset(presence_index==1.) %>%
  ggplot(aes(x = feature, fill = Group, group = Group)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0)) + 
  theme_few() + 
  scale_fill_ptol()+ 
  xlab('Component') + 
  ylab('Count')

p1 = tibble(feature = comp_species$ind %>% as.character, presence_index = comp_species$values, Group = comp_species$Group %>% as.character) %>%
  subset(presence_index==1.) %>%
  ggplot(aes(x = feature, fill = Group, group = Group)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0), aes(y = ..prop..)) + 
  theme_few() +
  scale_fill_ptol() + 
  xlab('Component')+ 
  ylab('Proportion')

pdf(paste0(fig_path, "species_group.pdf"), width = 10, height= 4)
ggarrange(p, p1, common.legend=T, legend = "right")
dev.off()

p = tibble(feature = comp_species$ind %>% as.character, presence_index = comp_species$values, major = comp_species$major %>% as.character %>% as.character) %>%
  subset(presence_index==1.) %>%
  ggplot(aes(x = feature, fill = major, group = major)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0)) + 
  theme_few() + 
  scale_fill_ptol() + 
  ylab('Count') + 
  xlab('Component')

p1 = tibble(feature = comp_species$ind %>% as.character, presence_index = comp_species$values, major = comp_species$major %>% as.character %>% as.character) %>%
  subset(presence_index==1.) %>%
  ggplot(aes(x = feature, fill = major, group = major)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0), aes(y = ..prop..)) + 
  theme_few() +
  scale_fill_ptol() + 
  xlab('Component')+ 
  ylab('Proportion')

pdf(paste0(fig_path, "species_major.pdf"), width = 10, height= 4)
ggarrange(p, p1, common.legend=T, legend = "right")
dev.off()

p = tibble(feature = comp_species$ind %>% as.character, presence_index = comp_species$values, phylum = comp_species$phylum %>% as.character %>% as.character) %>%
  subset(presence_index==1.) %>%
  ggplot(aes(x = feature, fill = phylum, group = phylum)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0)) + 
  theme_few() + 
  scale_fill_ptol() + 
  ylab('Count') + 
  xlab('Component')

p1 = tibble(feature = comp_species$ind %>% as.character, presence_index = comp_species$values, phylum = comp_species$phylum %>% as.character %>% as.character) %>%
  subset(presence_index==1.) %>%
  ggplot(aes(x = feature, fill = phylum, group = phylum)) + 
  geom_bar(position = position_dodge2(preserve = "single", padding = 0), aes(y = ..prop..)) + 
  theme_few() +
  scale_fill_ptol() + 
  xlab('Component')+ 
  ylab('Proportion')
pdf(paste0(fig_path, "species_phylum.pdf"), width = 10, height= 4)
ggarrange(p, p1, common.legend=T, legend = "right")
dev.off()

