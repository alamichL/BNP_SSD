homemade_nmfEstimateRank <- function(x, range, method = nmf.getOption("default.algorithm"),
                                     nrun = 30, model = NULL, ..., verbose = FALSE, stop = FALSE, with.silhouette = "both", filename_prefix = "rank") {
  homemade_nmf_cmp_save(
    x = x, range = range, method = method,
    nrun = nrun, model = model, ..., verbose = verbose, stop = stop, filename_prefix = filename_prefix
  )
  homemade_nmf_load_combine(
    x = x, range = range, method = method,
    nrun = nrun, model = model, ..., verbose = verbose, stop = stop, filename_prefix = filename_prefix
  )
}

homemade_nmf_cmp_save <- function(x, range, method = nmf.getOption("default.algorithm"),
                                  nrun = 30, model = NULL, ..., verbose = FALSE, stop = FALSE, filename_prefix = "rank") {
  k.rank <- 0
  sapply(range, function(r, ...) {
    k.rank <<- k.rank + 1L
    if (verbose) cat("Compute NMF rank=", r, " ... ")
    
    # restore RNG on exit (except after last rank)
    # => this ensures the methods use the same stochastic environment
    orng <- RNGseed()
    if (k.rank < length(range)) on.exit(RNGseed(orng), add = TRUE)
    
    res <- tryCatch({ # START_TRY
      
      res <- nmf(x, r, method, nrun = nrun, model = model, ...)
      # directly return the result if a valid NMF result
      saveRDS(object = res, file = here("rds",paste(filename_prefix, "_", r, ".rds", sep = "")))
    })
  })
}

homemade_nmf_load_combine <- function(x, range, method = nmf.getOption("default.algorithm"),
                                      nrun = 30, model = NULL, ..., verbose = FALSE, stop = FALSE, with.silhouette = "both", filename_prefix = "rank") {
  
  # fix method if passed NULL (e.g., from nmf('formula', 'numeric'))
  if (is.null(method)) {
    method <- nmf.getOption("default.algorithm")
  }
  
  # special handling of formula: get target data from the formula
  if (is(x, "formula")) {
    # dummy model to resolve formula
    dummy <- nmfModel(x, 0L, data = model)
    # retrieve target data
    V <- attr(dummy, "target")
  } else {
    V <- x
  }
  
  # remove duplicates and sort
  range <- sort(unique(range))
  
  # initiate the list of consensus matrices: start with single NA values
  c.matrices <- setNames(lapply(range, function(x) NA), as.character(range))
  fit <- setNames(lapply(range, function(x) NA), as.character(range))
  bootstrap.measures <- list()
  
  # combine function: take all the results at once and merge them into a big matrix
  comb <- function(...) {
    measures <- list(...)
    
    err <- which(sapply(measures, is.character))
    if (length(err) == length(measures)) { # all runs produced an error
      
      # build an warning using the error messages
      msg <- paste(paste("#", seq_along(range), " ", measures, sep = ""), collapse = "\n\t-")
      stop("All the runs produced an error:\n\t-", msg)
    } else if (length(err) > 0) { # some of the runs returned an error
      
      # simplify the results with no errors into a matrix
      measures.ok <- sapply(measures[-err], function(x) x)
      
      # build a NA matrix for all the results
      n <- nrow(measures.ok)
      tmp.res <- matrix(as.numeric(NA), n, length(range))
      rownames(tmp.res) <- rownames(measures.ok)
      
      # set the results that are ok
      tmp.res[, -err] <- measures.ok
      # set only the rank for the error results
      tmp.res["rank", err] <- range[err]
      # build an warning using the error messages
      msg <- paste(paste("#", err, measures[err], " ", sep = ""), collapse = "\n\t-")
      warning("NAs were produced due to errors in some of the runs:\n\t-", msg)
      
      # return full matrix
      tmp.res
    } else { # all the runs are ok
      sapply(measures, function(x) x)
    }
  }
  
  # 	measures <- foreach(r = range, .combine=comb, .multicombine=TRUE, .errorhandling='stop') %do% {
  k.rank <- 0
  measures <- sapply(range, function(r, ...) {
    k.rank <<- k.rank + 1L
    if (verbose) cat("Compute NMF rank=", r, " ... ")
    
    # restore RNG on exit (except after last rank)
    # => this ensures the methods use the same stochastic environment
    orng <- RNGseed()
    if (k.rank < length(range)) on.exit(RNGseed(orng), add = TRUE)
    
    res <- tryCatch(
      { # START_TRY
        
        # res <- nmf(x, r, method, nrun = nrun, model = model, ...)
        res <- readRDS(file = here("rds",paste(filename_prefix, "_", r, ".rds", sep = "")))
        
        # directly return the result if a valid NMF result
        if (!isNMFfit(res, recursive = FALSE)) {
          return(res)
        }
        
        # store the consensus matrix
        c.matrices[[as.character(r)]] <<- consensus(res)
        # store the fit
        fit[[as.character(r)]] <<- res
        
        # compute quality measures
        if (verbose) cat("+ measures ... ")
        measures <- summary(res, target = V, with.silhouette = with.silhouette)
        # measures <- summary(res, target = V, with.silhouette = "none")
        
        if (verbose) cat("OK\n")
        
        # return the measures
        measures
      } # END_TRY
      
      ,
      error = function(e) {
        mess <- if (is.null(e$call)) e$message else paste(e$message, " [in call to '", e$call[1], "']", sep = "")
        mess <- paste("[r=", r, "] -> ", mess, sep = "")
        if (stop) { # throw the error
          if (verbose) cat("\n")
          stop(mess, call. = FALSE)
        } # pass the error message
        if (verbose) message("ERROR")
        return(mess)
      }
    )
    
    # return the result
    res
  },
  ...,
  simplify = FALSE
  )
  
  measures <- do.call(comb, measures)
  
  # reformat the result into a data.frame
  measures <- as.data.frame(t(measures))
  
  # wrap-up result into a 'NMF.rank' S3 object
  res <- list(measures = measures, consensus = c.matrices, fit = fit)
  # if( conf.interval ) res$bootstrap.measure <- bootstrap.measures
  class(res) <- "NMF.rank"
  return(res)
}

is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

all_samples_in_table <- function(table){
  all_samples <- table$Sample_1 %>% # Create a vector with the names of all samples
    unique() %>%
    as.character() %>%
    c(table$Sample_2 %>%
        unique() %>%
        as.character()) %>%
    unique()
  return(all_samples)
}


create_table_with_sample_indices <- function(table) {
  all_samples <- table$Sample_1 %>% # Create a vector with the names of all samples
    unique() %>%
    as.character() %>%
    c(table$Sample_2 %>%
        unique() %>%
        as.character()) %>%
    unique()
  nsamples <- length(all_samples)
  
  table_with_sample_indices <- table %>%
    mutate(
      first_sample_index = as.numeric(factor(Sample_1, levels = all_samples)),
      second_sample_index = as.numeric(factor(Sample_2, levels = all_samples))
    ) # Create a index for each sample, preserving the original sample name
  return(table_with_sample_indices)
}


table_to_matrix <- function(table, val = "Dr_core") {
  table_with_sample_indices <- create_table_with_sample_indices(table)
  
  nsamples <- max(table_with_sample_indices$first_sample_index, table_with_sample_indices$second_sample_index)
  
  res <- matrix(data = NA, nrow = nsamples, ncol = nsamples)
  
  # Fill the matrix with entries from the table
  for (l in 1:nrow(table_with_sample_indices)) {
    res[
      table_with_sample_indices$first_sample_index[l],
      table_with_sample_indices$second_sample_index[l]
    ] <- table_with_sample_indices[[l, val]]
    
    res[table_with_sample_indices$second_sample_index[l], table_with_sample_indices$first_sample_index[l]] <- res[
      table_with_sample_indices$first_sample_index[l],
      table_with_sample_indices$second_sample_index[l]
    ] # Make the matrix symmetric
  }
  
  for (i in 1:nsamples) {
    res[i, i] <- 0 # The self divergence is 0
  }
  
  return(res)
  # return(table_with_sample_index %>% select(first_sample_index, second_sample_index, Dr))
}

barebone_heatmap <- function(mat) {
  # A function plotting a heatmap without anything fancy
  mat %>%
    as_tibble() %>%
    rowid_to_column(var = "sample1") %>%
    gather(sample2, val, -sample1) %>%
    mutate(sample2 = gsub("V", "", sample2) %>% as.numeric()) %>%
    ggplot(aes(x = sample1, y = sample2, fill = val)) +
    theme_bw() +
    geom_tile() +
    scale_color_viridis_c()
}

compute_NA_weight_matrix <- function(X) {
  # Returns a matrix just like X with 0 when an entry is NA and 1 otherwise
  X %>%
    apply(c(1, 2), function(x) ifelse(test = is.na(x), yes = 0, no = 1))
}

clean_NA_matrix <- function(X) {
  # replace NA by zeros
  X %>%
    apply(c(1, 2), function(x) ifelse(test = is.na(x), yes = 0, no = x))
}

#fonction qui n'impacte pas les NA
threshold_matrix <- function(mat, minthreshold = -Inf, maxthreshold = Inf) {
  mat %>%
    apply(c(1, 2), function(x) max(x, minthreshold)) %>%
    apply(c(1, 2), function(x) min(x, maxthreshold))
}

divergence_mat_to_scaled_similarity_mat <- function(divmat) {
  m <- min(1 - divmat, na.rm=T)
  scale <- 1 / (1 - m)
  divmat %>%
    apply(c(1, 2), function(x) (1 - x - m) * scale)
}

#We want to divide the previous function in two, one to convert distance to similarity and one to scale matrix
divergence_mat_to_similarity_mat <- function(divmat) {
  divmat %>%
    apply(c(1, 2), function(x) (1 - x))
}

scale_similarity_mat <- function(simmat) {
  (simmat-min(simmat,na.rm=T))/(max(simmat,na.rm=T)-min(simmat,na.rm=T))
}



table_multiple_species_to_matrix <- function(table, val = "Dr_core") {
  
  # Can't exactly call table_to_matrix multiple times for various species, as we want the sample indices to be common among species. So we use a custom function.
  
  table_with_sample_indices <- create_table_with_sample_indices(table) %>%
    mutate(species_index = Genome %>% as.factor() %>% as.numeric())
  
  nsamples <- max(table_with_sample_indices$first_sample_index, table_with_sample_indices$second_sample_index)
  
  nspecies <- length(unique(table$Genome))
  
  res <- lapply(1:nspecies, function(i) matrix(data = NA, nrow = nsamples, ncol = nsamples))
  for (l in 1:nrow(table_with_sample_indices)) {
    #print(l)
    res[[table_with_sample_indices$species_index[l]]][
      table_with_sample_indices$first_sample_index[l],
      table_with_sample_indices$second_sample_index[l]
    ] <- table_with_sample_indices[[l, val]]
    
    res[[table_with_sample_indices$species_index[l]]][table_with_sample_indices$second_sample_index[l], table_with_sample_indices$first_sample_index[l]] <- res[[table_with_sample_indices$species_index[l]]][
      table_with_sample_indices$first_sample_index[l],
      table_with_sample_indices$second_sample_index[l]
    ]
  }
  
  for (sp_index in 1:nspecies) {
    for (i in 1:nsamples) {
      res[[sp_index]][i, i] <- 0
    }
  }
  return(res %>% reduce(.f = rbind)) # Stack matrices together
}

mean_couvcore <- function(table,l){
  value=0.5*(table[l,"couv_core_Genome1"] + table[l,"couv_core_Genome2"])
  return(value[[1]])
}

product_couv_core<- function(table,l){
  value=table[l,"couv_core_Genome1"] * table[l,"couv_core_Genome2"]
  return(value[[1]])
}

nbsites <- function(table,l){
  value=table[l,"Nbsites_couv_core"]
  return(value[[1]])
}

mean_nbsites<- function(table,l){
  value=0.5*(table[l,"couv_core_Genome1"]/table[l,"Nbsites_couv_core"] + table[l,"couv_core_Genome2"]/table[l,"Nbsites_couv_core"])
  return(value[[1]])
}

geometric_mean<- function(table,l){
  value=sqrt(table[l,"couv_core_Genome1"]*table[l,"couv_core_Genome2"])/table[l,"Nbsites_couv_core"]
  return(value[[1]])
}


table_to_weight_matrix <- function(table,fun2){
  
  table <- create_table_with_sample_indices(table) %>%
    mutate(species_index = Genome %>% as.factor() %>% as.numeric())
  
  nsamples <- max(table$first_sample_index, table$second_sample_index)
  
  nspecies <- length(unique(table$Genome))
  
  W <- lapply(1:nspecies, function(i) matrix(data = 0, nrow = nsamples, ncol = nsamples))
  for (l in 1:nrow(table)) {
    #print(table_with_sample_indices$species_index[l])
    W[[table$species_index[l]]][
      table$first_sample_index[l],
      table$second_sample_index[l]]<-fun2(table,l)
    
    #Deuxième moitié matrice
    W[[table$species_index[l]]][table$second_sample_index[l], table$first_sample_index[l]] <- W[[table$species_index[l]]][
      table$first_sample_index[l],
      table$second_sample_index[l]
    ]
  }
  
  for (sp_index in 1:nspecies) {
    for (i in 1:nsamples) {
      W[[sp_index]][i, i] <- 0
    }
  }
  #return(W<-W %>% reduce(.f = rbind))
  return(W)
}

remove_diag <- function(list_matrix){
  for(i in 1:length(list_matrix)){
    list_matrix[[i]]<-t(matrix(t(list_matrix[[i]])[which(!is.na(list_matrix[[i]]))],nrow=dim(list_matrix[[i]])[2]-1,ncol=dim(list_matrix[[i]])[1]))
  }
  return(list_matrix %>% reduce(.f = rbind))
}


##############################################################################################
#Tensor functions
#############################################################################################


### Création tensor for parafac function
create_tensor_from_simmat <- function(simmat,num_samples,num_species){
  for_array <- as.data.frame(simmat)
  for_array <- split(for_array, factor(sort(rank(row.names(for_array))%%num_species)))
  new.arr <- array(as.vector(unlist(for_array)), dim = c(num_samples,num_samples,num_species))
  return(new.arr)
}

### imputation NA values with big rank
#Compute big rank ntf :
# num_species=3
# num_samples=236
#simmat=X
# rank=20
replace_na_values <- function(simmat, rank=20, num_samples, num_species){
 # cl <- makeCluster(20)
 # ce <- clusterEvalQ(cl, library(multiway))
 # clusterSetRNGStream(cl, 1)
  ntf <- parafac(simmat,nfac=rank,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=T)
 # stopCluster(cl)
  t <- replicate(num_species, diag(num_samples), simplify=F)
  for(i in 1:num_samples){
    for(j in 1:num_samples){
      for(k in 1:num_species){
        t[[k]][i,j] <-sum(ntf$A[i,]*ntf$B[j,]*ntf$C[k,])
      }
    }
  }
  t <- array(as.vector(unlist(t)), dim = c(num_samples,num_samples,num_species))
  #need to replace in original mat where NA value present
  for(i in 1:num_samples){
    for(j in 1:num_samples){
      for(k in 1:num_species){
        if(is.na(simmat[i,j,k])){
          simmat[i,j,k]<-t[i,j,k]
        }
      }
    }
  }
  return(simmat)
}

replace_na_values_binary <- function(simmat, rank=20, num_samples, num_species){
 # cl <- makeCluster(20)
 # ce <- clusterEvalQ(cl, library(multiway))
 # clusterSetRNGStream(cl, 1)
  ntf <- parafac(simmat,nfac=rank,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=T)
 # stopCluster(cl)
  t <- replicate(num_species, diag(num_samples), simplify=F)
  for(i in 1:num_samples){
    for(j in 1:num_samples){
      for(k in 1:num_species){
        t[[k]][i,j] <-sum(ntf$A[i,]*ntf$B[j,]*ntf$C[k,])
      }
    }
  }
  t <- array(as.vector(unlist(t)), dim = c(num_samples,num_samples,num_species))
  #need to replace in original mat where NA value present
  for(i in 1:num_samples){
    for(j in 1:num_samples){
      for(k in 1:num_species){
        if(is.na(simmat[i,j,k])){
          simmat[i,j,k]<-as.integer(t[i,j,k]>0.5)
        }
      }
    }
  }
  return(simmat)
}
	   
choose_rank <- function(simmat_clean, maxrank){
  criteria <- c()
  GCV <- c()
 # cl <- makeCluster(20)
 # ce <- clusterEvalQ(cl, library(multiway))
 # clusterSetRNGStream(cl, 1)
  for(i in 1:maxrank){
    print(i)
    ntf <- parafac(simmat_clean,nfac=i,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F)
    criteria[i] <- corcondia(simmat_clean,ntf)
    GCV[i] <- ntf$GCV
    #tmp <- c(i,corcondia(simmat_clean,ntf),ntf$GCV)
    #tmp
  }
  criteria <- as.data.frame(cbind(1:maxrank,criteria,GCV))
  colnames(criteria) <- c("V1","criteria","GCV")
 # stopCluster(cl)
  return(as.data.frame(criteria))
}


choose_rank_plot <- function(criteria_table){
  p <-ggplot(criteria_table,aes(x = V1, y = criteria)) +
    geom_point() +
    geom_line() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),) +
    scale_x_continuous(name="Rank",breaks=c(1:max(criteria_table$V1)),labels=c(1:max(criteria_table$V1)),limits = c(1,max(criteria_table$V1)))+
    ylab("Core consistency") +
    ggtitle("Bro and Kiers’s core consistency for each rank") +
    coord_cartesian(ylim=c(0, 100))
  return(p)
}

plot_GCV_multiway <- function(criteria_table){
  p <-ggplot(criteria_table,aes(x = V1, y = GCV)) +
    geom_point() +
    geom_line() +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),) +
    scale_x_continuous(name="Rank",breaks=c(1:max(criteria_table$V1)),labels=c(1:max(criteria_table$V1)),limits = c(1,max(criteria_table$V1)))+
    ylab("Cross validation") +
    ggtitle("GCV value of parafac function for each rank computed")
  return(p)
}

tensor_plot3D <- function(tensor,type="raw"){
  nb_spe <- dim(tensor)[3]
  nb_samp <- dim(tensor)[1]
  #Construction palette 
  ii <- cut(as.vector(tensor), breaks = seq(min(as.vector(tensor)), max(as.vector(tensor)), len = 15), 
            include.lowest = TRUE)
  colors <- colorRampPalette(c("#FFFFFF","#EAE99B", "#EC9839", "#C50000"))(14)[ii]
  #Ticks and breaks for species axis
  axz <- list(
    nticks = nb_spe, #nombre espèce
    range = c(1,nb_spe)) #nombre espèce
  #Create all coordinates
  x <- rep(rep(1:nb_samp,each=nb_samp),nb_spe) #50 = nombre échantillon
  y <- rep(rep(rep(1:nb_samp),nb_samp),nb_spe) #5 = espèce
  z <- rep(1:nb_spe,each=nb_samp*nb_samp) #2500 = échantillon*échantillon
  data <- as.data.frame(cbind(x=x,y=y,z=z,colour=colors,values=as.vector(tensor)))
  data %>%
    plot_ly(x = x, y = y, z = z, type='scatter3d',mode='markers', marker = list(color = ~colour, size = 5, colorscale=colors)) %>%
    layout(scene = list(zaxis=axz), title = paste0("Tensor visualization on ",type," data"))
}
# reconstruction tenseur avec output parafac
reconstruct_tensor <- function(parafac_result){
  t <- replicate(dim(parafac_result$C)[1], diag(dim(parafac_result$A)[1]), simplify=F)
  for(i in 1:dim(parafac_result$A)[1]){
    for(j in 1:dim(parafac_result$A)[1]){
      for(k in 1:dim(parafac_result$C)[1]){
        t[[k]][i,j] <-sum(parafac_result$A[i,]*parafac_result$B[j,]*parafac_result$C[k,])
      }  
    }
  }
  t <- array(as.vector(unlist(t)), dim = c(dim(parafac_result$A)[1],dim(parafac_result$A)[1],dim(parafac_result$C)[1]))
  return(t)
}

reconstruct_tensor_one_comp <- function(parafac_result, comp = 1){
  t <- replicate(dim(parafac_result$C)[1], diag(dim(parafac_result$A)[1]), simplify=F)
  for(i in 1:dim(parafac_result$A)[1]){
    for(j in 1:dim(parafac_result$A)[1]){
      for(k in 1:dim(parafac_result$C)[1]){
        t[[k]][i,j] <-sum(parafac_result$A[i,comp]*parafac_result$B[j,comp]*parafac_result$C[k,comp])
      }  
    }
  }
  t <- array(as.vector(unlist(t)), dim = c(dim(parafac_result$A)[1],dim(parafac_result$A)[1],dim(parafac_result$C)[1]))
  return(t)
}

remove_values_array <- function(X, index){
  X2 <- X; #Création d'une matrice temporaire
  X2[index] <- NA; #Remplacement des entrées choisies par les NAs
  return(X2)
}

distance_frobenius <- function(reconstruct, raw, index){
  return(mean((reconstruct[index] - raw[index])^2, na.rm=T))
}

cross_rank <- function(X,rank,p,save,dir,name){
  error <- c()
 # cl <- makeCluster(20)
 # ce <- clusterEvalQ(cl, library(multiway))
 # clusterSetRNGStream(cl, 1)
  for(r in rank){ #boucler plusieurs fois en prenant au hasard des valeurs différentes
    print(paste0("rang : ",r))
    #enlever une proportion p des entrées
    index <- sample(which(!is.na(X)), length(which(!is.na(X)))*p); #Choix de certaines entrées au hasard à remplacer par des NAs
    X2 <- remove_values_array(X, index)
    #calcul de la décomposition ntf avec un rang r
    #multiway::indscal(X = X2, nfac = r, type = "similarity")
    ntf <- parafac(X2,nfac=r,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F)
    if(save)  
      saveRDS(ntf,file=paste0(dir,"_",name,"_",r,".rds"))
    t <- reconstruct_tensor(ntf)
    #Calcul tenseur reconstruit
    error <- append(error, distance_frobenius(t, X, index)) #erreur moyenne de reconstruction pour rang r 
  }
 # stopCluster(cl)
  return(error)
}

cross_rank_par <- function(X,rank,p,save,dir,name){
  #error <- c()
#  cl <- makeCluster(20)
#  ce <- clusterEvalQ(cl, library(multiway))
#  clusterSetRNGStream(cl, 1)
  error <- foreach(r=rank, .combine=append, .export=c("remove_values_array","reconstruct_tensor","distance_frobenius"), .packages=c("multiway"))%dopar%{
    print(paste0("rang : ",r))
    #enlever une proportion p des entrées
    index <- sample(which(!is.na(X)), length(which(!is.na(X)))*p) #Choix de certaines entrées au hasard à remplacer par des NAs
    X2 <- remove_values_array(X, index)
    #calcul de la décomposition ntf avec un rang r
    ntf <- parafac(X2,nfac=r,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F)
    if(save)
      saveRDS(ntf,file=paste0(dir,name,r,".rds"))
    t <- reconstruct_tensor(ntf)
    #Calcul tenseur reconstruit
    tmp <- distance_frobenius(t, X, index)
    tmp #erreur moyenne de reconstruction pour rang r 
  }
  return(error)
}


zero_component <- function(X,rank,p){
  error <- c()
  nb_values <- prod(dim(X))
  #enlever une proportion p des entrées
  index <- sample(which(!is.na(X)), length(which(!is.na(X)))*p); #Choix de certaines entrées au hasard à remplacer par des NAs
  X2 <- remove_values_array(X, index)
  X_mean <- mean(X2,na.rm = TRUE)
  t <- array(rep(X_mean,nb_values), dim = c(dim(X)[1],dim(X)[1],dim(X)[3]))
  #Calcul tenseur reconstruit
  error <- distance_frobenius(t, X, index)
  return(rep(error,length(rank)))
}

oracle <- function(X,X_no_noise,rank,p){
  error <- c()
  nb_values <- prod(dim(X))
  #enlever une proportion p des entrées
  index <- sample(which(!is.na(X)), length(which(!is.na(X)))*p); #Choix de certaines entrées au hasard à remplacer par des NAs
  X2 <- remove_values_array(X, index)
  t <- X2
  t[which(is.na(X2))] <- X_no_noise[which(is.na(X2))]
  #Calcul tenseur reconstruit
  error <- distance_frobenius(t, X, index)
  return(rep(error,length(rank)))
}

cross_validation_plot <- function(cross_validation_table){
  cross_validation_table <- cross_validation_table %>% as.data.frame()
  #tracer les plots correspondants
 p1 <- cross_validation_table %>%
   group_by(rank) %>%
   summarize(min=t.test(error_mean, conf.level = 0.95)$conf.int[1], max=t.test(error_mean, conf.level = 0.95)$conf.int[2],
             error_mean = mean(error_mean), zero=zero,oracle=oracle) %>%
   ggplot(aes(x = rank)) +
   geom_line(aes(y=error_mean),alpha=0.4) +
   geom_line(aes(y=zero),alpha=0.4) +
   geom_line(aes(y=oracle),alpha=0.4) +
   geom_point(aes(y=error_mean),shape=1) +
   geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                 position=position_dodge(0.05)) +
   theme(panel.background = element_rect(fill = "white"),
         panel.grid = element_line(color = alpha("grey",0.2)),
         panel.border = element_rect(color = "grey", fill = NA, size = 1),
         legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=c(1:max(cross_validation_table$rank)),labels=c(1:max(cross_validation_table$rank)),limits = c(1,max(cross_validation_table$rank))) +
    ggtitle(paste0("Erreur moyenne par validation croisée, moyennée sur ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))

 # p2 <- cross_validation_table %>%
 #   group_by(rank) %>%
 #   summarise(error_mean = mean(error_mean),zero=zero,oracle=oracle) %>%
 #   ggplot(aes(x = rank)) +
 #   geom_line(aes(y=error_mean),alpha=0.4) +
 #   geom_line(aes(y=zero),alpha=0.4) +
 #   geom_point(aes(y=error_mean),shape=1) +
 #   theme(panel.background = element_rect(fill = "white"),
 #         panel.grid = element_line(color = alpha("grey",0.2)),
 #         panel.border = element_rect(color = "grey", fill = NA, size = 1),
 #         legend.position = "none") +
 #   ylab("Frobenius distance") +
 #   scale_x_continuous(name="Rank",breaks=c(1:max(cross_validation_table$rank)),labels=c(1:max(cross_validation_table$rank)),limits = c(1,max(cross_validation_table$rank))) +
 #   ggtitle(paste0("Erreur moyenne par validation croisée, moyennée sur ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))

  
 p3 <- cross_validation_table %>%
    ggplot(aes(x = rank, y = error_mean, color = as.factor(test)))+
    geom_line(alpha=0.4) +
    geom_point(shape=1) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),
          legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=c(1:max(cross_validation_table$rank)),labels=c(1:max(cross_validation_table$rank)),limits = c(1,max(cross_validation_table$rank))) +
    ggtitle(paste0("Erreur par validation croisée, pour ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  return(p1 + p3)
}

cross_validation_plot_label <- function(cross_validation_table, Max, by){
  cross_validation_table <- cross_validation_table %>% as.data.frame()
  cross_validation_table <- cross_validation_table[which(cross_validation_table$rank<=Max),]
  breaks = c(1,seq(5,  Max, by=by))
  #tracer les plots correspondants
  p1 <- cross_validation_table %>%
    group_by(rank) %>%
    summarize(min=t.test(error_mean, conf.level = 0.95)$conf.int[1], max=t.test(error_mean, conf.level = 0.95)$conf.int[2],
              error_mean = mean(error_mean), zero=zero,oracle=oracle) %>%
    ggplot(aes(x = rank)) +
    geom_line(aes(y=error_mean),alpha=0.4) +
    geom_line(aes(y=zero),alpha=0.4) +
    geom_line(aes(y=oracle),alpha=0.4) +
    geom_point(aes(y=error_mean),shape=1) +
    geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                  position=position_dodge(0.05)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),
          legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank", breaks=breaks, labels=breaks, limits = c(1,Max))
  
  p3 <- cross_validation_table %>%
    ggplot(aes(x = rank, y = error_mean, color = as.factor(test)))+
    geom_line(alpha=0.4) +
    geom_point(shape=1) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),
          legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=breaks, labels=breaks, limits = c(1,Max)) 
  
  return(p1 + p3)
}


cross_validation_plot_real_data <- function(cross_validation_table){
  cross_validation_table <- cross_validation_table %>% as.data.frame()
  #tracer les plots correspondants
  p1 <- cross_validation_table %>%
    group_by(rank) %>%
    summarise(sd=sd(error_mean), error_mean = mean(error_mean)) %>%
    ggplot(aes(x = rank)) +
    geom_line(aes(y=error_mean),alpha=0.4) +
    geom_point(aes(y=error_mean),shape=1) +
    geom_errorbar(aes(ymin=error_mean-sd, ymax=error_mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),
          legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=c(1,seq(5,max(cross_validation_table$rank), by=5)),labels=c(1,seq(5,max(cross_validation_table$rank), by=5)),limits = c(1,max(cross_validation_table$rank))) +
    ggtitle(paste0("Erreur moyenne par validation croisée, moyennée sur ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  # p2 <- cross_validation_table %>%
  #   group_by(rank) %>%
  #   summarise(error_mean = mean(error_mean),zero=zero,oracle=oracle) %>%
  #   ggplot(aes(x = rank)) +
  #   geom_line(aes(y=error_mean),alpha=0.4) +
  #   geom_line(aes(y=zero),alpha=0.4) +
  #   geom_point(aes(y=error_mean),shape=1) +
  #   theme(panel.background = element_rect(fill = "white"),
  #         panel.grid = element_line(color = alpha("grey",0.2)),
  #         panel.border = element_rect(color = "grey", fill = NA, size = 1),
  #         legend.position = "none") +
  #   ylab("Frobenius distance") +
  #   scale_x_continuous(name="Rank",breaks=c(1:max(cross_validation_table$rank)),labels=c(1:max(cross_validation_table$rank)),limits = c(1,max(cross_validation_table$rank))) +
  #   ggtitle(paste0("Erreur moyenne par validation croisée, moyennée sur ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  
  p3 <- cross_validation_table %>%
    ggplot(aes(x = rank, y = error_mean, color = as.factor(test)))+
    geom_line(alpha=0.4) +
    geom_point(shape=1) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
          legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=c(1,seq(5,max(cross_validation_table$rank), by=5)),labels=c(1,seq(5,max(cross_validation_table$rank), by=5)),limits = c(1,max(cross_validation_table$rank))) +
    ggtitle(paste0("Erreur par validation croisée, pour ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  return(p1 + p3)
}

log_cross_validation_plot_real_data <- function(cross_validation_table){
  cross_validation_table <- cross_validation_table %>% as.data.frame()
  #tracer les plots correspondants
  p1 <- cross_validation_table %>%
    group_by(rank) %>%
    summarise(sd=sd(log(error_mean)), error_mean = log(mean(error_mean))) %>%
    ggplot(aes(x = rank)) +
    geom_line(aes(y=error_mean),alpha=0.4) +
    geom_point(aes(y=error_mean),shape=1) +
    geom_errorbar(aes(ymin=error_mean-sd, ymax=error_mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),
          legend.position = "none") +
    ylab("log Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=c(1,seq(5,max(cross_validation_table$rank), by=5)),labels=c(1,seq(5,max(cross_validation_table$rank), by=5)),limits = c(1,max(cross_validation_table$rank))) +
    ggtitle(paste0("Erreur moyenne par validation croisée, moyennée sur ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  # p2 <- cross_validation_table %>%
  #   group_by(rank) %>%
  #   summarise(error_mean = mean(error_mean),zero=zero,oracle=oracle) %>%
  #   ggplot(aes(x = rank)) +
  #   geom_line(aes(y=error_mean),alpha=0.4) +
  #   geom_line(aes(y=zero),alpha=0.4) +
  #   geom_point(aes(y=error_mean),shape=1) +
  #   theme(panel.background = element_rect(fill = "white"),
  #         panel.grid = element_line(color = alpha("grey",0.2)),
  #         panel.border = element_rect(color = "grey", fill = NA, size = 1),
  #         legend.position = "none") +
  #   ylab("Frobenius distance") +
  #   scale_x_continuous(name="Rank",breaks=c(1:max(cross_validation_table$rank)),labels=c(1:max(cross_validation_table$rank)),limits = c(1,max(cross_validation_table$rank))) +
  #   ggtitle(paste0("Erreur moyenne par validation croisée, moyennée sur ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  
  p3 <- cross_validation_table %>%
    ggplot(aes(x = rank, y = error_mean, color = as.factor(test)))+
    geom_line(alpha=0.4) +
    geom_point(shape=1) +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = alpha("grey",0.2)),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 1),
          legend.position = "none") +
    ylab("Frobenius distance") +
    scale_x_continuous(name="Rank",breaks=c(1,seq(5,max(cross_validation_table$rank), by=5)),labels=c(1,seq(5,max(cross_validation_table$rank), by=5)),limits = c(1,max(cross_validation_table$rank))) +
    ggtitle(paste0("Erreur par validation croisée, pour ", length(1:round(1/max(cross_validation_table$p))), " essais\n",max(cross_validation_table$p)*100, " % valeurs tirées"))
  
  return(p1 + p3)
}

cross_validation <- function(X,X_no_noise, rank = seq(1,20), p = 0.1, save = FALSE, dir="~/Documents/",name="ntf_rank"){
  error_mean = c()
  # Initializes the progress bar
  # pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
  #                      max = round(1/p), # Maximum value of the progress bar
  #                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
  #                      width = 50,   # Progress bar width. Defaults to getOption("width")
  #                      char = "=")   # Character used to create the bar
  for(test in 1:30){ #boucle sur les rangs
    print(paste0("essai : ",test))
    error_mean <- append(error_mean,cross_rank(X,rank,p,save,dir,name))
    # Sets the progress bar to the current state
    #setTxtProgressBar(pb, test)
  }
  zero <- zero_component(X,rank,p)
  oracle <- oracle(X,X_no_noise,rank,p)
  table <- cbind(error_mean = error_mean, 
                 p = rep(p, length(error_mean)), 
                 rank = rep(rank, length(1:round(1/p))), 
                 test = rep(1:round(1/p),each=length(rank)),
                 zero = zero,
                 oracle = oracle)
  return(table)
}

cross_validation_par <- function(X,X_no_noise, rank = seq(1,20), p = 0.1, save = FALSE, dir="~/Documents/",name="ntf_rank"){
  #error_mean = c()
  # Initializes the progress bar
  # pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
  #                       max = round(1/p), # Maximum value of the progress bar
  #                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
  #                       width = 50,   # Progress bar width. Defaults to getOption("width")
  #                      char = "=")   # Character used to create the bar
  #
  error_mean <- foreach(test=1:round(1/p), .export='cross_rank_par', .combine=append)%dopar%{
    print(paste0("essai : ",test))#boucle sur les rangs
    tmp <-cross_rank_par(X,rank,p,save,dir,name)
    tmp
    # Sets the progress bar to the current state
    #setTxtProgressBar(pb, test)
  }
  zero <- zero_component(X,rank,p)
  oracle <- oracle(X,X_no_noise,rank,p)
  table <- cbind(error_mean = error_mean,
                 p = rep(p, length(error_mean)),
                 rank = rep(rank, length(1:round(1/p))),
                 test = rep(1:round(1/p),each=length(rank)),
                 zero = zero,
                 oracle = oracle)
  return(table)
}

raw_vs_reconstruct_points <- function(X,rank=5,species=5){
  cl <- makeCluster(20)
  ce <- clusterEvalQ(cl, library(multiway))
  clusterSetRNGStream(cl, 1)
  ntf <- parafac(X,nfac=rank,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F,parallel=T,cl=cl)
  stopCluster(cl)
  reconstruct <- reconstruct_tensor(ntf)
  for(i in 1:species){
    data = cbind(raw = as.vector(X[,,i]),reconstruct = as.vector(reconstruct[,,i]))
    g <- data %>% as.data.frame %>% filter(!is.na(raw),!is.na(reconstruct))%>%ggplot(aes(x=raw,y=reconstruct)) +
      geom_point() +
      theme(panel.background = element_rect(fill = "white"),
            panel.grid = element_line(color = alpha("grey",0.2)),
            panel.border = element_rect(color = "grey", fill = NA, size = 1),
            legend.position = "none") +
      geom_abline(a=1,col="red") +
      xlab("Raw data") +
      ylab("Reconstruct data") +
      ggtitle(paste0("Values of raw and reconstruct matrixes for species ",i))
    print(g)
  }
  
}

heatmap_each_3comp <- function(X){
  #compute ntf
  ntf <- parafac(X,nfac=5,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F)
  #get cluster assignment and value of coefficient
  A <- cbind(cluster = apply(ntf$A,1,which.max),ID = c(1:nrow(ntf$A))) %>% as.data.frame()
  A["coef"] = apply(ntf$A,1,max)
  #order by cluster and values
  A <- A %>% arrange(cluster, coef)
  # compute reconstruction for each component
  comp1 <- reconstruct_tensor_one_comp(ntf)
  comp2 <- reconstruct_tensor_one_comp(ntf,2)
  comp3 <- reconstruct_tensor_one_comp(ntf,3)
  comp4 <- reconstruct_tensor_one_comp(ntf,4)
  comp5 <- reconstruct_tensor_one_comp(ntf,5)
  for(i in 1:dim(comp1)[3]){
    #reorder matrix with preceding order found (reconstruct)
    species1 <- comp1[A$ID,A$ID,i] 
    species2 <- comp2[A$ID,A$ID,i]
    species3 <- comp3[A$ID,A$ID,i]
    species4 <- comp4[A$ID,A$ID,i]
    species5 <- comp5[A$ID,A$ID,i]
    # annotation for heatmap
    row_ha =rowAnnotation(cluster = as.factor(A$cluster), col=list(cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E")),show_annotation_name = TRUE)
    column_ha = HeatmapAnnotation(cluster = as.factor(A$cluster), col=list(cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E")), show_annotation_name = TRUE)
    # One heatmap per component
    h1 <- Heatmap(species1, name="Comp1", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
    h2 <- Heatmap(species2, name="Comp2", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
    h3 <- Heatmap(species3, name="Comp3", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
    h4 <- Heatmap(species4, name="Comp4", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
    h5 <- Heatmap(species5, name="Comp5", cluster_rows = F, cluster_columns = F, top_annotation = column_ha, right_annotation = row_ha)
    ht_list <- h1 + h2 + h3 + h4 +h5
    draw(ht_list, merge_legend = TRUE)
  }
}

create_one_group_struct_similarity_matrix = function(ngroups = 3, mydim){
  locations = 10 * 1:ngroups # 10 means well-separated
  x = rnorm(n = mydim[1], mean = 0) + sample(locations, size = mydim[1], replace = T)
  # check with density(x) %>% plot()
  sim_mat <- matrix(x, ncol = 1) %>% dist %>% as.matrix() %>% divergence_mat_to_scaled_similarity_mat
  # check with sim_mat %>% NMF::aheatmap(Rowv = NA, Colv = NA)
  # and also with library(NMF); sim_mat %>% NMF::nmf(rank = 1:(2*ngroups)) %>% plot
}
