##############################################################################################
#Tensor functions
#############################################################################################


### Create tensor for parafac function
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
  #enlever une proportion p des entrées
  index <- sample(which(!is.na(X)), length(which(!is.na(X)))*p); #Choix de certaines entrées au hasard à remplacer par des NAs
  X2 <- remove_values_array(X, index)
  for(r in rank){ #boucler plusieurs fois en prenant au hasard des valeurs différentes
    print(paste0("rang : ",r))
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
  index <- sample(which(!is.na(X)), length(which(!is.na(X)))*p) #Choix de certaines entrées au hasard à remplacer par des NAs
  X2 <- remove_values_array(X, index)
  error <- foreach(r=rank, .combine=append, .export=c("remove_values_array","reconstruct_tensor","distance_frobenius"), .packages=c("multiway"))%dopar%{
    print(paste0("rang : ",r))
    #enlever une proportion p des entrées
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

# raw_vs_reconstruct_points <- function(X,rank=5,species=5){
#   cl <- makeCluster(20)
#   ce <- clusterEvalQ(cl, library(multiway))
#   clusterSetRNGStream(cl, 1)
#   ntf <- parafac(X,nfac=rank,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F,parallel=T,cl=cl)
#   stopCluster(cl)
#   reconstruct <- reconstruct_tensor(ntf)
#   for(i in 1:species){
#     data = cbind(raw = as.vector(X[,,i]),reconstruct = as.vector(reconstruct[,,i]))
#     g <- data %>% as.data.frame %>% filter(!is.na(raw),!is.na(reconstruct))%>%ggplot(aes(x=raw,y=reconstruct)) +
#       geom_point() +
#       theme(panel.background = element_rect(fill = "white"),
#             panel.grid = element_line(color = alpha("grey",0.2)),
#             panel.border = element_rect(color = "grey", fill = NA, size = 1),
#             legend.position = "none") +
#       geom_abline(a=1,col="red") +
#       xlab("Raw data") +
#       ylab("Reconstruct data") +
#       ggtitle(paste0("Values of raw and reconstruct matrixes for species ",i))
#     print(g)
#   }
#   
# }

# heatmap_each_3comp <- function(X){
#   #compute ntf
#   ntf <- parafac(X,nfac=5,nstart=20,const=c("nonneg","nonneg","nonneg"),verbose=F)
#   #get cluster assignment and value of coefficient
#   A <- cbind(cluster = apply(ntf$A,1,which.max),ID = c(1:nrow(ntf$A))) %>% as.data.frame()
#   A["coef"] = apply(ntf$A,1,max)
#   #order by cluster and values
#   A <- A %>% arrange(cluster, coef)
#   # compute reconstruction for each component
#   comp1 <- reconstruct_tensor_one_comp(ntf)
#   comp2 <- reconstruct_tensor_one_comp(ntf,2)
#   comp3 <- reconstruct_tensor_one_comp(ntf,3)
#   comp4 <- reconstruct_tensor_one_comp(ntf,4)
#   comp5 <- reconstruct_tensor_one_comp(ntf,5)
#   for(i in 1:dim(comp1)[3]){
#     #reorder matrix with preceding order found (reconstruct)
#     species1 <- comp1[A$ID,A$ID,i] 
#     species2 <- comp2[A$ID,A$ID,i]
#     species3 <- comp3[A$ID,A$ID,i]
#     species4 <- comp4[A$ID,A$ID,i]
#     species5 <- comp5[A$ID,A$ID,i]
#     # annotation for heatmap
#     row_ha =rowAnnotation(cluster = as.factor(A$cluster), col=list(cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E")),show_annotation_name = TRUE)
#     column_ha = HeatmapAnnotation(cluster = as.factor(A$cluster), col=list(cluster = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A", "5" = "#66A61E")), show_annotation_name = TRUE)
#     # One heatmap per component
#     h1 <- Heatmap(species1, name="Comp1", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
#     h2 <- Heatmap(species2, name="Comp2", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
#     h3 <- Heatmap(species3, name="Comp3", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
#     h4 <- Heatmap(species4, name="Comp4", cluster_rows = F, cluster_columns = F, top_annotation = column_ha)
#     h5 <- Heatmap(species5, name="Comp5", cluster_rows = F, cluster_columns = F, top_annotation = column_ha, right_annotation = row_ha)
#     ht_list <- h1 + h2 + h3 + h4 +h5
#     draw(ht_list, merge_legend = TRUE)
#   }
# }
# 
# create_one_group_struct_similarity_matrix = function(ngroups = 3, mydim){
#   locations = 10 * 1:ngroups # 10 means well-separated
#   x = rnorm(n = mydim[1], mean = 0) + sample(locations, size = mydim[1], replace = T)
#   # check with density(x) %>% plot()
#   sim_mat <- matrix(x, ncol = 1) %>% dist %>% as.matrix() %>% divergence_mat_to_scaled_similarity_mat
#   # check with sim_mat %>% NMF::aheatmap(Rowv = NA, Colv = NA)
#   # and also with library(NMF); sim_mat %>% NMF::nmf(rank = 1:(2*ngroups)) %>% plot
# }
