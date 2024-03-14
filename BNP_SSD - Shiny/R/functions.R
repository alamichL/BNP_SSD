utils::globalVariables(c("."))

# functions
label_mandatory <- function(label) {
  tagList(
    label,
    span("*", class = "mandatory_star")
  )
}

inline <- function(x) {
  tags$div(style = "display:inline-block;", x)
}

hint <- function(x) HTML(paste0("<font color='grey'>", x, "</font>"))

zero_range <- function(x, tol = .Machine$double.eps^0.5) {
  if (length(x) == 1) {
    return(TRUE)
  }
  x <- range(na.omit(x)) / mean(na.omit(x))
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

center_scale <- function(dat, is.cens){
  if(!is.cens){
    res = list(data=(dat - mean(dat, na.rm = T))/sd(dat, na.rm = T), 
               mean=mean(dat, na.rm = T), sd=sd(dat, na.rm = T))
  }
  else{
    dat_uncensored <- sapply(1:length(dat[,1]), function(i){mean(as.numeric(dat[i,]), na.rm=T)})
    m <- mean(dat_uncensored)
    s <- sd(dat_uncensored)
    res = list(data=(dat - m)/s, mean=m, sd=s)
    # fit <- fitdistrplus::fitdistcens(censdata=dat,"norm")
    # res = list(data=(dat - fit$estimate[1])/fit$estimate[2],
    #            mean=fit$estimate[1], sd=fit$estimate[2])
  }
  res$centre = res$mean
  res$scale = res$sd
  return(res)
}

inv_centlog <- function(data, centre, scale, is.cs, is.log){
  if(is.cs)
    data <- scale*data + centre
  if(is.log)
    data <- exp(data)
  return(data)
}

clean_species <- function(data, is.cens, sp, conc_l, conc_u = NULL){
  if(!is.cens){
    for(spp in data[[sp]]){
      ind <- which(data[[sp]] == spp)
      if(length(ind)>1){
        data[[conc_l]][ind[1]] <- exp(mean(log(data[[conc_l]][ind])))
        data <- data[-ind[-1],]
      }
    }
  }
  else{
    for(spp in data[[sp]]){
      ind <- which(data[[sp]] == spp)
      if(length(ind)>1){
        x_r <- replace(data[[conc_u]], which(is.na(data[[conc_u]])), Inf)
        x_l <- replace(data[[conc_l]], which(is.na(data[[conc_l]])), -Inf)
        data[[conc_l]][ind[1]] <- min(x_l[ind])
        data[[conc_u]][ind[1]] <- max(x_r[ind])
        data <- data[-ind[-1],]
      }
    }
  }
  x_r <- replace(data[[conc_u]], which(is.infinite(data[[conc_u]])), NA)
  x_l <- replace(data[[conc_l]], which(is.infinite(data[[conc_l]])), NA)
  data[[conc_l]] <- x_l
  data[[conc_u]] <- x_r
  return(data)
}

prepare_data <- function(data, is.cens, is.log, is.center, conc_l,
                         conc_u = NULL, is.sp, sp_name){
  res <- list(original=data, log=!is.log, cs=!is.center)
  if(is.cens){
    if(is.sp){
      # Handle multiple values for a same species
      data <- clean_species(data, is.cens, sp_name, conc_l, conc_u) 
      res$original <- data
    }
    if(!is.log){
      data[[conc_l]] <- log(data[[conc_l]])
      data[[conc_u]] <- log(data[[conc_u]]) 
    }
    if(!is.center){
      dat <- data.frame(left=data[[conc_l]], right=data[[conc_u]])
      cs_dat <- center_scale(dat, is.cens)
      data[[conc_l]] <- cs_dat$data$left
      data[[conc_u]] <- cs_dat$data$right
      res$centre <- cs_dat$centre
      res$scale <- cs_dat$scale
    }
  }else{
    if(is.sp){
      # Handle multiple values for a same species
      data <- clean_species(data, is.cens, sp_name, conc_l)
      res$original <- data
    }
    if(!is.log){
      data[[conc_l]] <- log(data[[conc_l]])
    }
    if(!is.center){
      cs_dat <- center_scale(data[[conc_l]], is.cens)
      data[[conc_l]] <- cs_dat$data
      res$centre <- cs_dat$centre
      res$scale <- cs_dat$scale
    }
  }
  res$data <- data
  return(res)
}

pmixnorm = function(x, mus, sigmas, probs){
  # if(sum(probs)!=1) stop('check all the weights')
  if (length(unique(length(mus), length(sigmas), length(probs))) != 1){
    stop('Check your number of means, sds, unnormalized_weights')
  }
  
  mapply(function(mu_, sigma_, p_) p_*pnorm(q = x, mean = mu_, sd = sigma_), mus, sigmas, probs) %>% 
    (function(yy){ 
      if(!is.array(yy) || length(dn <- dim(yy)) < 2L) sum(yy)
      else rowSums(yy)
    } )
} 

qmixnorm = function(p, mus, sigmas, probs, interval = c(-10**6, 10**6)){
  # if(sum(probs)!=1) stop('check all the weights')
  if (length(unique(length(mus), length(sigmas), length(probs))) != 1){
    stop('Check your number of means, sds, unnormalized_weights')
  }
  
  sapply(p, FUN = function(pp){
    uniroot(f = function(x) {pmixnorm(x, mus, sigmas, probs)-pp}, 
            interval = interval, 
            extendInt = 'yes')$root
  })
}

compute_thinning_grid <- function(Nit, thinning_to = 200) {
  if (Nit <= 2 * thinning_to) { # Factor 2 to reduce the probability of having the same iterations selected twice
    it_retained <- 1:Nit
  }
  else {
    it_retained <- round(seq(1, Nit, length.out = thinning_to))
  }
  return(it_retained)
}

get_quantile <- function(fit, Q){
  thin <- compute_thinning_grid(length(fit$means)) 
  quantile <- sapply(thin, function(i){qmixnorm(Q, fit$means[[i]], fit$sigmas[[i]], fit$weights[[i]])})
  return(quantile)
}

####################
## Plot functions ##
plot_distributions <- function(fit, prep_dat, is.cens, text_size){
  if(!is.cens){
    m <- ncol(fit$qx)
    nbins <- length(hist(fit$data[,1], plot = FALSE)$breaks) - 1
    # nbins <- length(hist(as.data.frame(prep_dat$original)[,2], plot = FALSE)$breaks) - 1
    # gp <- ggplot2::ggplot(data.frame(xx = inv_centlog(fit$xx, prep_dat$centre, prep_dat$scale, prep_dat$cs, prep_dat$log), 
    #                                  infCI = fit$qx[, 2], supCI = fit$qx[, m], y =fit$qx[, 1]), ggplot2::aes_string(x = "xx")) +
    gp <- ggplot2::ggplot(data.frame(xx = fit$xx, infCI = fit$qx[, 2], supCI = fit$qx[, m], y =fit$qx[, 1]), 
                          ggplot2::aes_string(x = "xx")) +
      ggplot2::theme_classic() +
      ggplot2::geom_histogram(
        # data = data.frame(x = as.data.frame(prep_dat$original)[,2]), ggplot2::aes_string(x = "x", y = "..density.."),
        data = data.frame(x = fit$data[,1]), ggplot2::aes_string(x = "x", y = "..density.."),
        fill = grDevices::grey(0.9),
        colour = "black",
        bins = nbins
      ) +
      ggplot2::geom_line(ggplot2::aes_string(y = "y"), size = 1.) +
      ggplot2::geom_line(ggplot2::aes_string(y = "infCI"), colour = "blue", linetype = "dotted") +
      ggplot2::geom_line(ggplot2::aes_string(y = "supCI"), colour = "blue", linetype = "dotted") +
      ggplot2::xlab("Data") +
      ggplot2::ylab("Density")
  }else{
    # m <- ncol(fit$qx)
    # gp <- ggplot2::ggplot(data.frame(xx = inv_centlog(fit$xx, prep_dat$centre, prep_dat$scale, prep_dat$cs, prep_dat$log), 
    #                                  infCI = fit$qx[, 2], supCI = fit$qx[, m], y =fit$qx[, 1]), ggplot2::aes_string(x = "xx")) +
    #   ggplot2::theme_classic() +
    #   ggplot2::geom_line(ggplot2::aes_string(y = "y"), size = 1.) +
    #   ggplot2::geom_line(ggplot2::aes_string(y = "infCI"), colour = "blue", linetype = "dotted") +
    #   ggplot2::geom_line(ggplot2::aes_string(y = "supCI"), colour = "blue", linetype = "dotted") +
    #   ggplot2::xlab("Data") +
    #   ggplot2::ylab("Density")
    gp <- plot(fit)
  }
  gp <- gp +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12)
    )
  gp
}

boxplot_cpo <- function(x){
  y <- BNPdensity::cpo(x)
  bp <- ggplot2::ggplot() + ggplot2::geom_boxplot(ggplot2::aes(x="CPO", y=y), fill="gray", color="black") +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                                              axis.title.y=ggplot2::element_blank(),
                                              axis.text = ggplot2::element_text(size = 12),
                                              legend.text = ggplot2::element_text(size = 12))
  bp
}

plot_CDF_percentil <- function(fit, q, Q){
  gp <- BNPdensity:::plotCDF_censored(fit) +
    ggplot2::geom_vline(xintercept=mean(q), linetype="dotted")+
    ggplot2::geom_hline(yintercept=Q, linetype="dotted")
  gp
}
