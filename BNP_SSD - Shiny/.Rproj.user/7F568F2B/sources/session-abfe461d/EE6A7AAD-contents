library(compiler)
library(Rcpp)
library(survival)

#Some useful functions

save_the_plot_as_pdf = function(plot_, fname, ...){
  print(plot_)
  pdf(file = fname, ...)
  print(plot_)
  dev.off()
}


get_plotting_range = function(logdata){
  if (is.censored(logdata)){
    logddat = rowMeans(logdata, na.rm = T)
  }
  else{
    logddat = ddat
  }

  mindat = min(logddat)
  maxdat = max(logddat)
  span = maxdat-mindat

  seq(mindat-0.2*span, maxdat+0.2*span, length.out = 100) %>%
    return
}

centre_and_scale = function(ddat){
  if(is.censored(ddat)) {
    c_sc = ddat %>%
      (function(x) 0.5 * (x$left + x$right)) %>%
      (function(x) c('mean' = mean(x, na.rm = T), 'sd' = sd(x, na.rm = T)))
    res = list(data_ = (ddat - c_sc[['mean']])/c_sc[['sd']], mean = c_sc[['mean']], sd = c_sc[['sd']])
  }
  else {
    res = list(data_ = (ddat - mean(ddat, na.rm = T))/sd(ddat, na.rm = T), mean = mean(ddat, na.rm = T), sd = sd(ddat, na.rm = T))
  }
  res$centre = res$mean
  res$scale = res$sd
  return(res)
}

is.censored = function(ddat){
  if(is.null(ncol(ddat))) FALSE
  else if(ncol(ddat)==1) FALSE
  else if(ncol(ddat)==2) TRUE
  else stop('Wrong type/dim of data input')
}

make_grid_from_data_vec_no_na = function(data_vec_no_na, length_ = 100, Delta = 2){
  seq(min(data_vec_no_na)-Delta,
      max(data_vec_no_na)+Delta,
      length.out = length_) %>%
    c(data_vec_no_na) %>%
    unique() %>%
    sort()
}

make_grid_from_data = function(data, length_ = 100, Delta = 2){
  if(data %>% is.censored){
    xs = data %>% rowMeans(na.rm = T) %>%
      make_grid_from_data_vec_no_na(length_ = length_, Delta = Delta)
  }
  else{
    xs = data %>%
      na.omit() %>%
      as.numeric() %>%
      make_grid_from_data_vec_no_na(length_ = length_, Delta = Delta)
  }

}


source('scripts/mixnorm.R')

dscaled_student = function(x, df, ncp, scale, log = FALSE) {
  1 / scale * dt(
    x = (x - ncp) / scale, df = df, ncp = 0, log = log
  )
}

pscaled_student = function(x, df, ncp, scale, log = FALSE) {
  pt(
    q = (x - ncp) / scale, df = df, ncp = 0, log.p = log
  )
}


########Get quantiles

source('scripts/get_median_quantile_BNP_density.R')

get_quantiles = function(fit, ps = seq(0.01,0.99, length.out = 100), par_ = T){
  fit$data_ = centre_and_scale(fit$data_) %>% .$data_
  get_quantiles_noscale(fit, ps, par_) * fit$scale + fit$centre
}

get_quantiles_noscale = function(fit, ps = seq(0.01,0.99, length.out = 100), par_ = T) {
  #Do not use quantiles that are too small or too big or you risk a convergence problem with the root finding algorithms
  if (fit$qmethod == 'normal_ssd' |
      fit$qmethod == 'logistic_ssd')
    fit %>% quantile(probs = ps) %>% .$quantiles %>% unlist()
  else if (fit$qmethod == 'normal_mixture_ssd') {
    search_interval = c(min(fit$data_) - 5, max(fit$data_) + 5)
    ps %>%
      lapply(function(x) {
        function(y) {
          pmixnorm(
            x = y, mus = fit$mu, sigmas = fit$sigma, probs = fit$lambda
          ) - x
        }
      }) %>%
      lapply(function(x)
        uniroot(x, interval = search_interval, extendInt = "yes")$root) %>%
      unlist
  }
  else if (fit$qmethod == 'normal_kernel_mixture_ssd') {
    search_interval = c(min(fit$data_) - 5, max(fit$data_) + 5)
    ps %>%
      lapply(function(x) {
        function(y) {
          fit$pmixkernel(x = y) - x
        }
      }) %>%
      lapply(function(x)
        uniroot(x, interval = search_interval, extendInt = "yes")$root) %>%
      unlist
  }

  else if (fit$qmethod %>% grepl('DP_mixture_ssd',.)) {
    search_interval = c(min(fit$data_) - 10 ^ 6, max(fit$data_) + 10 ^ 6)


    select_pars = function(parname = 'mu') {
      if (parname %in% c('m1', 'k0', 'psi1', 'ncluster', 'alpha')) {
        fit$save.state$thetasave[1,] %>% names() %>% grepl(pattern = parname, x = .)
      }
      #Acceptable parnames: mu, sigma, m1, k0, psi1, ncluster, alpha
      else if (parname %in% c('mu', 'sigma')) {
        fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = parname, x = .)
      }
      else
        stop('no param with this name')
    }

    select_sigmas = function() {
      fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = 'sigma', x = .)
    }

    ppred = function(iter) {
      function(y) {
        1 / (fit$prior$alpha + fit$nrec) * (
          fit$prior$alpha * pscaled_student(
            x = y,
            df = 2 * fit$prior$nu1,
            ncp = fit$save.state$thetasave[iter,select_pars('m1')],
            scale = sqrt(
              (1 + 1 / fit$save.state$thetasave[iter,select_pars('k0')]) /
                (fit$save.state$thetasave[iter,select_pars('psi1')] *
                   fit$prior$nu1)
            )
          ) +
            pmixnorm(
              y ,
              mus = fit$save.state$randsave[iter,select_pars('mu')][1:fit$nrec],
              sigmas = fit$save.state$randsave[iter,select_pars('sigma')][1:fit$nrec] %>% sqrt,##Yeah, well, seems like its a sigma**2 in the end
              probs = rep(1, fit$nrec)
            )
        )
      }
    }

    cdF_list = seq_len(fit$mcmc$nsave) %>%
      lapply(ppred)

    ps %>%
      lapply(function(ps_i) {
        lapply(cdF_list, function(cdF_j) {
          function(y) {
            cdF_j(y) - ps_i
          }
        }) %>%
          mclapply(function(cdf_j)
          {uniroot(cdf_j, interval = search_interval, extendInt = "yes")$root},
          mc.cores = ifelse(test = par_, yes = detectCores(), no = 1)) %>%
          unlist %>%
          data.frame(ps = ps_i, val = .)
      }) %>%
      Reduce(rbind, .) %>%
      dplyr::group_by(ps) %>%
      dplyr::summarise(qs = median(val)) %>%
      .$qs
  }
  # else if (fit$qmethod %>% grepl('BNP_mixture_ssd',.)) {
  else if (fit$qmethod %>% grepl('BNP',.)) {
    search_interval = c(min(fit$data_) - 10 ^ 6, max(fit$data_) + 10 ^ 6)

    #
    #     ppred = function(iter) {
    #       function(y) {
    #         pmixnorm(y ,
    #                  mus = fit$means[[iter]],
    #                  sigmas = fit$sigmas[[iter]],
    #                  probs = fit$weights[[iter]])
    #       }
    #     }
    #
    #     cdF_list = seq_along(fit$means) %>%
    #       lapply(ppred)
    #
    #
    #     ps %>%
    #       lapply(function(ps_i) {
    # lapply(cdF_list, function(cdF_j) {
    #   function(y) {
    #     cdF_j(y) - ps_i
    #   }
    # }) %>%
    #   mclapply(function(cdf_j)
    #     {uniroot(cdf_j, interval = search_interval, extendInt = "yes")$root},
    #     mc.cores = ifelse(test = par_, yes = detectCores(), no = 1)) %>%
    #   unlist %>%
    #   data.frame(ps = ps_i, val = .)
    #       }) %>%
    #       Reduce(rbind, .) %>%
    #       dplyr::group_by(ps) %>%
    #       dplyr::summarise(qs = median(val)) %>%
    #       .$qs

    if (fit$qmethod %>% grepl('semi',.)){
      ps %>%
        lapply(function(ps_i) {
          get_median_quantile_semi_BNPdensity(fit, q = ps_i, search_interval = search_interval) %>%
            data.frame(ps = ps_i, val = .)
        }) %>%
        Reduce(rbind, .) %>%
        dplyr::group_by(ps) %>% #Why is this added ?
        dplyr::summarise(qs = median(val)) %>%
        .$qs
    }

    else{
      ps %>%
        lapply(function(ps_i) {
          get_median_quantile_full_BNPdensity(fit, q = ps_i, search_interval = search_interval) %>%
            data.frame(ps = ps_i, val = .)
        }) %>%
        Reduce(rbind, .) %>%
        dplyr::group_by(ps) %>% #Same here, why is this added ?
        dplyr::summarise(qs = median(val)) %>%
        .$qs
    }
  }
}

source('scripts/etx.R')

get_quantiles_and_CI = function(fit, p = 0.05, par_ = T, nbootsamples = 10**4){
  fit$data_ = centre_and_scale(fit$data_) %>% .$data_
  get_quantiles_and_CI_noscale(fit, p, par_, nbootsamples) * fit$scale + fit$centre
}


get_quantiles_and_CI_fitdistcens_object = function(fit, p = 0.05, nbootsamples = 10**4){
  fit %>%
    bootdistcens(niter = nbootsamples) %>%
    quantile(probs = p)  %>%
    (function(xx) c(xx$quantiles$`p=0.05`, xx$quantCI$`p=0.05`)) %>%
    setNames(c('HC','C_inf','C_sup')) %>%
    t %>%
    data.frame %>% mutate(p = p) %>%
    return()
}

get_quantiles_and_CI_noscale = function(fit, p = 0.05, par_ = T, nbootsamples = 10**4) {
  #Not for censored data for the moment
  #Do not use quantiles that are too small or too big or you risk a convergence problem with the root finding algorithms
  if (fit$qmethod == 'logistic_ssd'){
    if(is.censored(fit$data_)){
      get_quantiles_and_CI_fitdistcens_object(fit = fit, p = p, nbootsamples = nbootsamples)
    }
    else{
      fit %>%
        bootdist(niter = nbootsamples) %>%
        quantile(probs = p)  %>%
        (function(xx) c(xx$quantiles$`p=0.05`, xx$quantCI$`p=0.05`)) %>%
        setNames(c('HC','C_inf','C_sup')) %>%
        t %>%
        data.frame %>% mutate(p = p) %>%
        return()
    }}
  # fit %>% quantile(probs = p) %>% .$quantiles %>% unlist() %>% t %>% data.frame() %>% mutate(p = p) %>% return() #finish with bootstrap


  else if (fit$qmethod == 'normal_ssd'){
    if(is.censored(fit$data_)){
      get_quantiles_and_CI_fitdistcens_object(fit = fit, p = p, nbootsamples = nbootsamples)
    }
    else fit$data_ %>%
      etx_on_log_fast(hcx = p) %>%
      unlist() %>%
      t %>%
      data.frame() %>%
      mutate(p = p)
  }

  else if (fit$qmethod == 'normal_mixture_ssd') {
    search_interval = c(min(fit$data_, na.rm = T) - 5, max(fit$data_, na.rm = T) + 5)
    ps %>%
      lapply(function(x) {
        function(y) {
          pmixnorm(
            x = y, mus = fit$mu, sigmas = fit$sigma, probs = fit$lambda
          ) - x
        }
      }) %>%
      lapply(function(x)
        uniroot(x, interval = search_interval, extendInt = "yes")$root) %>%
      unlist #finish with bootstrap ?

  }
  else if (fit$qmethod == 'normal_kernel_mixture_ssd') {
    search_interval = c(min(fit$data_, na.rm = T) - 5, max(fit$data_, na.rm = T) + 5)
    ndat = length(fit$data_)

    make_fit = function(data_){
      #The output function is not vectorise
      hn = 1.06 * sd(data_) * ndat^(-1/5)
      return(function(x) 1/ndat*sum(pnorm(q = (x-data_)/hn, mean = 0, sd = 1)))
    }

    pnorm_minus_p = function(x){
      make_fit(fit$data_)(x)-p
    }

    q = uniroot(pnorm_minus_p, interval = search_interval, extendInt = "yes")$root

    bootstrap_ci = sample(fit$data_, size = nbootsamples*ndat, replace = T) %>%
      split(.,ceiling(seq_along(.)/ndat)) %>%
      lapply(function(x) make_fit(x)) %>%
      lapply(function(fit) {
        function(x) {
          fit(x) - p
        }
      }) %>%
      lapply(function(xfun)
        uniroot(xfun, interval = search_interval, extendInt = "yes")$root) %>%
      unlist %>%
      quantile(probs = c(0.025, 0.975) )

    return(data.frame(HC = q) %>%
             cbind(bootstrap_ci %>%
                     setNames(c('C_inf','C_sup')) %>%
                     t) %>% mutate(p = p)
    )
  }

  else if (fit$qmethod %>% grepl('DP_mixture_ssd',.)) {
    search_interval = c(min(fit$data_, na.rm = T) - 10 ^ 6, max(fit$data_, na.rm = T) + 10 ^ 6)


    select_pars = function(parname = 'mu') {
      if (parname %in% c('m1', 'k0', 'psi1', 'ncluster', 'alpha')) {
        fit$save.state$thetasave[1,] %>% names() %>% grepl(pattern = parname, x = .)
      }
      #Acceptable parnames: mu, sigma, m1, k0, psi1, ncluster, alpha
      else if (parname %in% c('mu', 'sigma')) {
        fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = parname, x = .)
      }
      else
        stop('no param with this name')
    }

    select_sigmas = function() {
      fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = 'sigma', x = .)
    }

    ppred = function(iter) {
      function(y) {
        1 / (fit$prior$alpha + fit$nrec) * (
          fit$prior$alpha * pscaled_student(
            x = y,
            df = 2 * fit$prior$nu1,
            ncp = fit$save.state$thetasave[iter,select_pars('m1')],
            scale = sqrt(
              (1 + 1 / fit$save.state$thetasave[iter,select_pars('k0')]) /
                (fit$save.state$thetasave[iter,select_pars('psi1')] *
                   fit$prior$nu1)
            )
          ) +
            pmixnorm(
              y ,
              mus = fit$save.state$randsave[iter,select_pars('mu')][1:fit$nrec],
              sigmas = fit$save.state$randsave[iter,select_pars('sigma')][1:fit$nrec] %>% sqrt,##Yeah, well, seems like its a sigma**2 in the end
              probs = rep(1, fit$nrec)
            )
        )
      }
    }

    cdF_list = from = 1:fit$mcmc$nsave %>%
      lapply(ppred)


    p %>%
      lapply(function(ps_i) {
        lapply(cdF_list, function(cdF_j) {
          function(y) {
            cdF_j(y) - ps_i
          }
        }) %>%
          mclapply(function(cdf_j)
          {uniroot(cdf_j, interval = search_interval, extendInt = "yes")$root},
          mc.cores = ifelse(test = par_, yes = detectCores(), no = 1)) %>%
          unlist %>%
          data.frame(ps = ps_i, val = .)
      }) %>%
      Reduce(rbind, .) %>%
      dplyr::group_by(ps) %>%
      do(q=data.frame(quantile(.$val, probs = c(0.5,0.025,0.975)) %>% t)) %>%
      (function(x){
        data.frame(p = x$ps) %>%
          cbind(x$q %>% Reduce(rbind,.) %>% setNames(nm = c('HC','C_inf','C_sup')))
      })
  }
  # else if (fit$qmethod %>% grepl('BNP_mixture_ssd',.)) {
  else if (fit$qmethod %>% grepl('BNP',.)) {
    search_interval = c(min(fit$data_, na.rm = T) - 10 ^ 6, max(fit$data_, na.rm = T) + 10 ^ 6)


    if (fit$qmethod %>% grepl('semi',.)){
      p %>%
        lapply(function(ps_i) {
          get_median_quantile_and_CI_semi_BNP_density(fit, q = ps_i, search_interval = search_interval) %>%
            t %>%
            data.frame() %>%
            setNames(nm = c('HC','C_inf','C_sup')) %>%
            mutate(p = ps_i)
        }) %>%
        Reduce(rbind, .) #%>%
      # dplyr::group_by(ps) %>% #Why is this added ?
      # dplyr::summarise(qs = median(val)) %>%
      # .$qs
    }

    else{
      p %>%
        lapply(function(ps_i) {
          get_median_quantile_and_CI_full_BNP_density(fit, q = ps_i, search_interval = search_interval) %>%
            t %>%
            data.frame() %>%
            setNames(nm = c('HC','C_inf','C_sup')) %>%
            mutate(p = ps_i)
        }) %>%
        Reduce(rbind, .)
    }
  }
}

#Test quantile functions
#
# qnorm(0.05, mean = 40, sd = 20)
# fit = fit_normal_ssd_on_log_data(rnorm(10**5, mean = 40, sd = 20))
# fit$data_ %>% quantile(0.05)
# fit %>% get_quantiles(0.05)
# fit %>% get_quantiles_and_CI()
#
# fit = fit_logis_ssd_on_log_data(rnorm(10**5, mean = 40, sd = 20))
# fit$data_ %>% quantile(0.05)
# fit %>% get_quantiles(0.05)
# fit %>% get_quantiles_and_CI()
#
#
# fit = fit_kernel_mixture_ssd_on_log_data(rnorm(10**4))
# fit$data_ %>% quantile(0.05)
# fit %>% get_quantiles(0.05)
# fit %>% get_quantiles_and_CI()
#
# fit = fit_DP_mixture_ssd_on_log_data(rnorm(10**4))
# fit$data_ %>% quantile(0.05)
# fit %>% get_quantiles(0.05)
# fit %>% get_quantiles_and_CI()

# fit = fit_DP_mixture_ssd_on_log_data(rnorm(10**2))
# fit$data_ %>% quantile(0.05)
# fit %>% get_quantiles(0.05)
# fit %>% get_quantiles_and_CI()

# fit = fit_BNP_mixture_ssd_on_log_data(rnorm(10**2), Nit = 200)
# fit$data_ %>% quantile(0.05)
# fit %>% get_quantiles(0.05)
# fit %>% get_quantiles_and_CI(p = 0.05)


#####Get density


source('scripts/get_dens_DP_package.R')
source('scripts/get_dens_BNPdensity.R')
source('scripts/get_median_dens_BNPdensity.R')

get_dens = function(fit, xs = seq(-5,5, length.out = 100), median = F){
  fit_ = fit
  scaled_data = centre_and_scale(fit_$data_)
  fit_$data_ = scaled_data$data_
  get_dens_no_scale(fit = fit_, xs = (xs - scaled_data$centre)/scaled_data$scale  , median = median)/scaled_data$scale
}

get_dens_no_scale = cmpfun(function(fit, xs = seq(-5,5, length.out = 100), median = F) {
  if (fit$qmethod == 'normal_ssd')
    dnorm(x = xs, mean = fit$estimate['mean'], sd = fit$estimate['sd'])
  else  if (fit$qmethod == 'logistic_ssd')
    dlogis(
      x = xs, location = fit$estimate['location'], scale = fit$estimate['scale']
    )
  else if (fit$qmethod == 'normal_mixture_ssd') {
    dmixnorm(
      x = xs, mus = fit$mu, sigmas = fit$sigma, probs = fit$lambda
    )
  }
  else if (fit$qmethod == 'normal_kernel_mixture_ssd') {
    fit$dmixkernel(x = xs)
  }

  else if ((fit$qmethod %>% grepl('DP_mixture',.))&!(fit$qmethod %>% grepl('BNP',.))) {

    if(median) get_median_dens_DP_package(fit = fit, xs = xs)

    else get_dens_DP_package(fit = fit, xs = xs)

  }

  else if (fit$qmethod %>% grepl('BNP',.)) {

    if (fit$qmethod %>% grepl('semi',.)){
      if(median) get_median_dens_semi_BNPdensity(fit, xs)
      else get_dens_semi_BNPdensity(fit, xs)
    }

    else {
      if(median) get_median_dens_full_BNPdensity(fit, xs)
      else get_dens_full_BNPdensity(fit, xs)
    }
  }
})


#Test density function
#
# fit = fit_normal_ssd_on_log_data(rnorm(10**4, mean = 40, sd = 6))
#
# hist(fit$data_, freq = F)
# seq(20,60,length.out = 100) %>%
#   (function(x){
#     lines(x, get_dens(fit, x), col = 'red')
#   })
# curve(dnorm(x, mean = 40, sd = 6), col = 'blue', add = T)

# fit = fit_logis_ssd_on_log_data(rnorm(10**4, mean = 40, sd = 6))
#
# hist(fit$data_, freq = F)
# seq(20,60,length.out = 100) %>%
#   (function(x){
#     lines(x, get_dens(fit, x), col = 'red')
#   })
#
# fit = fit_kernel_mixture_ssd_on_log_data(rnorm(10**4, mean = 40, sd = 6))
#
# hist(fit$data_, freq = F)
# seq(20,60,length.out = 100) %>%
#   (function(x){
#     lines(x, get_dens(fit, x), col = 'red')
#   })
#
# fit = fit_DP_mixture_ssd_on_log_data(rnorm(10**2, mean = 40, sd = 6))
#
# hist(fit$data_, freq = F)
# seq(20,60,length.out = 100) %>%
#   (function(x){
#     lines(x, get_dens(fit, x), col = 'red')
#   })
#
# fit = fit_BNP_mixture_ssd_on_log_data(rnorm(10**2, mean = 40, sd = 6), Nit = 200)
#
# hist(fit$data_, freq = F)
# seq(20,60,length.out = 100) %>%
#   (function(x){
#     lines(x, get_dens(fit, x), col = 'red')
#   })

#####Get CDF


source('scripts/get_CDF_BNPdensity.R')

get_cdf = function(fit, xs = seq(-5,5, length.out = 100)){
  fit_ = fit
  scaled_data = centre_and_scale(fit_$data_)
  fit_$data_ = scaled_data$data_
  get_cdf_no_scale(fit = fit_, xs = (xs - scaled_data$centre)/scaled_data$scale)
}

get_cdf_no_scale = function(fit, xs = seq(-5,5, length.out = 100)) {
  if (fit$qmethod == 'normal_ssd')
    pnorm(q = xs, mean = fit$estimate['mean'], sd = fit$estimate['sd'])
  else  if (fit$qmethod == 'logistic_ssd')
    plogis(
      q = xs, location = fit$estimate['location'], scale = fit$estimate['scale']
    )
  else if (fit$qmethod == 'normal_mixture_ssd') {
    pmixnorm(
      x = xs, mus = fit$mu, sigmas = fit$sigma, probs = fit$lambda
    )
  }
  else if (fit$qmethod == 'normal_kernel_mixture_ssd') {
    fit$pmixkernel(x = xs)
  }

  else if (fit$qmethod %>% grepl('DP_mixture_ssd',.)) {
    search_interval = c(min(fit$data_) - 10 ^ 6, max(fit$data_) + 10 ^ 6)


    select_pars = function(parname = 'mu') {
      if (parname %in% c('m1', 'k0', 'psi1', 'ncluster', 'alpha')) {
        fit$save.state$thetasave[1,] %>% names() %>% grepl(pattern = parname, x = .)
      }
      #Acceptable parnames: mu, sigma, m1, k0, psi1, ncluster, alpha
      else if (parname %in% c('mu', 'sigma')) {
        fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = parname, x = .)
      }
      else
        stop('no param with this name')
    }

    select_sigmas = function() {
      fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = 'sigma', x = .)
    }

    ppred = function(iter) {
      function(y) {
        1 / (fit$prior$alpha + fit$nrec) * (
          fit$prior$alpha * pscaled_student(
            x = y,
            df = 2 * fit$prior$nu1,
            ncp = fit$save.state$thetasave[iter,select_pars('m1')],
            scale = sqrt(
              (1 + 1 / fit$save.state$thetasave[iter,select_pars('k0')]) /
                (fit$save.state$thetasave[iter,select_pars('psi1')] *
                   fit$prior$nu1)
            )
          ) +
            pmixnorm(
              y ,
              mus = fit$save.state$randsave[iter,select_pars('mu')][1:fit$nrec],
              sigmas = fit$save.state$randsave[iter,select_pars('sigma')][1:fit$nrec] %>% sqrt,##Yeah, well, seems like its a sigma**2 in the end
              probs = rep(1, fit$nrec)
            )
        )
      }
    }

    set.seed(0)
    cdF_list = 1:fit$mcmc$nsave %>% sample(size = min(200, fit$mcmc$nsave)) %>% lapply(ppred)


    cdF_list %>%
      mclapply(function(cdf_j)
        cdf_j(xs), mc.cores = detectCores()) %>%
      Reduce(cbind,.) %>%
      rowMeans

  }

  else if (fit$qmethod %>% grepl('BNP',.)) {

    if (fit$qmethod %>% grepl('semi',.)) get_CDF_semi_BNPdensity(fit, xs)

    else get_CDF_full_BNPdensity(fit, xs)

  }
}


comparison_CDFplot = function(fit_list, xs = NULL, par_ = T, n_splitted_plot = 4) {

  if(is.null(xs)){
    xs = fit_list[[1]]$data_ %>%
      make_grid_from_data(length_ = 100, Delta = 1)
  }

  if(is.censored(fit_list[[1]]$data_)){
    survdata <- Surv(time = fit_list[[1]]$data_$left, time2 = fit_list[[1]]$data_$right,
                     type = "interval2")
  }

  else{survdata <- Surv(time = fit_list[[1]]$data_, time2 = fit_list[[1]]$data_,
                        type = "interval2")}

  survfitted <- survfit(survdata ~ 1)

  emp = data.frame(xs = survfitted$time, ps = 1-survfitted$surv)


  res = fit_list %>%
    mclapply(function(x) {
      get_cdf(x, xs = xs) %>%
        data.frame(
          xs = xs, type = x$qmethod, ps = .
        )
    }, mc.cores = ifelse(par_, yes = floor(detectCores()), no = 1)) %>%
    Reduce(rbind,.) %>%
    mutate(type =  factor(type,levels=sort(as.character( unique(type) ))) )


  if(n_splitted_plot>1) {

    create_cv_fun = function(res){
      types = res$type %>% unique()
      cv_array = mapply(FUN = function(x,y) y, types, 1:n_splitted_plot) %>% setNames(types)
      function(x) cv_array[as.numeric(x)]
    }

    cv_fun = create_cv_fun(res)

    res = res %>%
      group_by(type) %>%
      mutate(plot_id = cv_fun(type))

  }

  plot_ =  res %>%
    ggplot(aes(x = xs, y = ps, colour = type)) +
    geom_line(size = 1.5) +
    geom_step(data = emp, aes(x = xs, y = ps, colour = NULL), size = 1.5) +
    # annotate('segment', x = 0, xend = 3, y= 0, yend = 3, colour = 'red')
    ylab('') +
    xlab('Data')

  if(n_splitted_plot>1) plot_ = plot_ + facet_wrap(~plot_id)




  if(is.censored(fit_list[[1]]$data_)){
    plot_ + ggtitle('fitted CDF and Turnbull estimator')
  }

  else{plot_ + ggtitle('fitted CDF and empirical estimator')}

}

#Test CDF function

# rnorm(n = 10**3, mean = 40, sd = 6) %>%
#   (function(x){
#     list(fit_normal_ssd_on_log_data,
#          fit_logis_ssd_on_log_data,
#          fit_kernel_mixture_ssd_on_log_data,
#          fit_DP_mixture_ssd_on_log_data,
#          function(x) fit_BNP_mixture_ssd_on_log_data(x, Nit = 200)
#          ) %>%
#       lapply(FUN = function(fitfun) fitfun(x))
#   }) %>%
#   comparison_CDFplot()


comparison_qqplot = function(fit_list, ps = seq(0.01,0.99, length.out = 100)) {
  empirical_quantiles = fit_list[[1]]$data_ %>% quantile(probs = ps) %>% unlist()

  fit_list %>%
    lapply(function(x) {
      get_quantiles(x, ps) %>%
        data.frame(
          ps = ps, type = x$qmethod, qs = ., empirical_quantiles = empirical_quantiles
        )
    }) %>%
    Reduce(rbind,.) %>%
    mutate(model =  factor(model,levels=sort(as.character( unique(model) ))) ) %>%
    ggplot(aes(x = empirical_quantiles, y = qs, colour = type)) +
    geom_point(size = 4) +
    # annotate('segment', x = 0, xend = 3, y= 0, yend = 3, colour = 'red')
    geom_abline(intercept = 0, slope = 1) +
    ylab('Fitted quatiles') +
    xlab('Empirical quantiles')

}


comparison_densplot = function(fit_list, xs = NULL, par_ = T,
                               logsc = F, n_splitted_plot = 4) {

  if(is.null(xs)){
    xs = fit_list[[1]]$data_ %>%
      make_grid_from_data(length_ = 100)
  }

  res = fit_list %>%
    mclapply(function(x) {
      get_dens(fit = x, xs = xs) %>%
        data.frame(
          xs = xs, type = x$qmethod, dens = .)
    }, mc.cores = ifelse(par_, yes = floor(detectCores()), no = 1)) %>%
    Reduce(rbind,.) %>%
    mutate(type =  factor(type,levels=sort(as.character( unique(type) ))) )

  if(n_splitted_plot>1) {

    create_cv_fun = function(res){
      types = res$type %>% unique()
      cv_array = mapply(FUN = function(x,y) y, types, 1:n_splitted_plot) %>% setNames(types)
      function(x) cv_array[as.numeric(x)]
    }

    cv_fun = create_cv_fun(res)

    res = res %>%
      group_by(type) %>%
      mutate(plot_id = cv_fun(type))
  }



  if(fit_list[[1]]$data_ %>% is.censored){
    p = data.frame(dat = fit_list[[1]]$data_) %>%
      ggplot(aes(x = dat)) +
      geom_line(data = res, aes(x = xs, y = dens, colour = type), size = 1.5) +
      ylab('Density') +
      xlab('Data')
  }
  else{
    p = data.frame(dat = fit_list[[1]]$data_) %>%
      ggplot(aes(x = dat)) +
      geom_histogram(aes(y = ..density..)) +
      geom_line(data = res, aes(x = xs, y = dens, colour = type), size = 1.5) +
      ylab('Density') +
      xlab('Data')
  }

  if(logsc) p = p + scale_y_log10()

  if(n_splitted_plot>1) p = p + facet_wrap(~plot_id, scales = 'free_y')

  p = p +
    theme_few() +
    scale_colour_few()

  return(p)

}

comparison_plot = function(fit_list, ps = seq(0.01,0.99, length.out = 100), xs = seq(-5,5, length.out = 100)){

  if(is.censored(fit_list[[1]]$data_)){
    grid.arrange( fit_list %>%
                    comparison_CDFplot(xs = xs))
  }
  else{
    grid.arrange(fit_list %>%
                   comparison_qqplot(ps = ps),
                 fit_list %>%
                   comparison_densplot(xs = xs),
                 fit_list %>%
                   comparison_CDFplot(xs = xs))
  }
}

censor_code_rl = function(left, right) {
  test_ = function(k) {
    if (is.na(left[[k]]) & is.na(right[[k]]))
      NA #for safety
    else if (is.na(left[[k]]))
      2 #left-censored
    else if (is.na(right[[k]]))
      0 #right-censored
    else if (left[[k]] == right[[k]])
      1 #non-censored
    else
      3 #interval-censored
  }
  sapply(seq_along(left), FUN = test_)
}

compute_cpo = function(fit, par_ = T) {

  ncores = ifelse(par_, yes = detectCores(), no = 1)

  if (fit$qmethod %>% grepl('DP_mixture_ssd',.)) {
    select_pars = function(parname = 'mu') {
      if (parname %in% c('m1', 'k0', 'psi1', 'ncluster', 'alpha')) {
        fit$save.state$thetasave[1,] %>%
          names() %>%
          grepl(pattern = parname, x = .)
      }
      #Acceptable parnames: mu, sigma, m1, k0, psi1, ncluster, alpha
      else if (parname %in% c('mu', 'sigma')) {
        fit$save.state$randsave[1,] %>%
          names() %>%
          grepl(pattern = parname, x = .)
      }
      else
        stop('no param with this name')
    }

    dpred = function(iter) {
      function(y) {
        1 / (fit$prior$alpha + fit$nrec) * (
          fit$prior$alpha * dscaled_student(
            x = y,
            df = 2 * fit$prior$nu1,
            ncp = fit$save.state$thetasave[iter,select_pars('m1')],
            scale = sqrt(
              (1 + 1 / fit$save.state$thetasave[iter,select_pars('k0')]) /
                (fit$save.state$thetasave[iter,select_pars('psi1')] * fit$prior$nu1)
            )
          ) +
            dmixnorm(
              y ,
              mus = fit$save.state$randsave[iter,select_pars('mu')][1:fit$nrec],
              sigmas = fit$save.state$randsave[iter,select_pars('sigma')][1:fit$nrec] %>% sqrt,
              probs = rep(1, fit$nrec)
            )
        )
      }
    }
    set.seed(0)
    itsample = 1:fit$mcmc$nsave
    d_list = itsample %>% lapply(dpred)
  }

  else if (fit$qmethod %>% grepl('BNP_mixture_ssd',.)) {
    dpred = function(iter) {
      function(y) {
        dmixnorm(
          y ,
          mus = fit$means[[iter]],
          sigmas = fit$sigmas[[iter]],
          probs = fit$weights[[iter]]
        )
      }
    }
    d_list = seq(fit$means) %>% lapply(dpred)

  }

  if (is.censored(fit$data_)) {
    if (fit$qmethod %>% grepl('DP_mixture_ssd',.)) {
      ppred = function(iter) {
        function(y) {
          1 / (fit$prior$alpha + fit$nrec) * (
            fit$prior$alpha * pscaled_student(
              x = y,
              df = 2 * fit$prior$nu1,
              ncp = fit$save.state$thetasave[iter,select_pars('m1')],
              scale = sqrt(
                (1 + 1 / fit$save.state$thetasave[iter,select_pars('k0')]) /
                  (fit$save.state$thetasave[iter,select_pars('psi1')] *
                     fit$prior$nu1)
              )
            ) +
              pmixnorm(
                y ,
                mus = fit$save.state$randsave[iter,select_pars('mu')][1:fit$nrec],
                sigmas = fit$save.state$randsave[iter,select_pars('sigma')][1:fit$nrec] %>% sqrt,##Yeah, well, seems like its a sigma**2 in the end
                probs = rep(1, fit$nrec)
              )
          )
        }
      }

      cdF_list = itsample %>% lapply(ppred)

    }

    else if (fit$qmethod %>% grepl('BNP_mixture_ssd',.)) {
      ppred = function(iter) {
        function(y) {
          pmixnorm(
            y ,
            mus = fit$means[[iter]],
            sigmas = fit$sigmas[[iter]],
            probs = fit$weights[[iter]]
          )
        }
      }


      cdF_list = seq(fit$means) %>% lapply(ppred)

    }

    pred_lik_cens = function(xleft, xright, iter) {
      c_code = censor_code_rl(xleft, xright)

      res = seq_along(c_code) #initialisation

      if(any(c_code == 1)) res[c_code == 1] = d_list[[iter]](xleft[c_code == 1])
      if(any(c_code == 2)) res[c_code == 2] = cdF_list[[iter]](xright[c_code == 2])
      if(any(c_code == 0)) res[c_code == 0] = 1 - cdF_list[[iter]](xleft[c_code == 0])
      if(any(c_code == 3)) res[c_code == 3] = cdF_list[[iter]](xright[c_code == 3]) - cdF_list[[iter]](xleft[c_code == 3])

      res
    }

    fit$cpo = seq_along(d_list) %>%
      mclapply(function(pred_lik_cens_j) {
        pred_lik_cens(fit$data_$left, fit$data_$right, pred_lik_cens_j)
      }, mc.cores = ncores) %>%
      Reduce(cbind,.) %>%
      rowMeans

  }

  else{
    fit$cpo = d_list %>%
      mclapply(function(f_j) {
        f_j(fit$data_)
      }, mc.cores = ncores) %>%
      Reduce(cbind,.) %>%
      (function(x) 1/x) %>%
      rowMeans %>%
      (function(x) 1/x)

  }

  return(fit)
}

compare_quantiles_and_CI = function(fit_list, p = 0.05, numit = 500, nbootsamples = 10**4){
  fit_list %>%
    lapply(FUN = function(x){
      if(x$qmethod %>% grepl('BNP',.)) thin_fit(x, numit)
      else x
    }) %>%
    lapply(FUN = function(fit){
      fit %>%
        get_quantiles_and_CI(p = p, nbootsamples = nbootsamples) %>%
        mutate(model = fit$qmethod)
    }) %>%
    Reduce(rbind,.)
}

plot_compare_cpo_nc = function(fit_listnc, plot_outliers = T) {
  data_to_plot = fit_listnc %>%
    lapply(
      FUN = function(x) {
        data.frame(x = x$data_, cpo = x$cpo, model = x$qmethod)
      }
    ) %>%
    Reduce(rbind,.) %>%
    mutate(model =  factor(model,levels = sort(as.character(unique(
      model
    )))))


  plot_ = data_to_plot %>%
    ggplot(aes(y = cpo, x = model, fill = model)) +
    geom_boxplot() +
    xlab('') +
    scale_x_discrete(breaks = NULL)

  if(plot_outliers) plot_ + scale_y_log10()
  else {
    ylim1 = boxplot.stats(data_to_plot$cpo)$stats[c(1, 5)]
    ylim1 = ylim1 + c(-0.05, 0.05) * diff(ylim1) / 2
    plot_ +
      # scale_y_log10(limits = ylim1*1.05) +
      scale_y_log10(limits = ylim1) +
      ggtitle('Boxplot without the outliers')
  }

}

plot_compare_given_quantiles = function(fit_list, numit = 500, nbootsamples = 10**4) {
  fit_list %>%
    compare_quantiles_and_CI(p = c(0.05), numit = numit, nbootsamples = nbootsamples) %>%
    ggplot(aes(
      x = model, y = HC, ymin = C_inf, ymax = C_sup, color = model
    )) +
    geom_errorbar() +
    geom_point(size = 3.5, colour = 'white') +
    geom_point(size = 2.5) +
    xlab('') +
    ylab('5th percentile')  +
    scale_x_discrete(breaks = NULL) +
    theme_few() +
    scale_colour_few()
}

plot_compare_cpo_cens = function(cpo, plot_outliers = T) {
  data_to_plot = cpo %>%
    mutate(cens_code = censor_code_rl(left = x.left, right = x.right)) %>%
    mutate(type = ifelse(cens_code == 1, 'Non-censored', 'Censored')) %>%
    (function(x) {
      x %>% rbind(x,x %>% mutate(type = 'Both'))
    }) %>%
    mutate(model =  factor(model,levels = sort(as.character(unique(
      model
    )))))

  plot_ = data_to_plot %>% ggplot(aes(
    y = cpo, x = factor(type), fill = model
  )) +
    geom_boxplot() +
    xlab('Type of data')

  if(plot_outliers) plot_ + scale_y_log10()
  else {
    ylim1 = boxplot.stats(data_to_plot$cpo)$stats[c(1, 5)]
    ylim1 = ylim1 + c(-0.05, 0.05) * diff(ylim1) / 2
    plot_ +
      # scale_y_log10(limits = ylim1*1.05) +
      scale_y_log10(limits = ylim1) +
      ggtitle('Boxplot without the outliers')
  }


}

thin_fit = function(fit, numit = 10**6){
  if(fit$qmethod %>% grepl('DP_mixture_ssd',.)){
    stop('do it if we need it some day')
  }

  else if (fit$qmethod %>% grepl('BNP_mixture_ssd',.)|fit$qmethod %>% grepl('DP_mixt_ssdBNPdens',.)) {
    selected_iterations = seq(from = 1, to = length(fit$means), by = ceiling(length(fit$means)/numit))

    fit$means = fit$means[selected_iterations]
    fit$sigmas = fit$sigmas[selected_iterations]
    fit$weights = fit$weights[selected_iterations]

    return(fit)
  }

  else fit
}

plot_BNPdensity_fit = function(fit, xs = NULL){
  if(is.null(xs)){
    xs = fit$data_ %>%
      make_grid_from_data(length_ = 100)
  }

  data.frame(x = fit$data_) %>%
    ggplot(aes(x = x, y = ..density..)) +
    geom_histogram() +
    geom_line(data = data.frame(xs = xs, dens = fit %>% get_dens(xs)), aes(x = xs, y = dens))

}

plot_BNPdensity_CDFfit = function(fit, xs = NULL, names = F, clustering = F) {

  if(is.null(xs)){
    xs = fit$data_ %>%
      make_grid_from_data(length_ = 100, Delta = 1)
  }

  if(names) xs = c(xs, max(xs) + 2)

  if(is.censored(fit_list[[1]]$data_)){
    survdata <- Surv(time = fit$data_$left, time2 = fit_list[[1]]$data_$right,
                     type = "interval2")
  }

  else{survdata <- Surv(time = fit$data_, time2 = fit$data_,
                        type = "interval2")}

  survfitted <- survfit(survdata ~ 1)

  emp = data.frame(xs = survfitted$time, ps = 1-survfitted$surv)

  plot_ =  fit %>%
    (function(x){
      get_cdf(x, xs = xs) %>%
        data.frame(
          xs = xs, type = x$qmethod, ps = .
        )
    }) %>%
    ggplot(aes(x = xs, y = ps)) +
    geom_line(size = 1.5, colour = 'red') +
    geom_step(data = emp, aes(x = xs, y = ps, colour = NULL), size = 1.) +
    # annotate('segment', x = 0, xend = 3, y= 0, yend = 3, colour = 'red')
    ylab('') +
    xlab('Data')

  if(clustering){
    library(mcclust.ext)
    fit.draw = fit$Allocs %>% Reduce(rbind,.)
    psm=comp.psm(fit.draw)
    fit.VI=minVI(psm, cls.draw = fit.draw, method=("greedy"))

    plot_ = plot_ +
      geom_point(data = data.frame(x = fit$data_, y = get_cdf(fit, xs = fit$data_), clust = fit.VI$cl),
                 aes(x = x, y = y, colour = factor(clust)), size = 1.5) +
      scale_colour_manual(name = 'cluster index', values = c("blue", "green"))
  }

  if(names){
    if(clustering){
      plot_ = plot_ +
        geom_text(data = data.frame(x = fit$data_,
                                    y = get_cdf(fit, xs = fit$data_),
                                    text = fit$data_ %>% names,
                                    clust = fit.VI$cl),
                  aes(x = x, y = y, label = text,
                      colour = factor(clust)),
                  nudge_x = 2, check_overlap = F, size = 3
        ) +
        scale_colour_manual(name = 'cluster index', values = c("blue", "green"))
    }
    else{
      plot_ = plot_ +
        geom_text(data = data.frame(x = fit$data_, y = get_cdf(fit, xs = fit$data_), text = fit$data_ %>% names),
                  aes(x = x, y = y, label = text),
                  nudge_x = 2, check_overlap = F, size = 3
        )
    }

  }

  if(is.censored(fit_list[[1]]$data_)){
    plot_ + ggtitle('fitted CDF and Turnbull estimator')
  }

  else{plot_ + ggtitle('fitted CDF and empirical estimator')}

}
