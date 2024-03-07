get_dens_DP_package = function(fit, xs = seq(-5,5, length.out = 100)){

  select_pars = function(parname = 'mu') {
    if (parname %in% c('m1', 'k0', 'psi1', 'ncluster', 'alpha')) {
      fit$save.state$thetasave[1,] %>% names() %>% grepl(pattern = parname, x =
                                                           .)
    }
    #Acceptable parnames: mu, sigma, m1, k0, psi1, ncluster, alpha
    else if (parname %in% c('mu', 'sigma')) {
      fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = parname, x =
                                                          .)
    }
    else
      stop('no param with this name')
  }

  select_sigmas = function() {
    fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = 'sigma', x =
                                                        .)
  }

  if(!is.null(fit$prior$alpha)){
    fit$save.state$thetasave = fit$save.state$thetasave %>%
      cbind(alpha = rep(fit$prior$alpha, fit$mcmc$nsave))
  }

  dpred = function(iter) {
    function(y) {
      1 / (fit$save.state$thetasave[iter,select_pars('alpha')] + fit$nrec) * (
        fit$save.state$thetasave[iter,select_pars('alpha')] * dscaled_student(
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


  d_list = seq_len(fit$mcmc$nsave) %>% lapply(dpred)

  d_list %>%
    mclapply(function(f_j)
      f_j(xs), mc.cores = detectCores()) %>%
    Reduce(cbind,.) %>%
    rowMeans

}


library(Rcpp)

sourceCpp('scripts/col_rowMedianRcpp.cpp')

get_median_dens_DP_package = function(fit, xs = seq(-5,5, length.out = 100)){

  select_pars = function(parname = 'mu') {
    if (parname %in% c('m1', 'k0', 'psi1', 'ncluster', 'alpha')) {
      fit$save.state$thetasave[1,] %>% names() %>% grepl(pattern = parname, x =
                                                           .)
    }
    #Acceptable parnames: mu, sigma, m1, k0, psi1, ncluster, alpha
    else if (parname %in% c('mu', 'sigma')) {
      fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = parname, x =
                                                          .)
    }
    else
      stop('no param with this name')
  }

  select_sigmas = function() {
    fit$save.state$randsave[1,] %>% names() %>% grepl(pattern = 'sigma', x =
                                                        .)
  }

  if(!is.null(fit$prior$alpha)){
    fit$save.state$thetasave = fit$save.state$thetasave %>%
      cbind(alpha = rep(fit$prior$alpha, fit$mcmc$nsave))
  }

  dpred = function(iter) {
    function(y) {
      1 / (fit$save.state$thetasave[iter,select_pars('alpha')] + fit$nrec) * (
        fit$save.state$thetasave[iter,select_pars('alpha')] * dscaled_student(
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


  d_list = seq_len(fit$mcmc$nsave) %>% lapply(dpred)

  d_list %>%
    mclapply(function(f_j)
      f_j(xs), mc.cores = detectCores()) %>%
    Reduce(cbind,.) %>%
    rowMedianRcpp()

}
