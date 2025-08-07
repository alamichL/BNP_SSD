### All the functions to fit an SSD to data
library(dplyr)
library(tidyr)
library(parallel)
library(ggplot2)
library(fitdistrplus)
library(gridExtra)
library(BNPdensity)
library(compiler)


centre_and_scale <- function(ddat) {
  if (is.censored(ddat)) {
    c_sc <- ddat %>%
      (function(x) 0.5 * (x$left + x$right)) %>%
      (function(x) c("mean" = mean(x, na.rm = T), "sd" = sd(x, na.rm = T)))
    return(list(data_ = (ddat - c_sc[["mean"]]) / c_sc[["sd"]], mean = c_sc[["mean"]], sd = c_sc[["sd"]]))
  } else {
    return(list(data_ = (ddat - mean(ddat, na.rm = T)) / sd(ddat, na.rm = T), mean = mean(ddat, na.rm = T), sd = sd(ddat, na.rm = T)))
  }
}

is.censored <- function(ddat) {
  if (is.null(ncol(ddat))) {
    FALSE
  } else if (ncol(ddat) == 1) {
    FALSE
  } else if (ncol(ddat) == 2) {
    TRUE
  } else {
    stop("Wrong type/dim of data input")
  }
}

## Parametric SSD

### Log-normal SSD

fit_normal_ssd_on_log_data <- function(ddat) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_normal_ssd_on_log_data_dont_scale(centred_and_scaled_dd$data_) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}


fit_normal_ssd_on_log_data_dont_scale <- function(ddat) {
  if (is.censored(ddat)) {
    fitdistcens(
      censdata = data.frame(
        left = ddat[, 1] %>%
          ifelse(test = is.nan(.), yes = NA, no = .),
        right = ddat[, 2] %>%
          ifelse(test = is.nan(.), yes = NA, no = .)
      ),
      distr = "norm"
    ) %>%
      (function(x) {
        x$qmethod <- "normal_ssd"
        x$data_ <- ddat
        x %>% return()
      })
  } else {
    fitdist(data = ddat, distr = "norm") %>%
      (function(x) {
        x$qmethod <- "normal_ssd"
        x$data_ <- ddat
        x %>% return()
      })
  }
}

# Exemple:
#
# ex_c_list %>%
#   lapply(FUN = get_log_dat) %>%
#   lapply(FUN = fit_normal_ssd_on_log_data) %>%
#   lapply(FUN = plot)


### Log-logistic SSD

fit_logis_ssd_on_log_data <- function(ddat) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_logis_ssd_on_log_data_dont_scale(centred_and_scaled_dd$data_) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}

fit_logis_ssd_on_log_data_dont_scale <- function(ddat) {
  if (is.censored(ddat)) {
    fitdistcens(
      censdata = data.frame(
        left = ddat[, 1] %>%
          ifelse(test = is.nan(.), yes = NA, no = .),
        right = ddat[, 2] %>%
          ifelse(test = is.nan(.), yes = NA, no = .)
      ),
      distr = "logis"
    ) %>%
      (function(x) {
        x$qmethod <- "logistic_ssd"
        x$data_ <- ddat
        x %>% return()
      })
  } else {
    fitdist(data = ddat, distr = "logis") %>%
      (function(x) {
        x$qmethod <- "logistic_ssd"
        x$data_ <- ddat
        x %>% return()
      })
  }
}

# Exemple:
#
# ex_c_list %>%
#   lapply(FUN = get_log_dat) %>%
#   lapply(FUN = fit_logis_ssd_on_log_data) %>%
#   lapply(FUN = plot)
#

## Non parametric SSD

# Kernel non parametric

# Wang uses normal kernel with (asymptotic mean integrated squared error, amise) bandwith given by $h_n = 1.06 \hat{\sigma} n^{-\frac{1}{5}}$ (Silverman 1986) où je suppose que $\hat{\sigma}$ est la variance estimée à partir de la variance empirique.

fit_kernel_mixture_ssd_on_log_data <- function(ddat) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_kernel_mixture_ssd_on_log_data_dont_scale(ddat = centred_and_scaled_dd$data_) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}

fit_kernel_mixture_ssd_on_log_data_dont_scale <- function(ddat) {
  # if(is.censored(ddat))    stop('no kernel mixture for censored data yet')
  # else {
  #
  #   K = function(x){
  #     1/sqrt(2*pi) * exp(-x^2/2)
  #   }
  #
  #   n = length(ddat)
  #
  #   hn = 1.06 * sd(ddat) * n^(-1/5)
  #
  #   f_nonvec = function(x){
  #     1/(n*hn) * sum(sapply(X = ddat, FUN = function(xi) K((x-xi)/hn)))
  #   }
  #
  #   return(function(x) sapply(x, f_nonvec))
  # }
  fit_kernel_mixture_ssd_on_log_data_make_object(ddat)
}

source("scripts/mixnorm.R")

fit_kernel_mixture_ssd_on_log_data_make_object <- function(ddat) {
  if (is.censored(ddat)) {
    stop("no kernel mixture for censored data yet")
  } else {
    res <- list()

    K <- function(x) {
      1 / sqrt(2 * pi) * exp(-x^2 / 2)
    }



    n <- length(ddat)

    hn <- 1.06 * sd(ddat) * n^(-1 / 5)

    dmixkernel_nonvec <- function(x) {
      1 / (n * hn) * sum(sapply(X = ddat, FUN = function(xi) K((x - xi) / hn)))
    }

    pmixkernel_nonvec <- function(x) {
      1 / (n) * sum(sapply(X = ddat, FUN = function(xi) pnorm(q = (x - xi) / hn, mean = 0, sd = 1)))
    }

    rmixkernel <- function(ndata) {
      rmixnorm(n = ndata, mus = ddat, sigmas = hn, probs = 1 / n * rep(1, n))
    }

    res %>%
      (function(x) {
        x$dmixkernel <- function(x) sapply(x, dmixkernel_nonvec)
        x$pmixkernel <- function(x) sapply(x, pmixkernel_nonvec)
        x$rmixkernel <- rmixkernel
        x$qmethod <- "normal_kernel_mixture_ssd"
        x$data_ <- ddat
        x %>%
          return()
      }) %>%
      return()
  }
}


# Exemple
#
# ex_c_list %>%
#   lapply(get_log_dat) %>%
#   lapply(function(x){
#     kernel_ssd = x %>% fit_kernel_mixture_ssd_on_log_data
#     hist(x, freq = F)
#
#     seq(min(x)-1, max(x)+1, length.out = 1000) %>%
#       lines(., kernel_ssd(.))
#
#   })



# BNP-SSD

## BNPdensity


source("scripts/MixNRMI1dens_mvcpp.R")
source("scripts/MixNRMI2dens_mvcpp.R")

source("scripts/MixNRMI1_2denscens_mvcpp_faster.R")
### DP mixture

fit_semi_DP_mixture_ssd_on_log_data <- function(ddat, Nit = 100, Pbi = 0.2, label = "", ...) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_semi_DP_mixture_ssd_on_log_data_dont_scale(ddat = centred_and_scaled_dd$data_, Nit = Nit, Pbi = Pbi, label = label, ...) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}

fit_semi_DP_mixture_ssd_on_log_data_dont_scale <- function(ddat, Nit = 100, Pbi = 0.2, label = "", ...) {
  if (is.censored(ddat)) {
    MixNRMI1denscens_mvcpp_faster(
      xleft = ddat$left, xright = ddat$right,
      Beta = 1, Gama = 0, # To get DP mixture
      Nit = Nit, Pbi = Pbi, ...
    ) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("semi_DP_mixt_ssdBNPdens", label, sep = "_")
        x
      }) %>%
      return()
  } else {
    MixNRMI1dens_mvcpp(ddat,
      Beta = 1, Gama = 0, # To get DP mixture
      Nit = Nit, Pbi = Pbi, ...
    ) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("semi_DP_mixt_ssdBNPdens", label, sep = "_")
        x
      }) %>%
      return()
  }
}

fit_DP_mixture_ssd_on_log_data <- function(ddat, mu.pz0 = 1, sigma.pz0 = 1, Nit = 100, Pbi = 0.2, label = "", ...) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_DP_mixture_ssd_on_log_data_dont_scale(ddat = centred_and_scaled_dd$data_, mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0, Nit = Nit, Pbi = Pbi, label = label, ...) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}


fit_DP_mixture_ssd_on_log_data_dont_scale <- function(ddat, mu.pz0 = 1, sigma.pz0 = 1, Nit = 100, Pbi = 0.2, label = "", ...) {
  if (is.censored(ddat)) {
    MixNRMI2denscens_mvcpp_faster(
      xleft = ddat$left, xright = ddat$right,
      mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0,
      Beta = 1, Gama = 0, # To get DP mixture
      Nit = Nit, Pbi = Pbi, ...
    ) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("DP_mixt_ssdBNPdens", label, sep = "_")
        x
      }) %>%
      return()
  } else {
    MixNRMI2dens_mvcpp(ddat,
      mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0,
      Beta = 1, Gama = 0, # To get DP mixture
      Nit = Nit, Pbi = Pbi, ...
    ) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("DP_mixt_ssdBNPdens", label, sep = "_")
        x
      }) %>%
      return()
  }
}


### BNP mixture

fit_semi_BNP_mixture_ssd_on_log_data <- function(ddat, Nit = 1000, Pbi = 0.2, label = "", ...) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_semi_BNP_mixture_ssd_on_log_data_dont_scale(ddat = centred_and_scaled_dd$data_, Nit = Nit, Pbi = Pbi, label = label, ...) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}

fit_semi_BNP_mixture_ssd_on_log_data_dont_scale <- function(ddat, Nit = 1000, Pbi = 0.2, label = "", ...) {
  if (is.censored(ddat)) {
    MixNRMI1denscens_mvcpp_faster(xleft = ddat$left, xright = ddat$right, Nit = Nit, Pbi = Pbi, ...) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("semi_BNP_mixture_ssd", label, sep = "_")
        x
      }) %>%
      return()
  } else {
    MixNRMI1dens_mvcpp(ddat, Nit = Nit, Pbi = Pbi, ...) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("semi_BNP_mixture_ssd", label, sep = "_")
        x
      }) %>%
      return()
  }
}

fit_BNP_mixture_ssd_on_log_data <- function(ddat, mu.pz0 = 1, sigma.pz0 = 1, Nit = 100000, Pbi = 0.2, label = "", ...) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_BNP_mixture_ssd_on_log_data_dont_scale(ddat = centred_and_scaled_dd$data_, mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0, Nit = Nit, Pbi = Pbi, label = label, ...) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}


fit_BNP_mixture_ssd_on_log_data_dont_scale <- function(ddat, mu.pz0 = 1, sigma.pz0 = 1, Nit = 100000, Pbi = 0.2, label = "", ...) {
  if (is.censored(ddat)) {
    MixNRMI2denscens_mvcpp_faster(xleft = ddat$left, xright = ddat$right, mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0, Nit = Nit, Pbi = Pbi, ...) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("BNP_mixture_ssd", label, sep = "_")
        x
      }) %>%
      return()
  } else {
    MixNRMI2dens_mvcpp(ddat, mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0, Nit = Nit, Pbi = Pbi, ...) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("BNP_mixture_ssd", label, sep = "_")
        x
      }) %>%
      return()
  }
}

source("scripts/MixNRMI2dens_no_resamp.R")
source("scripts/MixNRMI1_2denscens_no_resamp.R")


fit_on_censored_or_noncensored <- function(ddat, label, ...) {
  if (is.censored(ddat)) {
    MixNRMI2cens(xleft = ddat$left, xright = ddat$right, extras = T, ...) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("BNP_mixture_ssd", label, sep = "_")
        x
      }) %>%
      return()
  } else {
    MixNRMI2cens(xleft = ddat, xright = ddat, extras = T, ...) %>%
      (function(x) {
        x$data_ <- ddat
        x$qmethod <- paste("BNP_mixture_ssd", label, sep = "_")
        x
      }) %>%
      return()
  }
}



fit_BNP_mixture_ssd_on_log_data_truncated_normal <- function(ddat, mu.pz0 = 1, sigma.pz0 = 1,
                                                             Nit = 1000, Pbi = 0.2, ...) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_BNP_mixture_ssd_on_log_data_truncated_normal_no_scale(
    ddat = centred_and_scaled_dd$data_, mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0,
    Nit = Nit, Pbi = Pbi, ...
  ) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}


fit_BNP_mixture_ssd_on_log_data_truncated_normal_no_scale <- function(ddat, mu.pz0 = 1, sigma.pz0 = 1,
                                                                      Nit = 1000, Pbi = 0.2, ...) {
  fit_on_censored_or_noncensored(ddat,
    mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0, distr.pz0 = 10,
    Nit = Nit, Pbi = Pbi, label = paste("mu", mu.pz0,
      "_sigma", sigma.pz0,
      "_truncated_normal_Gama0.4",
      sep = ""
    ), Gama = 0.4
  )
}



fit_BNP_mixture_ssd_on_log_data_uniform <- function(ddat, mu.pz0 = 0.1, sigma.pz0 = 1.5,
                                                    Nit = 1000, Pbi = 0.2, ...) {
  original_data <- ddat
  centred_and_scaled_dd <- centre_and_scale(ddat)
  fit_BNP_mixture_ssd_on_log_data_uniform_no_scale(
    ddat = centred_and_scaled_dd$data_, mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0,
    Nit = Nit, Pbi = Pbi, ...
  ) %>%
    (function(x) {
      x$data_ <- original_data
      x$centre <- centred_and_scaled_dd$mean
      x$scale <- centred_and_scaled_dd$sd
      x %>% return()
    })
}

fit_BNP_mixture_ssd_on_log_data_uniform_no_scale <- function(ddat, mu.pz0 = 0.1, sigma.pz0 = 1.5,
                                                             Nit = 1000, Pbi = 0.2, ...) {
  fit_on_censored_or_noncensored(ddat,
    mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0, distr.pz0 = 9,
    Nit = Nit, Pbi = Pbi, label = paste("mu", mu.pz0,
      "_sigma", sigma.pz0,
      "_uniform_Gama0.4",
      sep = ""
    ), Gama = 0.4, ...
  )
}



Var <- function(x, ...) {
  y <- length(x)
  var(x, ...) * (y - 1) / y
}


Sd <- function(x, ...) {
  y <- length(x)
  sd(x, ...) * sqrt((y - 1) / y)
}
