library(BNPdensity)
source('scripts/MvInv_cpp.R')

censor_code_rl = function(left, right) {
  test_ = function(k) {
    if (is.na(left[[k]]) & is.na(right[[k]]))
      NA #for safety
    if (is.na(left[[k]]))
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


rk = function (n, distr = NULL, mu = NULL, sigma = NULL)
{
  if (is.null(distr)) {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
  }
  else if (distr == 1) {
    a <- ifelse(is.null(mu), 0, mu)
    b <- ifelse(is.null(sigma), 1, sigma)
    rk <- rnorm(n, mean = a, sd = b)
  }
  else if (distr == 4) {
    a <- ifelse(is.null(mu), 0, mu)
    b <- ifelse(is.null(sigma), 1 / sqrt(2), sigma / sqrt(2))
    rk <- a + b * sample(c(-1,+1), size = n, replace = TRUE) *
      rexp(n)
  }
  else if (distr == 5) {
    a <- ifelse(is.null(mu), exp(1 / 2), log(mu / sqrt(1 + (sigma / mu) ^ 2)))
    b <- ifelse(is.null(sigma), exp(1) * (exp(1) - 1), sqrt(log(1 +
                                                                  (sigma /
                                                                     mu) ^ 2)))
    rk <- rlnorm(n, meanlog = a, sdlog = b)
  }
  else if (distr == 2) {
    a <- ifelse(is.null(mu), 1, mu ^ 2 / sigma ^ 2)
    b <- ifelse(is.null(sigma), 1, mu / sigma ^ 2)
    rk <- rgamma(n, shape = a, rate = b)
  }
  else if (distr == 3) {
    a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu / sigma) ^ 2 -
                  mu)
    b <-
      ifelse(is.null(sigma), 1 / sqrt(12), (mu * (1 - mu) / sigma ^ 2 -
                                              1) * (1 - mu))
    if (any(c(a, b) <= 0))
      stop(paste("\nNegative Beta parameters:\n a =", a,
                 ";\t b =", b))
    rk <- rbeta(n, shape1 = a, shape2 = b)
  }
  else {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
  }
  return(rk)
}


dk = function (x, distr = NULL, mu = NULL, sigma = NULL)
{
  if (is.null(distr)) {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
  }
  else if (distr == 1) {
    a <- ifelse(is.null(mu), 0, mu)
    b <- ifelse(is.null(sigma), 1, sigma)
    dk <- dnorm(x, mean = a, sd = b)
  }
  else if (distr == 4) {
    a <- ifelse(is.null(mu), 0, mu)
    b <- ifelse(is.null(sigma), 1 / sqrt(2), sigma / sqrt(2))
    dk <- exp(-abs(x - a) / b) / (2 * b)
  }
  else if (distr == 5) {
    a <- ifelse(is.null(mu), exp(1 / 2), log(mu / sqrt(1 + (sigma / mu) ^ 2)))
    b <- ifelse(is.null(sigma), exp(1) * (exp(1) - 1), sqrt(log(1 +
                                                                  (sigma /
                                                                     mu) ^ 2)))
    dk <- dlnorm(x, meanlog = a, sdlog = b)
  }
  else if (distr == 2) {
    a <- ifelse(is.null(mu), 1, mu ^ 2 / sigma ^ 2)
    b <- ifelse(is.null(sigma), 1, mu / sigma ^ 2)
    dk <- dgamma(x, shape = a, rate = b)
  }
  else if (distr == 3) {
    a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu / sigma) ^ 2 -
                  mu)
    b <-
      ifelse(is.null(sigma), 1 / sqrt(12), (mu * (1 - mu) / sigma ^ 2 -
                                              1) * (1 - mu))
    if (any(c(a, b) <= 0))
      stop(paste("\nNegative Beta parameters:\n a =", a,
                 ";\t b =", b))
    dk <- dbeta(x, shape1 = a, shape2 = b)
  }
  else {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
  }
  return(dk)
}

pk = function (q, distr = NULL, mu = NULL, sigma = NULL)
{
  if (is.null(distr)) {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
  }
  else if (distr == 1) {
    a <- ifelse(is.null(mu), 0, mu)
    b <- ifelse(is.null(sigma), 1, sigma)
    pk <- pnorm(q, mean = a, sd = b)
  }
  else if (distr == 2) {
    stop('pk for exponential not yet done')
    #         a <- ifelse(is.null(mu), 0, mu)
    #         b <- ifelse(is.null(sigma), 1/sqrt(2), sigma/sqrt(2))
    #         pk <- exp(-abs(x - a)/b)/(2 * b)
  }
  else if (distr == 3) {
    a <- ifelse(is.null(mu), exp(1 / 2), log(mu / sqrt(1 + (sigma / mu) ^ 2)))
    b <- ifelse(is.null(sigma), exp(1) * (exp(1) - 1), sqrt(log(1 +
                                                                  (sigma /
                                                                     y) ^ 2)))
    pk <- plnorm(q, meanlog = a, sdlog = b)
  }
  else if (distr == 4) {
    a <- ifelse(is.null(mu), 1, mu ^ 2 / sigma ^ 2)
    b <- ifelse(is.null(sigma), 1, mu / sigma ^ 2)
    pk <- pgamma(q, shape = a, rate = b)
  }
  else if (distr == 5) {
    a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu / sigma) ^ 2 -
                  mu)
    b <-
      ifelse(is.null(sigma), 1 / sqrt(12), (mu * (1 - mu) / sigma ^ 2 -
                                              1) * (1 - mu))
    if (any(c(a, b) <= 0))
      stop(paste("\nNegative Beta parameters:\n a =", a,
                 ";\t b =", b))
    pk <- pbeta(q, shape1 = a, shape2 = b)
  }
  else {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4 or 5")
  }
  return(pk)
}



dkcens2 = function(xleft, xright, c_code_filters, distr = NULL, mu = NULL, sigma = NULL) {
  # c_code = censor_code_rl(xleft, xright)

  res = seq_along(xleft) #initialisation

  res[c_code_filters[['1']]] = dk(x = xleft[c_code_filters[['1']]], distr, mu, sigma)
  res[c_code_filters[['2']]] = pk(xright[c_code_filters[['2']]], distr, mu, sigma)
  res[c_code_filters[['0']]] = 1 - pk(xleft[c_code_filters[['0']]], distr, mu, sigma)
  res[c_code_filters[['3']]] = pk(xright[c_code_filters[['3']]], distr, mu, sigma) -
    pk(xleft[c_code_filters[['3']]], distr, mu, sigma)

  return(res)
}

dkcens2_1val = function(xleft, xright, c_code, distr = NULL, mu = NULL, sigma = NULL) {
  if(c_code == 1) dk(x = xleft, distr, mu, sigma)
  else if(c_code == 2) pk(xright, distr, mu, sigma)
  else if(c_code == 0) 1 - pk(xleft, distr, mu, sigma)
  else if(c_code == 3) pk(xright, distr, mu, sigma) - pk(xleft, distr, mu, sigma)
}

rfystarcens2 = function (v, v2, xleft, xright, censor_code, distr.k, sigma.k, distr.p0, mu.p0, sigma.p0)
{
  alpha <-
    p0(v, distr = distr.p0, mu = mu.p0, sigma = sigma.p0) / p0(v2,
                                                               distr = distr.p0, mu = mu.p0, sigma = sigma.p0)
  Prod <- 1
  for (i in seq_along(xleft)) {
    fac <-
      dkcens2_1val(
        xleft = xleft[i], xright = xright[i], c_code = censor_code[i], distr = distr.k, mu = v, sigma = sigma.k
      ) / dkcens2_1val(
        xleft = xleft[i], xright = xright[i], c_code = censor_code[i],
        distr = distr.k, mu = v2, sigma = sigma.k
      )
    Prod <- Prod * fac
  }
  f <- alpha * Prod
  return(f)
}

gs4cens2 = function (ystar, xleft, xright, censor_code, idx, distr.k, sigma.k, distr.p0, mu.p0, sigma.p0)
{
  r <- length(ystar)
  nstar <- as.numeric(table(idx))
  for (j in seq(r)) {
    id <- which(!is.na(match(idx, j)))
    xjleft <- xleft[id]
    xjright <- xright[id]
    xbar <-
      0.5 * sum(xjleft + xjright, na.rm = T) / nstar[j]####We estimate the cluster mean by choosing centre interval censored data and removing left/right censored data
    y2star <-
      rk(1, distr = distr.k, mu = xbar, sigma = sigma.k / sqrt(nstar[j]))
    f.ratio <-
      rfystarcens2(
        v = y2star, v2 = ystar[j], xleft = xjleft, xright = xjright, censor_code = censor_code[id], distr.k = distr.k,
        sigma.k = sigma.k, distr.p0 = distr.p0, mu.p0 = mu.p0,
        sigma.p0 = sigma.p0
      )
    k.ratio <-
      dk(
        ystar[j], distr = distr.k, mu = xbar, sigma = sigma.k / sqrt(nstar[j])
      ) / dk(
        y2star,
        distr = distr.k, mu = xbar, sigma = sigma.k /
          sqrt(nstar[j])
      )
    q2 <- min(1, f.ratio * k.ratio)
    ystar[j] <- ifelse(runif(1) <= q2, y2star, ystar[j])
  }
  return(ystar)
}

fcondYXAcens2 = function (xleft, xright, censor_code_filters, distr, Tau, J, sigma)
{
  K <- matrix(NA, nrow = length(Tau), ncol = length(xleft))
  for (i in seq(Tau)) {
    K[i,] <-
      dkcens2(
        xleft = xleft, xright = xright, c_code_filters = censor_code_filters, distr = distr, mu = Tau[i], sigma = sigma
      ) *
      J[i]
  }
  pK <- prop.table(K, margin = 2)
  y <- apply(pK, 2, function(x)
    sample(Tau, size = 1, prob = x))
  return(y)
}


gs5cens2 = function (sigma, xleft, xright, censor_code, y, distr = 1, asigma = 1, bsigma = 2, delta = 4)
{
  sigmaStar <- rgamma(1, shape = delta, rate = delta / sigma)
  sigmaT <- sigma
  qgammas <- sigmaT / sigmaStar
  Qgammas <- sigmaStar / sigmaT
  Term2 <- qgammas ^ (2 * delta - 1) * exp(-delta * (qgammas -
                                                       Qgammas))
  Kgamma <- Qgammas ^ (asigma - 1) * exp(-bsigma * (sigmaStar -
                                                      sigmaT))
  Prod <- 1
  for (i in seq_along(xleft)) {
    Prod <-
      Prod * (
        dkcens2_1val(
          xleft = xleft[i], xright = xright[i], c_code = censor_code[i], distr = distr, mu = y[i], sigma = sigmaStar
        ) / dkcens2_1val(
          xleft = xleft[i], xright = xright[i], c_code = censor_code[i],
          distr = distr, mu = y[i], sigma = sigmaT
        )
      )
  }
  q3 <- min(1, Kgamma * Prod * Term2)
  sigma <- ifelse(runif(1) <= q3, sigmaStar, sigmaT)
  return(sigma)
}

rfyzstarcens2 = function (v, v2, z, z2, xleft, xright, censor_code, distr.k, distr.py0, mu.py0, sigma.py0, distr.pz0, mu.pz0, sigma.pz0)
{
  alpha <- p0(v, distr = distr.py0, mu = mu.py0, sigma = sigma.py0)/p0(v2,
                                                                       distr = distr.py0, mu = mu.py0, sigma = sigma.py0) *
    p0(z, distr = distr.pz0, mu = mu.pz0, sigma = sigma.pz0)/p0(z2,
                                                                distr = distr.pz0, mu = mu.pz0, sigma = sigma.pz0)
  Prod <- 1
  for (i in seq_along(xleft)) {
    fac <- dkcens2_1val(xleft = xleft[i], xright = xright[i], c_code = censor_code[i],
                        distr = distr.k, mu = v, sigma = z)/dkcens2_1val(xleft = xleft[i], xright = xright[i], c_code = censor_code[i],
                                                                         distr = distr.k, mu = v2, sigma = z2)
    Prod <- Prod * fac
  }
  f <- alpha * Prod
  return(f)
}


gsYZstartcens2 = function (ystar, zstar, nstar, rstar, idx, xleft, xright, censor_code, delta, kappa, distr.k,
                           distr.py0, mu.py0, sigma.py0, distr.pz0, mu.pz0, sigma.pz0)
{
  for (j in seq(rstar)) {
    flag <- 1
    while (flag == 1) {
      id <- which(!is.na(match(idx, j)))
      xjleft <- xleft[id]
      xjright <- xright[id]
      xbar <-
        0.5 * sum(xjleft + xjright, na.rm = T) / nstar[j]####We estimate the cluster mean by choosing centre interval censored data and removing left/right censored data

      z2star <- rk(1, distr = distr.pz0, mu = zstar[j],
                   sigma = zstar[j]/sqrt(delta))
      y2star <- rk(1, distr = distr.py0, mu = xbar, sigma = kappa *
                     z2star/sqrt(nstar[j]))
      f.ratio <- rfyzstarcens2(v = y2star, v2 = ystar[j], z = z2star, z2 = zstar[j],
                               xleft = xjleft, xright = xjright, censor_code = censor_code[id],
                               distr.k = distr.k, distr.py0 = distr.py0,
                               mu.py0 = mu.py0, sigma.py0 = sigma.py0, distr.pz0 = distr.pz0,
                               mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0)
      k.ratioNum <- dk(zstar[j], distr = distr.pz0, mu = z2star,
                       sigma = z2star/sqrt(delta))
      k.ratioDen <- dk(z2star, distr = distr.pz0, mu = zstar[j],
                       sigma = zstar[j]/sqrt(delta))
      k.ratio <- k.ratioNum/k.ratioDen
      k.ratioNum <- dk(ystar[j], distr = distr.py0, mu = xbar,
                       sigma = kappa * zstar[j]/sqrt(nstar[j]))
      k.ratioDen <- dk(y2star, distr = distr.py0, mu = xbar,
                       sigma = kappa * z2star/sqrt(nstar[j]))
      k.ratio <- k.ratio * k.ratioNum/k.ratioDen
      q2 <- min(1, f.ratio * k.ratio)
      if (is.na(q2)) {
        flag <- 1
      }
      else {
        flag <- 0
        if (runif(1) <= q2) {
          ystar[j] <- y2star
          zstar[j] <- z2star
        }
      }
    }
  }
  return(list(ystar = ystar, zstar = zstar))
}

fcondYZXAcens2 = function (xleft, xright, censor_code_filters, distr, Tauy, Tauz, J)
{
  K <- matrix(NA, nrow = length(Tauy), ncol = length(xleft))
  for (i in seq(Tauy)) {
    K[i, ] <- dkcens2(xleft, xright, c_code_filters = censor_code_filters, distr = distr, mu = Tauy[i], sigma = Tauz[i]) *
      J[i]
  }
  if (any(is.na(K)))
    print(K, Tauy, Tauz, J)
  pK <- prop.table(K, margin = 2)
  j <- apply(pK, 2, function(x) sample(length(Tauy), size = 1,
                                       prob = x))
  return(matrix(c(y = Tauy[j], z = Tauz[j]), nrow = length(xleft),
                ncol = 2))
}

fcondXA2cens2 = function (xleft, xright, censor_code_filters, distr, Tauy, Tauz, J)
{
  pJ <- J/sum(J)
  K <- matrix(NA, nrow = length(Tauy), ncol = length(xleft))
  for (i in seq(Tauy)) {
    K[i, ] <- dkcens2(xleft = xleft, xright = xright, c_code_filters = censor_code_filters, distr = distr, mu = Tauy[i], sigma = Tauz[i])
  }
  fcondXA2cens <- apply(K, 2, function(x) sum(x * pJ))
  return(fcondXA2cens)
}


cens_data_check = function(xleft, xright){
  if(any(xright<xleft, na.rm = T)) stop('in censored data, left bound not always smaller than right bound')
  if(any(mapply(FUN = function(xileft, xiright){is.na(xileft)&is.na(xiright)}, xleft, xright))) stop('in censored data, there is an NA NA')
}


print_progress = function(j, Nit, tInit){
  if (floor(j / 1000) == ceiling(j / 1000))
    cat("MCMC iteration", j, "of", Nit, ' in ', proc.time()[3]-tInit[3]," seconds\n")
}

MixNRMI2denscens_mvcpp_faster =  function (xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1, Beta = 0,
                                             Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2, mu.pz0 = 3,
                                             sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2, Meps = 0.01,
                                             Nx = 100, Nit = 500, Pbi = 0.1, epsilon = NULL, printtime = TRUE)
{
  #
  #   ########for debug
  #   probs = c(0.025, 0.5, 0.975); Alpha = 1; Beta = 0;
  #   Gama = 0.4; distr.k = 1; distr.py0 = 1; distr.pz0 = 2; mu.pz0 = 3;
  #   sigma.pz0 = sqrt(10); delta = 4; kappa = 2; Delta = 2; Meps = 0.01;
  #   Nx = 100; Nit = 500; Pbi = 0.1; epsilon = NULL; printtime = TRUE
  #   #######

  if (is.null(distr.k))
    stop("Argument distr.k is NULL. Should be provided. See help for details.")
  if (is.null(distr.py0))
    stop("Argument distr.py0 is NULL. Should be provided. See help for details.")

  cens_data_check(xleft, xright)

  tInit <- proc.time()
  xpoint = as.numeric(
    na.omit(
      0.5*(xleft+xright)
    )
  )

  censor_code = censor_code_rl(xleft, xright)
  censor_code_filters = lapply(0:3, FUN = function(x) censor_code==x)
  names(censor_code_filters) = 0:3


  n <- length(xleft)
  y <- seq(n)
  y[seq(n / 2)] <- mean(xpoint[1:(length(xpoint) / 2)])
  y[-seq(n / 2)] <- mean(xpoint[(length(xpoint) / 2):length(xpoint)])
  z <- rep(1, n)
  u <- 1
  if (is.null(epsilon))
    epsilon <- sd(xpoint) / 4
  xx <- seq(min(xpoint) - epsilon, max(xpoint) + epsilon, length = Nx)
  Fxx <- matrix(NA, nrow = Nx, ncol = Nit)
  fx <- matrix(NA, nrow = n, ncol = Nit)
  R <- seq(Nit)
  U <- seq(Nit)
  Nmt <- seq(Nit)

  means = vector(mode = "list", length = Nit)
  sigmas = vector(mode = "list", length = Nit)
  weights = vector(mode = "list", length = Nit)

  mu.py0 = mean(xpoint)
  sigma.py0 = sd(xpoint)

  for (j in seq(Nit)) {
    print_progress(j, Nit, tInit)
    tt <- comp2(y, z)
    ystar <- tt$ystar
    zstar <- tt$zstar
    nstar <- tt$nstar
    rstar <- tt$rstar
    idx <- tt$idx
    if (Gama != 0)
      u <- gs3(
        u, n = n, r = rstar, alpha = Alpha, beta = Beta,
        gama = Gama, delta = Delta
      )
    JiC <- MvInv_cpp(
      eps = Meps, u = u, alpha = Alpha, beta = Beta,
      gama = Gama, N = 50001
    )
    Nm <- length(JiC)
    TauyC <-
      rk(Nm, distr = distr.py0, mu = mu.py0, sigma = sigma.py0)
    TauzC <-
      rk(Nm, distr = distr.pz0, mu = mu.pz0, sigma = sigma.pz0)
    tt <- gsYZstartcens2(ystar = ystar, zstar = zstar, nstar = nstar, rstar = rstar, idx = idx,
                         xleft = xleft, xright = xright, censor_code = censor_code, delta = delta,
                         kappa = kappa, distr.k = distr.k, distr.py0 = distr.py0,
                         mu.py0 = mu.py0, sigma.py0 = sigma.py0, distr.pz0 = distr.pz0,
                         mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0
    )
    ystar <- tt$ystar
    zstar <- tt$zstar
    tt <- gsHP(ystar, rstar, distr.py0)
    mu.py0 <- tt$mu.py0
    sigma.py0 <- tt$sigma.py0
    Jstar <- rgamma(rstar, nstar - Gama, Beta + u)
    Tauy <- c(TauyC, ystar)
    Tauz <- c(TauzC, zstar)
    J <- c(JiC, Jstar)
    tt <- fcondYZXAcens2(xleft = xleft, xright = xright, censor_code_filters = censor_code_filters,
                         distr = distr.k, Tauy = Tauy, Tauz = Tauz, J = J)
    y <- tt[, 1]
    z <- tt[, 2]
    Fxx[, j] <- fcondXA2(xx, distr = distr.k, Tauy, Tauz,
                         J)
    fx[, j] <- fcondXA2cens2(xleft = xleft, xright = xright, censor_code_filters = censor_code_filters, distr = distr.k, Tauy, Tauz, J)
    R[j] <- rstar
    U[j] <- u
    Nmt[j] <- Nm
    means[[j]] = Tauy
    sigmas[[j]] = Tauz
    weights[[j]] = J / sum(J)
  }
  biseq <- seq(floor(Pbi * Nit))
  Fxx <- Fxx[,-biseq]
  qx <- as.data.frame(t(apply(Fxx, 1, quantile, probs = probs)))
  names(qx) <- paste("q", probs, sep = "")
  qx <- cbind(mean = apply(Fxx, 1, mean), qx)
  R <- R[-biseq]
  U <- U[-biseq]
  means <- means[-biseq]
  sigmas <- sigmas[-biseq]
  weights <- weights[-biseq]
  cpo <- 1 / apply(1 / fx[,-biseq], 1, mean)
  if (printtime) {
    cat(" >>> Total processing time (sec.):\n")
    print(procTime <- proc.time() - tInit)
  }
  return(
    list(
      xx = xx, qx = qx, cpo = cpo, R = R, U = U, Nm = Nmt,
      Nx = Nx, Nit = Nit, Pbi = Pbi, procTime = procTime,
      means = means, sigmas = sigmas, weights = weights
    )
  )
}

MixNRMI1denscens_mvcpp_faster = function (xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1, Beta = 0,
                                            Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5, bsigma = 0.5,
                                            delta = 3, Delta = 2, Meps = 0.01, Nx = 100, Nit = 500, Pbi = 0.1,
                                            epsilon = NULL, printtime = TRUE){


  #   #for debug
  #     probs = c(0.025, 0.5, 0.975); Alpha = 1; Beta = 0;
  #     Gama = 0.4; distr.k = 1; distr.p0 = 1; asigma = 0.5; bsigma = 0.5;
  #     delta = 3; Delta = 2; Meps = 0.01; Nx = 100; Nit = 500; Pbi = 0.1;
  #     epsilon = NULL; printtime = TRUE
  #   ##

  if (is.null(distr.k))
    stop("Argument distr.k is NULL. Should be provided. See help for details.")
  if (is.null(distr.p0))
    stop("Argument distr.p0 is NULL. Should be provided. See help for details.")

  cens_data_check(xleft, xright)


  tInit <- proc.time()
  xpoint = as.numeric(
    na.omit(
      0.5*(xleft+xright)
    )
  )

  censor_code = censor_code_rl(xleft, xright)
  censor_code_filters = lapply(0:3, FUN = function(x) censor_code==x)
  names(censor_code_filters) = 0:3

  n <- length(xleft)
  y <- seq(n)
  y[seq(n / 2)] <- mean(xpoint[1:(length(xpoint) / 2)])
  y[-seq(n / 2)] <- mean(xpoint[(length(xpoint) / 2):length(xpoint)])
  u <- 1
  sigma <- 1
  if (is.null(epsilon))
    epsilon <- sd(xpoint) / 4
  xx <- seq(min(xpoint) - epsilon, max(xpoint) + epsilon, length = Nx)
  Fxx <- matrix(NA, nrow = Nx, ncol = Nit)
  fx <- matrix(NA, nrow = n, ncol = Nit)
  R <- seq(Nit)
  S <- seq(Nit)
  U <- seq(Nit)
  Nmt <- seq(Nit)
  means = vector(mode = "list", length = Nit)
  sigmas = vector(mode = "list", length = Nit)
  weights = vector(mode = "list", length = Nit)
  mu.p0 = mean(xpoint)
  sigma.p0 = sd(xpoint)
  for (j in seq(Nit)) {
    print_progress(j, Nit, tInit)

    tt <- comp1(y)
    ystar <- tt$ystar
    nstar <- tt$nstar
    r <- tt$r
    idx <- tt$idx
    if (Gama != 0)
      u <- gs3(
        u, n = n, r = r, alpha = Alpha, beta = Beta,
        gama = Gama, delta = Delta
      )
    JiC <- MvInv_cpp(
      eps = Meps, u = u, alpha = Alpha, beta = Beta,
      gama = Gama, N = 50001
    )
    Nm <- length(JiC)
    TauiC <- rk(Nm, distr = distr.p0, mu = mu.p0, sigma = sigma.p0)
    ystar <- gs4cens2(ystar = ystar, xleft = xleft, xright = xright, censor_code = censor_code,
                      idx =  idx, distr.k = distr.k, sigma.k = sigma,
                      distr.p0 = distr.p0, mu.p0 = mu.p0, sigma.p0 = sigma.p0
    )
    Jstar <- rgamma(r, nstar - Gama, Beta + u)
    Tau <- c(TauiC, ystar)
    J <- c(JiC, Jstar)
    tt <- gsHP(ystar, r, distr.p0)
    mu.p0 <- tt$mu.py0
    sigma.p0 <- tt$sigma.py0
    y <-
      fcondYXAcens2(
        xleft = xleft, xright = xright, censor_code_filters = censor_code_filters,
        distr = distr.k, Tau = Tau, J = J, sigma = sigma)
    sigma <- gs5cens2(
      sigma = sigma, xleft = xleft, xright = xright, censor_code = censor_code,
      y = y, distr = distr.k, asigma = asigma,
      bsigma = bsigma, delta = delta
    )
    Fxx[, j] <- fcondXA(
      xx, distr = distr.k, Tau = Tau, J = J,
      sigma = sigma
    )
    fx[, j] <- fcondYXAcens2(
      xleft = xleft, xright = xright, censor_code_filters = censor_code_filters, distr = distr.k, Tau = Tau, J = J, sigma = sigma
    )

    R[j] <- r
    S[j] <- sigma
    U[j] <- u
    Nmt[j] <- Nm
    means[[j]] = Tau
    sigmas[[j]] = sigma
    weights[[j]] = J / sum(J)
  }
  biseq <- seq(floor(Pbi * Nit))
  Fxx <- Fxx[,-biseq]
  qx <- as.data.frame(t(apply(Fxx, 1, quantile, probs = probs)))
  names(qx) <- paste("q", probs, sep = "")
  qx <- cbind(mean = apply(Fxx, 1, mean), qx)
  R <- R[-biseq]
  S <- S[-biseq]
  U <- U[-biseq]
  means <- means[-biseq]
  sigmas <- sigmas[-biseq]
  weights <- weights[-biseq]
  cpo <- 1 / apply(1 / fx[,-biseq], 1, mean)
  if (printtime) {
    cat(" >>> Total processing time (sec.):\n")
    print(procTime <- proc.time() - tInit)
  }
  return(
    list(
      xx = xx, qx = qx, cpo = cpo, R = R, S = S, U = U,
      Nm = Nmt, Nx = Nx, Nit = Nit, Pbi = Pbi, procTime = procTime,
      means = means, sigmas = sigmas, weights = weights
    )
  )
}

environment(MixNRMI1denscens_mvcpp_faster) <- environment(MixNRMI1)

environment(pk) <- environment(MixNRMI1)

environment(dkcens2_1val) <- environment(MixNRMI1)
environment(dkcens2) <- environment(MixNRMI1)
environment(rfystarcens2) <- environment(MixNRMI1)
environment(fcondYXAcens2) <- environment(MixNRMI1)
environment(gs5cens2) <- environment(MixNRMI1)
environment(gs4cens2) <- environment(MixNRMI1)
environment(rfyzstarcens2) <- environment(MixNRMI2)
environment(fcondYZXAcens2) <- environment(MixNRMI2)
environment(fcondXA2cens2) <- environment(MixNRMI2)
environment(gsYZstartcens2) <- environment(MixNRMI2)

environment(MixNRMI2denscens_mvcpp_faster) <- environment(MixNRMI2)
