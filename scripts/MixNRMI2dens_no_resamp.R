source('scripts/MvInv_cpp.R')
source('scripts/distr_k.R')

MixNRMI2dens_no_resamp = function (x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Beta = 0,
                         Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2, mu.pz0 = 3,
                         sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2, Meps = 0.01,
                         Nx = 100, Nit = 500, Pbi = 0.1, epsilon = NULL, printtime = TRUE)
{

#   rk : distr = 1 : normal
#   rk : distr = 2 : gamma
#   rk : distr = 3 : beta
#   rk : distr = 4 : exp
#   rk : distr = 5 : lnorm
#   rk : distr = 6 : half-cauchy

  ###########For debug

#   probs = c(0.025, 0.5, 0.975); Alpha = 1; Beta = 0;
#   Gama = 0.4; distr.k = 1; distr.py0 = 1; distr.pz0 = 2; mu.pz0 = 1;
#   sigma.pz0 = 1; delta = 4; kappa = 2; Delta = 2; Meps = 0.01;
#   Nx = 100; Nit = 500; Pbi = 0.1; epsilon = NULL; printtime = TRUE
#   lapply(list('comp2', 'gs3', 'gsYZstar', 'gsHP', 'p0', 'rk', 'dk', 'fcondYZXA', 'fcondXA2'),
#          function(funname) assign(x = funname, value = getAnywhere(funname)$objs[[1]], pos = .GlobalEnv))
  ##########


  if (is.null(distr.k))
    stop("Argument distr.k is NULL. Should be provided. See help for details.")
  if (is.null(distr.py0))
    stop("Argument distr.py0 is NULL. Should be provided. See help for details.")
  tInit <- proc.time()
  n <- length(x)
  y <- x

  #Initialisation of the means of the clusters y, first half to the mean of half the dataset, second half to the mean of the rest
  #note : kind of strange, since the dataset is not sorted.
  #Plus the last value of y remains equal to the last value of x
  for (i in 1:(n / 2)) {
    y[i] <- mean(x[1:(n / 2)])
  }
  for (i in (n / 2):n) {
    y[i] <- mean(x[(n / 2):n])
  }

  z <- rep(1, n) #set all variance clusters to 1
  u <- 1 #Should be the latent variable

  #Epsilon is used to define the range on which the predictive density will be computed.
  #The grid goes over the edges of the dataset by epsilon
  if (is.null(epsilon))
    epsilon <- sd(x) / 4

  xx <- seq(min(x) - epsilon, max(x) + epsilon, length = Nx) #Grid on which the predictive density will be computed
  Fxx <- matrix(NA, nrow = Nx, ncol = Nit) #Predictive density on the grid
  fx <- matrix(NA, nrow = n, ncol = Nit) #Predictive density on the data to compute CPOs

  R <- seq(Nit) #Number of components
  U <- seq(Nit) #Latent variable
  Nmt <- seq(Nit)

  means = vector(mode = "list", length = Nit)
  sigmas = vector(mode = "list", length = Nit)
  weights = vector(mode = "list", length = Nit)
  Js = vector(mode = "list", length = Nit) #jumps


  #Location and scale of the prior on the means
  mu.py0 = mean(x)
  sigma.py0 = sd(x)

  for (j in seq(Nit)) {
  # for (j in 1:1000) {
    if (floor(j / 100) == ceiling(j / 100))
      cat("MCMC iteration", j, "of", Nit, "in", proc.time()[3] - tInit[3],"s \n")

    # This function computes the distinct observations (couples)
    # and their frequencies in a bivariate numeric vector. see ?comp2
    tt <- comp2(y, z)

    ystar <- tt$ystar
    zstar <- tt$zstar
    nstar <- tt$nstar
    rstar <- tt$rstar
    idx <- tt$idx #allocation variable for each data point

    if (Gama != 0)#Gama = 0 is the Dirichlet process case
      u <- gs3( #This function simulates from the conditional posterior distribution of the latent U.
        ut = u, n = n, r = rstar, alpha = Alpha, beta = Beta,
        gama = Gama, delta = Delta
      )

    JiC <- MvInv_cpp(#Determines the jump heights of an increasing additive process by inverting the M(v) function.
      eps = Meps, u = u, alpha = Alpha, beta = Beta,
      gama = Gama, N = 50001
    )

    Nm <- length(JiC) #Number of components
    TauyC <- #random generation of means
      rk(Nm, distr = distr.py0, mu = mu.py0, sigma = sigma.py0)
    TauzC <- #random generation of sigmas
      rk(Nm, distr = distr.pz0, mu = mu.pz0, sigma = sigma.pz0)

#     #Jointly resampling Ystar and Zstar function
#     tt <- gsYZstar(ystar = ystar, zstar = zstar, nstar = nstar, rstar = rstar, idx = idx, x = x, delta = delta,
#       kappa = kappa, distr.k = distr.k, distr.py0 = distr.py0,
#       mu.py0 = mu.py0, sigma.py0 = sigma.py0, distr.pz0 = distr.pz0,
#       mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0
#     )
#
#     ystar <- tt$ystar
#     zstar <- tt$zstar

    #Updates the hyper-parameters of py0
    tt <- gsHP(ystar, rstar, distr.py0)
    mu.py0 <- tt$mu.py0
    sigma.py0 <- tt$sigma.py0


    Jstar <- rgamma(rstar, nstar - Gama, Beta + u)
    Tauy <- c(TauyC, ystar)
    Tauz <- c(TauzC, zstar)
    J <- c(JiC, Jstar)

    #Conditional posterior distribution of the bivariate latents (Y,Z)
    tt <- fcondYZXA(x = x, distr = distr.k, Tauy = Tauy, Tauz = Tauz, J = J)
    y <- tt[, 1]
    z <- tt[, 2]

    #Conditional density evaluation in the fully nonparametric model
    Fxx[, j] <- fcondXA2(xx, distr = distr.k, Tauy, Tauz,
                         J)
    fx[, j] <- fcondXA2(x, distr = distr.k, Tauy, Tauz, J)

    #save parameters of iteration
    R[j] <- rstar
    U[j] <- u
    Nmt[j] <- Nm
    means[[j]] = Tauy
    sigmas[[j]] = Tauz
    weights[[j]] = J / sum(J)
    Js[[j]] = J

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
  Js <- Js[-biseq]
  cpo <- 1 / apply(1 / fx[,-biseq], 1, mean)
  if (printtime) {
    cat(" >>> Total processing time (sec.):\n")
    print(procTime <- proc.time() - tInit)
  }
  return(
    list(
      xx = xx, qx = qx, cpo = cpo, R = R, U = U, Nm = Nmt,
      Nx = Nx, Nit = Nit, Pbi = Pbi, procTime = procTime, Js = Js,
      means = means, sigmas = sigmas, weights = weights
    )
  )
}


environment(MixNRMI2dens_no_resamp) <- environment(MixNRMI2)
unlockBinding(sym = 'rk', environment(MixNRMI2dens_no_resamp))
assign(x = 'rk', value = rk, pos = environment(MixNRMI2dens_no_resamp))
# unlockBinding(sym = 'dk', environment(MixNRMI2dens_no_resamp)) #unnecessary without the resampling step
# assign(x = 'dk', value = dk, pos = environment(MixNRMI2dens_no_resamp))
