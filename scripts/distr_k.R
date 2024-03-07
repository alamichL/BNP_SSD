source('scripts/halfcauchy.R')
source('scripts/halfnorm.R')
source('scripts/halfstudent.R')
library(msm)

rk = function (n, distr = NULL, mu = NULL, sigma = NULL)
{
  if (is.null(distr)) {
    stop("Argument \"distr\" should be defined numeric with possible values 1 (normal), 2 (gamma), 3 (beta), 4 (exponential), 5 (lognormal), 6 (half-Cauchy), 7 (half-normal), 8 (half-student) to 8")
  }
  else if (distr == 1) {
    rk <- rnorm(n, mean = mu, sd = sigma)
  }
  else if (distr == 4) {
    a <- ifelse(is.null(mu), 0, mu)
    b <- ifelse(is.null(sigma), 1/sqrt(2), sigma/sqrt(2))
    rk <- a + b * sample(c(-1, +1), size = n, replace = TRUE) *
      rexp(n)
  }
  else if (distr == 5) {
    a <- ifelse(is.null(mu), exp(1/2), log(mu/sqrt(1 + (sigma/mu)^2)))
    b <- ifelse(is.null(sigma), exp(1) * (exp(1) - 1), sqrt(log(1 +
                                                                  (sigma/mu)^2)))
    rk <- rlnorm(n, meanlog = a, sdlog = b)
  }
  else if (distr == 2) {
    a <- ifelse(is.null(mu), 1, mu^2/sigma^2)
    b <- ifelse(is.null(sigma), 1, mu/sigma^2)
    rk <- rgamma(n, shape = a, rate = b)
  }
  else if (distr == 3) {
    a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu/sigma)^2 -
                  mu)
    b <- ifelse(is.null(sigma), 1/sqrt(12), (mu * (1 - mu)/sigma^2 -
                                               1) * (1 - mu))
    if (any(c(a, b) <= 0))
      stop(paste("\nNegative Beta parameters:\n a =", a,
                 ";\t b =", b))
    rk <- rbeta(n, shape1 = a, shape2 = b)
  }
  else if (distr == 6) {
    rk <- rhalfcauchy(n, location =  ifelse(is.null(mu), 0, mu),
                      scale = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 7) {
    rk <- rhalfnorm(n, mean =  ifelse(is.null(mu), 0, mu),
                      sd = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 8) {
    rk <- rhalft(n, df = 10, mean =  ifelse(is.null(mu), 0, mu),
                 sd = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 9) {
    rk <- runif(n, min =  ifelse(is.null(mu), 0, mu),
                 max = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 10) {
    rk <- rtnorm(n, mean = ifelse(is.null(mu), 0, mu),
                 sd = ifelse(is.null(sigma), 1, sigma),
                 lower = 0.1)
  }
  else {
    stop("Argument \"distr\" should be defined numeric with possible values 1 to 8")
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
    b <- ifelse(is.null(sigma), 1/sqrt(2), sigma/sqrt(2))
    dk <- exp(-abs(x - a)/b)/(2 * b)
  }
  else if (distr == 5) {
    a <- ifelse(is.null(mu), exp(1/2), log(mu/sqrt(1 + (sigma/mu)^2)))
    b <- ifelse(is.null(sigma), exp(1) * (exp(1) - 1), sqrt(log(1 +
                                                                  (sigma/mu)^2)))
    dk <- dlnorm(x, meanlog = a, sdlog = b)
  }
  else if (distr == 2) {
    a <- ifelse(is.null(mu), 1, mu^2/sigma^2)
    b <- ifelse(is.null(sigma), 1, mu/sigma^2)
    dk <- dgamma(x, shape = a, rate = b)
  }
  else if (distr == 3) {
    a <- ifelse(is.null(mu), 0.5, (1 - mu) * (mu/sigma)^2 -
                  mu)
    b <- ifelse(is.null(sigma), 1/sqrt(12), (mu * (1 - mu)/sigma^2 -
                                               1) * (1 - mu))
    if (any(c(a, b) <= 0))
      stop(paste("\nNegative Beta parameters:\n a =", a,
                 ";\t b =", b))
    dk <- dbeta(x, shape1 = a, shape2 = b)
  }
  else if (distr == 6) {
    dk <- dhalfcauchy(x, location =  ifelse(is.null(mu), 0, mu),
                      scale = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 7) {
    dk <- dhalfnorm(x, mean =  ifelse(is.null(mu), 0, mu),
                      sd = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 8) {
    dk <- dhalft(x, df = 10, mean =  ifelse(is.null(mu), 0, mu),
                      sd = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 9) {
    dk <- dunif(x, min =  ifelse(is.null(mu), 0, mu),
                max = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 10) {
    dk <- dtnorm(x, mean = ifelse(is.null(mu), 0, mu),
                 sd = ifelse(is.null(sigma), 1, sigma),
                 lower = 0.1)
  }
  else {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4,5 or 6")
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
  else if (distr == 6) {
    pk <- phalfcauchy(x, location =  ifelse(is.null(mu), 0, mu),
                      scale = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 7) {
    pk <- phalfnorm(x, mean =  ifelse(is.null(mu), 0, mu),
                      sd = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 8) {
    pk <- phalft(x, df = 10, mean =  ifelse(is.null(mu), 0, mu),
                 sd = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 9) {
    pk <- punif(x, min =  ifelse(is.null(mu), 0, mu),
                 max = ifelse(is.null(sigma), 1, sigma))
  }
  else if (distr == 10) {
    pk <- ptnorm(q, mean = ifelse(is.null(mu), 0, mu),
                 sd = ifelse(is.null(sigma), 1, sigma),
                 lower = 0.1)
  }
  else {
    stop("Argument \"distr\" should be defined numeric with possible values 1,2,3,4,5 or 6")
  }
  return(pk)
}

