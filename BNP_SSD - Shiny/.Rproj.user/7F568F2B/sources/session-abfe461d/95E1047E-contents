#' Compute the optimal clustering from an MCMC sample
#'
#' Summarizes the posterior on all possible clusterings by an optimal
#' clustering where optimality is defined as minimizing the posterior
#' expectation of a specific loss function, the Variation of Information or
#' Binder's loss function. Computation can be lengthy for large datasets,
#' because of the large size of the space of all clusterings.
#'
#'
#' @param fit The fitted object, obtained from one of the MixNRMIx functions
#' @param loss_type Defines the loss function to be used in the expected
#' posterior loss minimization. Can be one of "VI" (Variation of Information),
#' "B" (Binder's loss), "NVI" (Normalized Variation of Information) or "NID"
#' (Normalized Information Distance). Defaults to "VI".
#' @return A vector of integers with the same size as the data, indicating the
#' allocation of each data point.
#' @export compute_optimal_clustering
compute_optimal_clustering <- function(fit, loss_type = "VI") {
  if (!requireNamespace("GreedyEPL", quietly = TRUE)) {
    stop("Package GreedyEPL is needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  fit.draw <- Reduce(rbind, fit$Allocs)
  fit_VI <- GreedyEPL::MinimiseEPL(sample_of_partitions = fit.draw, pars = list("loss_type" = loss_type))
  return(fit_VI$decision)
}

# clustering = compute_optimal_clustering(out)

plot_clustering_and_CDF_noncensored <- function(fit, clustering, label_vector = NULL) {
  data <- fit$data

  grid <- grid_from_data(data)

  if (is_semiparametric(fit)) {
    cdf <- get_CDF_semi_BNPdensity(fit = fit, xs = data)
  }
  else {
    cdf <- get_CDF_full_BNPdensity(fit = fit, xs = data)
  }
  p <- ggplot(data.frame(data = data, cluster_id = clustering, cdf = cdf)) +
    theme_bw() +
    geom_step(aes_string(x = "data", y = "ecdf(data)(data)")) +
    geom_point(aes_string(x = "data", y = "cdf", colour = "factor(cluster_id)")) +
    viridis::scale_colour_viridis(discrete = TRUE) +
    theme(legend.position = "none") +
    ylab("CDF") +
    xlab("Data")

  if (!is.null(label_vector)) {
    p + geom_text(
      data = data.frame(
        txt = label_vector, x = data,
        y = cdf + 0.05,
        cluster_id = clustering
      ),
      aes_string(x = "x", y = "y", colour = "factor(cluster_id)", label = "txt")
    )
  }
  else {
    return(p)
  }
}

decide_abscissa <- function(censored_data, clustering) {
  df <- cbind(censored_data, data.frame(
    cluster_id = clustering,
    loc = rowMeans(censored_data),
    is_censored = is.na(rowMeans(censored_data))
  ))
  df$loc <- unlist(mapply(FUN = function(is_cens, cluster_id, loc) {
    if (is_cens) {
      mean(df$loc[df$cluster_id == cluster_id], na.rm = TRUE)
    }
    else {
      loc
    }
  }, df$is_censored, df$cluster_id, df$loc, SIMPLIFY = FALSE))
  return(df)
}

plot_clustering_and_CDF_censored <- function(fit, clustering, label_vector = NULL) {
  data <- fit$data

  # grid <- grid_from_data(data)
  grid <- decide_abscissa(data, clustering)$loc

  Survival_object <- survival::survfit(formula = survival::Surv(data$left, data$right, type = "interval2") ~ 1)

  if (is_semiparametric(fit)) {
    cdf <- get_CDF_semi_BNPdensity(fit = fit, xs = grid[!is.na(grid)])
  }
  else {
    cdf <- get_CDF_full_BNPdensity(fit = fit, xs = grid[!is.na(grid)])
  }
  p <- ggplot2::ggplot(
    data = data.frame(data = grid[!is.na(grid)], CDF = cdf, cluster_id = clustering[!is.na(grid)]),
    aes_string(x = "data", y = "CDF")
  ) +
    geom_point(aes_string(colour = "factor(cluster_id)")) +
    theme_classic() +
    geom_step(
      data = data.frame(
        x = c(Survival_object$time, max(grid)),
        y = c(1 - Survival_object$surv, 1)
      ),
      aes_string(x = "x", y = "y")
    ) +
    viridis::scale_colour_viridis(discrete = TRUE) +
    theme(legend.position = "none") +
    ylab("CDF") +
    xlab("Data")

  if (!is.null(label_vector)) {
    p + geom_text(
      data = data.frame(
        txt = label_vector[!is.na(grid)],
        x = grid[!is.na(grid)],
        y = cdf + 0.05,
        cluster_id = clustering[!is.na(grid)]
      ),
      aes_string(x = "x", y = "y", colour = "factor(cluster_id)", label = "txt")
    )
  }
  else {
    return(p)
  }
}



#' Plot the clustering and the Cumulative Distribution Function
#'
#' This is a function to visualize the clustering induced by the BNP model. The
#' data points are plotted with a color reflecting their cluster.
#'
#'
#' @param fit The fitted object, obtained from one of the MixNRMIx functions
#' @param clustering A vector of integers with the same length as the data,
#' representing the allocation variable for data each point.
#' @param label_vector A vector of data labels to be plotted, to provide some
#' identification to each point.
#' @return A plot of the Cumulative Distribution Function (or Turnbull estimate
#' for censored data) with data points whose color denotes the cluster
#' allocation. For censored data, right or left censored data points are not
#' represented, while interval censored data points are represented at the
#' middle of the censoring interval.
#' @export plot_clustering_and_CDF
plot_clustering_and_CDF <- function(fit, clustering, label_vector = NULL) {
  if (is_censored(fit$data)) {
    plot_clustering_and_CDF_censored(fit, clustering, label_vector = label_vector)
  }
  else {
    plot_clustering_and_CDF_noncensored(fit, clustering, label_vector = label_vector)
  }
}
