# need to make this an S4 class (npsp::grid.par)
setClass("grid.par")


# class representing a distribution p(n, x, theta) on a grid.
# The idea is to completely abstract away the respective densities and only use
# these pivots for later evaluation of the predictive/marginal distribution and
# the poserior.
#
# We need p(theta | x, n), i.e., the posterior and p(x | n), i.e., the marginal/
# predictive distribution. Both are only required conditional on n.
#
# Here, (x, n) is the stage-wise sufficient statistic for theta
#' @export
setClass("StageWiseDataModel", representation(
    marginal_pdfs  = "matrix",     # #n-by-#x matrix
    marginal_cdfs  = "matrix",     # #n-by-#x matrix
    marginal_grid  = "grid.par",
    posterior_pdfs = "array",      # #n-by-#x-by-#theta array
    posterior_cdfs = "array",      # #n-by-#x-by-#theta array
    posterior_grid = "grid.par",
    x              = "numeric",
    n              = "numeric",
    theta          = "numeric")
)

#' @export
setGeneric("StageWiseDataModel", function(dist, prior, ...) standardGeneric("StageWiseDataModel"))

# constructor: precompute pdf/cdf for marginal and posterior and store as arrays
#' @export
setMethod("StageWiseDataModel", signature("DataDistribution", "Prior"),
          function(dist, prior,
                   support_n = c(1, 1000), support_x = c(-5, 5),
                   dims = c(NA_real_, NA_real_, NA_real_)) {

              n_max  <- max(support_n)
              n_min  <- min(support_n)
              if (is.na(dims[1])) dims[1] <- ceiling((n_max - n_min + 1) / 5)
              n      <- seq(n_min, n_max, length.out = dims[1])
              dn     <- n[2] - n[1]


              if(is(prior, "ContinuousPrior")) {
                  x_min  <- sqrt(n_max) * prior@support[1] + support_x[1]
                  x_max  <- sqrt(n_max) * prior@support[2] + support_x[2]
              } else{
                  x_min  <- sqrt(n_max) * min(prior@theta) + support_x[1]
                  x_max  <- sqrt(n_max) * max(prior@theta) + support_x[2]
              }
              if (is.na(dims[2])) dims[2] <- ceiling(3 * (x_max - x_min))
              x      <- seq(x_min, x_max, length.out = dims[2])
              dx     <- x[2] - x[1]


              if(is(prior, "ContinuousPrior")) {
                  theta_min <- prior@support[1]
                  theta_max <- prior@support[2]
                  if (is.na(dims[3])) dims[3] <- 50 * ceiling(theta_max - theta_min)
                  theta  <- seq(theta_min, theta_max, length.out = dims[3])
                  dtheta <- theta[2] - theta[1]

                  prior_pdf <- function(th) prior@pdf(th)

              } else {
                  theta   <- prior@theta
                  mass    <- prior@mass
                  if(length(theta) == 1) {
                      theta <- c(prior@theta - .Machine$double.eps, prior@theta + .Machine$double.eps)
                      mass  <- rep(prior@mass, 2)
                  }

                  dims[3] <- length(theta)
                  dtheta  <- 1 #/ length(theta) # todo: check

                  prior_pdf <- function(th) {
                      res <- rep(0, length(th))
                      for(i in 1:length(theta)) {
                          for(j in 1:length(th)) if(th[j] == theta[i]) res[j] <- mass[i]
                      }
                      return(res)
                  }
              }


                pdfs <- expand.grid(
                    n = n,
                    x = x,
                    theta = theta
                ) %>%
                mutate(
                    joint_pdf = log(probability_density_function(dist, x, n, theta)) + log(prior_pdf(theta))
                ) %>%
                group_by(n) %>% # normalize
                mutate(
                    joint_pdf = joint_pdf - log(sum(exp(joint_pdf))) - log(dx * dtheta)
                ) %>%
                ungroup() %>%
                group_by(n, x) %>%
                mutate(
                    marginal_pdf  = exp(log(sum(exp(joint_pdf))) + log(dtheta)),
                    posterior_pdf = ifelse(marginal_pdf == 0, 0, exp(joint_pdf) / marginal_pdf)
                ) %>%
                ungroup()

                marginal_grid <- npsp::grid.par(dims[1:2],
                                                min = c(min(n), min(x)),
                                                max = c(max(n), max(x)))

                marginal_pdfs <- pdfs %>%
                    distinct(n, x, marginal_pdf) %>%
                    pull(marginal_pdf) %>%
                    matrix(nrow = length(n), byrow = FALSE) # n-by-x matrix of marginals

                marginal_pdfs %>%
                    apply(., 1, function(y) cumsum(y) * dx) %>%
                    t(.) -> marginal_cdfs

                posterior_grid <- npsp::grid.par(dims,
                                                 min = c(min(n), min(x), min(theta)),
                                                 max = c(max(n), max(x), max(theta)))

                posterior_pdfs <- pdfs %>%
                    pull(posterior_pdf) %>%
                    array(dim = dims) # n-by-x-by-theta array of posteriors

                if(is(prior, "PointMassPrior") && length(prior@theta) == 1)  posterior_pdfs <- 2 * posterior_pdfs

                posterior_cdfs <- posterior_pdfs
                for (i in 1:length(n)) {
                    for (j in 1:length(x)) {
                        posterior_cdfs[i, j, ] <- cumsum(posterior_pdfs[i, j, ]) * dtheta
                    }
                }

                new("StageWiseDataModel", marginal_pdfs = marginal_pdfs, marginal_cdfs = marginal_cdfs,
                    marginal_grid = marginal_grid, posterior_pdfs = posterior_pdfs,
                    posterior_cdfs = posterior_cdfs, posterior_grid = posterior_grid,
                    n = n, x = x, theta = theta)

})


#' @export
marginal_pdf <- function(m, n, x) { # vectorized over x, not n
    npsp::interp(m@marginal_grid, m@marginal_pdfs, cbind(n, x))$y
}

#' @export
posterior_pdf <- function(m, theta, n, x) { # vectorized over theta, not n, x
    npsp::interp(m@posterior_grid, m@posterior_pdfs, cbind(n, x, theta))$y
}
