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
    theta          = "numeric"))


# constructor: precompute pdf/cdf for marginal and posterior and store as arrays
#' @export
StageWiseDataModel <- function(data_pdf, prior_pdf, dims = c(100, 100, 100)) {
    n_max   <- attr(data_pdf, "support_n")[2]
    n      <- seq(attr(data_pdf, "support_n")[1], n_max, length.out = dims[1])
    min_x  <- sqrt(n_max) * attr(prior_pdf, "support")[1] + attr(data_pdf, "rel_support_x")[1]
    max_x  <- sqrt(n_max) * attr(prior_pdf, "support")[2] + attr(data_pdf, "rel_support_x")[2]
    x      <- seq(min_x,max_x, length.out = dims[2])
    theta  <- seq(attr(prior_pdf, "support")[1], attr(prior_pdf, "support")[2], length.out = dims[3])
    dn     <- n[2] - n[1]
    dx     <- x[2] - x[1]
    dtheta <- theta[2] - theta[1]

    pdfs <- expand.grid(n = n, x = x, theta = theta) %>%
        mutate(
            joint_pdf = log(data_pdf(n, x, theta)) + log(prior_pdf(theta))
        ) %>%
        group_by(n) %>% # normalize
        mutate(
            joint_pdf = joint_pdf - log(sum(exp(joint_pdf))) - log(dx * dtheta)
        ) %>%
        ungroup() %>%
        group_by(n, x) %>%
        mutate(
            marginal_pdf  = exp(log(sum(exp(joint_pdf))) + log(dtheta)),
            posterior_pdf = exp(joint_pdf) / ifelse(marginal_pdf == 0, 1, marginal_pdf)
        ) %>%
        ungroup()


    marginal_grid <- npsp::grid.par(dims[1:2], c(min(n), min(x)), c(max(n), max(x)))

    marginal_pdfs <- pdfs %>%
        distinct(n, x, marginal_pdf) %>%
        pull(marginal_pdf) %>%
        matrix(nrow = length(n), byrow = FALSE) # n-by-x matrix of marginals

    marginal_pdfs %>%
        apply(., 1, function(y) cumsum(y) * dx) %>%
        t(.) -> marginal_cdfs

    posterior_grid <- npsp::grid.par(dims, c(min(n), min(x), min(theta)), c(max(n), max(x), max(theta)))

    posterior_pdfs <- pdfs %>%
        pull(posterior_pdf) %>%
        array(dim = dims) # n-by-x-by-theta array of posteriors


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

}


#' @export
marginal_pdf <- function(m, n, x) { # vectorized over x, not n
    npsp::interp(m@marginal_grid, m@marginal_pdfs, cbind(n, x))$y}

#' @export
posterior_pdf <- function(m, theta, n, x) { # vectorized over theta, not n, x
    npsp::interp(m@posterior_grid, m@posterior_pdfs, cbind(n, x, theta))$y}
