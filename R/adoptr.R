#' Adaptive Optimal Two-Stage Designs
#'
#' The \pkg{adoptr} package provides functionality to explore custom optimal
#' two-stage designs for one or two-arm superiority tests.
#' Currently, only (asymptotically) normal test statistics are supported.
#' \pkg{adoptr} is programmed in an object-oriented way.
#' A description on object-oriented usage of \code{R} can be found
#' \href{http://adv-r.had.co.nz/OO-essentials.html}{here}.
#'
#'
#' @section Quickstart:
#'
#' For a sample workflow and a quick demo of the capabilities, see
#' \href{https://kkmann.github.io/adoptr/articles/adoptr.html}{here}.
#'
#' A variety of examples is presented in the validation package
#' \pkg{adoptrValidation} and can be seen
#' \href{https://kkmann.github.io/adoptrValidation/}{here}.
#'
#'
#' @section Designs:
#'
#' \pkg{adoptr} currently supports \code{\link{TwoStageDesign}},
#' \code{\link{GroupSequentialDesign}}, and \code{\link{OneStageDesign}}.
#'
#'
#' @section Data distributions:
#'
#' Currently, the only implemented data distribution is \code{\link{Normal}}.
#'
#'
#' @section Priors:
#'
#' Both \code{\link{ContinuousPrior}} and \code{\link{PointMassPrior}} are
#' supported for the single parameter of a \code{\link{DataDistribution}}.
#' An example on working with priors is provided
#' \href{https://kkmann.github.io/adoptr/articles/working-with-priors.html}{here}.
#'
#'
#' @section Scores:
#'
#' \pkg{adoptr} provides the score types \code{\link{UnconditionalScore}} and
#' \code{\link{ConditionalScore}}. The conditional scores
#' \code{\link{ConditionalPower}} and \code{\link{ConditionalScore}} are
#' already implemented. Unconditional scores that are expectations of
#' conditional scores can be created via \code{\link{expected}} and are
#' represented by the class \code{\link{IntegralScore}}.
#' For an example how to work with scores, see \href{https://kkmann.github.io/adoptr/articles/score-and-constraints-arithmetic.html}{here}.
#'
#' @import methods
#' @docType package
#' @name adoptr
NULL
