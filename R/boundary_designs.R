#' Boundary designs
#'
#' The optimization method \code{\link{minimize}} is based on the package
#' \code{nloptr}. This requires upper and lower boundaries for optimization.
#' Such boundaries can be computed via \code{lower_boundary_design}
#' respectively \code{upper_boundary_design}.
#' They are implemented by default in \code{\link{minimize}}.
#' Note that \code{\link{minimize}} allows the user to define its own
#' boundary designs, too.
#'
#' @param initial_design The initial design
#' @param n1 bound for the first-stage sample size n1
#' @param n2_pivots bound for the second-stage sample size n2
#' @param c1_buffer shift of the early-stopping boundaries from the initial ones
#' @param c2_buffer shift of the final decision boundary from the initial one
#' @param ... optional arguments
#'
#' The values \code{c1f} and \code{c1e} from the initial design are shifted
#' to \code{c1f - c1_buffer} and \code{c1e - c1_buffer} in
#' \code{get_lower_boundary_design}, respectively, to \cr
#' \code{c1f + c1_buffer} and \code{c1e + c1_buffer} in
#' \code{get_upper_boundary_design}.
#' This is handled analogously with \code{c2_pivots} and \code{c2_buffer}.
#'
#' @examples
#' initial_design <- TwoStageDesign(
#'   n1    = 25,
#'   c1f   = 0,
#'   c1e   = 2.5,
#'   n2    = 50,
#'   c2    = 1.96,
#'   order = 7L
#'   )
#' get_lower_boundary_design(initial_design)
#'
#' @rdname boundary-designs
#' @export
setGeneric("get_lower_boundary_design",
           function(initial_design, ...) standardGeneric("get_lower_boundary_design"))



#' @rdname boundary-designs
#' @export
setGeneric("get_upper_boundary_design",
           function(initial_design, ...) standardGeneric("get_upper_boundary_design"))



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("OneStageDesign"),
          function(initial_design, n1 = 1, c1_buffer = 2, ...) {
              lb_design <- OneStageDesign(n1, min(0, initial_design@c1f - c1_buffer))
              lb_design@tunable <- initial_design@tunable
              return(lb_design)

})


#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("GroupSequentialDesign"),
          function(
              initial_design,
              n1        = 1,
              n2_pivots = 1,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              lb_design <- GroupSequentialDesign(
                  n1,
                  initial_design@c1f - c1_buffer,
                  initial_design@c1e - c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots - c2_buffer,
                  order = length(initial_design@c2_pivots)
              )
              lb_design@tunable <- initial_design@tunable
              return(lb_design)

})



#' @rdname boundary-designs
#' @export
setMethod("get_lower_boundary_design", signature("TwoStageDesign"),
          function(
              initial_design,
              n1        = 1,
              n2_pivots = 1,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              lb_design <- TwoStageDesign(
                  n1,
                  initial_design@c1f - c1_buffer,
                  initial_design@c1e - c1_buffer,
                  rep(n2_pivots, length(initial_design@c2_pivots)),
                  initial_design@c2_pivots - c2_buffer,
                  order = length(initial_design@c2_pivots)
                  )
              lb_design@tunable <- initial_design@tunable
              return(lb_design)

})





#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("OneStageDesign"),
          function(initial_design, n1 = 5 * initial_design@n1, c1_buffer = 2, ...) {
              ub_design <- OneStageDesign(n1, initial_design@c1f + c1_buffer)
              ub_design@tunable <- initial_design@tunable
              return(ub_design)
})



#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("GroupSequentialDesign"),
          function(
              initial_design,
              n1        = 5 * initial_design@n1,
              n2_pivots = 5 * initial_design@n2_pivots,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              ub_design <- GroupSequentialDesign(
                  n1,
                  initial_design@c1f + c1_buffer,
                  initial_design@c1e + c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots + c2_buffer,
                  order = length(initial_design@c2_pivots)
            )
            ub_design@tunable <- initial_design@tunable
            return(ub_design)
})


#' @rdname boundary-designs
#' @export
setMethod("get_upper_boundary_design", signature("TwoStageDesign"),
          function(
              initial_design,
              n1        = 5 * initial_design@n1,
              n2_pivots = 5 * initial_design@n2_pivots,
              c1_buffer = 2,
              c2_buffer = 2,
              ...
          ) {
              ub_design <- TwoStageDesign(
                  n1,
                  initial_design@c1f + c1_buffer,
                  initial_design@c1e + c1_buffer,
                  n2_pivots,
                  initial_design@c2_pivots + c2_buffer,
                  order = length(initial_design@c2_pivots)
            )
              ub_design@tunable <- initial_design@tunable
              return(ub_design)

})
